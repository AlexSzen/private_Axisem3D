// axisem.cpp
// created by Kuangdai on 26-Mar-2016 
// main of AxiSEM3D

#include "axisem.h"
#include "MultilevelTimer.h"

#include "XMPI.h"
#include "eigenc.h"
#include "eigenp.h"

int axisem_main(int argc, char *argv[]) {
    
    try {
        
        // variable sets
        PreloopVariables pl;
        SolverVariables sv;
        
        // initialize mpi
        XMPI::initialize(argc, argv);
        
        //////// spectral-element constants
        SpectralConstants::initialize(nPol);  
        
        //////// input parameters 
        int verbose;
        Parameters::buildInparam(pl.mParameters, verbose);
        
        //////// preloop timer
        MultilevelTimer::initialize(Parameters::sOutputDirectory + "/develop/preloop_timer.txt", 4);
        if (pl.mParameters->getValue<bool>("DEVELOP_DIAGNOSE_PRELOOP")) {
            MultilevelTimer::enable();
        }
        
        //////// exodus model and attenuation parameters 
        MultilevelTimer::begin("Build Exodus", 0);
        ExodusModel::buildInparam(pl.mExodusModel, *(pl.mParameters), pl.mAttParameters, verbose);
        MultilevelTimer::end("Build Exodus", 0);
        
        //////// fourier field 
        MultilevelTimer::begin("Build NrField", 0);
        NrField::buildInparam(pl.mNrField, *(pl.mParameters), verbose);
        MultilevelTimer::end("Build NrField", 0);
		
		MultilevelTimer::begin("Build Sources", 0);
		SourceCollection::buildInparam(pl.mSources, *(pl.mParameters), verbose);
		double srcLat = pl.mSources->getLatitude();
		double srcLon = pl.mSources->getLongitude();
		double srcDep = pl.mSources->getDepth();
		MultilevelTimer::end("Build Sources", 0);
        //////// 3D models 
        MultilevelTimer::begin("Build 3D Models", 0);
        Volumetric3D::buildInparam(pl.mVolumetric3D, *(pl.mParameters), pl.mExodusModel, 
            srcLat, srcLon, srcDep, verbose);
        Geometric3D::buildInparam(pl.mGeometric3D, *(pl.mParameters), verbose);
        OceanLoad3D::buildInparam(pl.mOceanLoad3D, *(pl.mParameters), verbose);
        MultilevelTimer::end("Build 3D Models", 0);
		
        
        //////// mesh, phase 1
        // define mesh
        MultilevelTimer::begin("Mesh Definition", 0);
        pl.mMesh = new Mesh(pl.mExodusModel, pl.mNrField, 
            srcLat, srcLon, srcDep, *(pl.mParameters), verbose);
        pl.mMesh->setVolumetric3D(pl.mVolumetric3D);
        pl.mMesh->setGeometric3D(pl.mGeometric3D);
        pl.mMesh->setOceanLoad3D(pl.mOceanLoad3D);
        MultilevelTimer::end("Mesh Definition", 0);
        
        // build unweighted local mesh 
        MultilevelTimer::begin("Build Unweighted Mesh", 0);
        pl.mMesh->buildUnweighted();
        MultilevelTimer::end("Build Unweighted Mesh", 0);
        
        //////// static variables in solver, mainly FFTW
        bool disableWisdomFFTW = pl.mParameters->getValue<bool>("DEVELOP_DISABLE_FFTW_WISDOM");
        MultilevelTimer::begin("Initialize FFTW", 0);
        initializeSolverStatic(pl.mMesh->getMaxNr(), disableWisdomFFTW); 
        MultilevelTimer::end("Initialize FFTW", 0);
        //////// dt
        MultilevelTimer::begin("Compute DT", 0);
        double dt = pl.mParameters->getValue<double>("TIME_DELTA_T");
        if (dt < tinyDouble) {
            dt = pl.mMesh->getDeltaT();
        }
        double dt_fact = pl.mParameters->getValue<double>("TIME_DELTA_T_FACTOR");
        if (dt_fact < tinyDouble) {
            dt_fact = 1.0;
        }
        dt *= dt_fact;
        MultilevelTimer::end("Compute DT", 0);
        
        //////// attenuation
        MultilevelTimer::begin("Build Attenuation", 0);
        AttBuilder::buildInparam(pl.mAttBuilder, *(pl.mParameters), pl.mAttParameters, dt, verbose);
        MultilevelTimer::end("Build Attenuation", 0);
        
        //////// mesh, phase 2
        MultilevelTimer::begin("Build Weighted Mesh", 0);
        pl.mMesh->setAttBuilder(pl.mAttBuilder);
        pl.mMesh->buildWeighted();
        MultilevelTimer::end("Build Weighted Mesh", 0);
        
        //////// mesh test 
        // test positive-definiteness and self-adjointness of stiffness and mass matrices
        // better to turn with USE_DOUBLE 
        // pl.mMesh->test();
        // XMPI::barrier();
        // exit(0);
		
		//////// source time functions
		MultilevelTimer::begin("Build Source Time Functions", 0);
		STFCollection::buildInparam(pl.mSTFs, *(pl.mParameters), dt, verbose);
		MultilevelTimer::end("Build Source Time Functions", 0);

        //////// receivers
        MultilevelTimer::begin("Build Receivers", 0);
        ReceiverCollection::buildInparam(pl.mReceivers, *(pl.mParameters), 
            srcLat, srcLon, srcDep, pl.mSTFs->getSize(), verbose);
        MultilevelTimer::end("Build Receivers", 0);    
		
		//////// Kernels 
		MultilevelTimer::begin("Build Kernels", 0);
		Kernels::buildInparam(pl.mKernels, *(pl.mParameters), pl.mSTFs->getSize(), verbose);
		MultilevelTimer::end("Build Kernels", 0);

        //////// computational domain
        MultilevelTimer::begin("Computational Domain", 0);
        sv.mDomain = new Domain();
        
        // release mesh
        MultilevelTimer::begin("Release Mesh", 1);
        pl.mMesh->release(*(sv.mDomain));
        MultilevelTimer::end("Release Mesh", 1);

		MultilevelTimer::begin("Release Sources", 1);
		pl.mSources->release(*(sv.mDomain), *(pl.mMesh));
		MultilevelTimer::end("Release Sources", 1);
		
		// release stfs 
		MultilevelTimer::begin("Release STFs", 1);
		pl.mSTFs->release(*(sv.mDomain));
		MultilevelTimer::end("Release STFs", 1);      

		// dump mesh quantities 
		MultilevelTimer::begin("Dump mesh quantities",1);
		pl.mMesh->dumpFields(*(sv.mDomain), *(pl.mParameters));
		MultilevelTimer::end("Dump mesh quantities",1);
		//throw std::runtime_error("Need to stop now bro");
        // release receivers
        MultilevelTimer::begin("Release Receivers", 1);
        pl.mReceivers->release(*(sv.mDomain), *(pl.mMesh));
        sv.mDomain->initializeRecorders();
        MultilevelTimer::end("Release Receivers", 1);
		
		//release kernels
		MultilevelTimer::begin("Release Kerner", 1);
		pl.mKernels->release(*(sv.mDomain), *(pl.mMesh), *(pl.mParameters), pl.mSTFs->getDeltaT());
		sv.mDomain->initializeKerner();
		MultilevelTimer::end("Release Kerner", 1);
        
        // verbose domain 
        MultilevelTimer::begin("Verbose", 1);
        if (verbose) {
            XMPI::cout << sv.mDomain->verbose();
        }
        MultilevelTimer::end("Verbose", 1);
        MultilevelTimer::end("Computational Domain", 0);
        
        MultilevelTimer::finalize();
        
        //////////////////////// PREPROCESS DONE ////////////////////////
        
        //////// Newmark
        int infoInt = pl.mParameters->getValue<int>("OPTION_LOOP_INFO_INTERVAL");
        int stabInt = pl.mParameters->getValue<int>("OPTION_STABILITY_INTERVAL");
        sv.mNewmark = new Newmark(sv.mDomain, infoInt, stabInt);
        
        //////// final preparations
        // finalize preloop variables before time loop starts
        pl.finalize();
        // forbid matrix allocation in time loop
        #ifndef NDEBUG
            Eigen::internal::set_is_malloc_allowed(false);
        #endif
            
        //////// GoGoGo
        XMPI::barrier();
        sv.mNewmark->solve(verbose);

        //////// finalize solver
        // solver 
        sv.mDomain->finalizeKerner();
		sv.mDomain->finalizeRecorders();
        sv.finalize();
        // static variables in solver
        finalizeSolverStatic();
        
        // finalize mpi 
        XMPI::finalize();
        
    } catch (const std::exception &e) {
        // print exception
        XMPI::cout.setp(XMPI::rank());
        XMPI::printException(e);
        
        // abort program
        // TODO 
        // MPI_Abort is necessary here. Otherwise, if an exception
        // is thrown from one of the procs, deadlock will occur.
        // But the problem is, how we free memories on other procs?!
        XMPI::abort();
    }
    
    return 0;
}

#include "SolverFFTW.h"
#include "SolverFFTW_1.h"
#include "SolverFFTW_3.h"
#include "SolverFFTW_N3.h"
#include "SolverFFTW_N6.h"
#include "SolverFFTW_N9.h"
#include "PreloopFFTW.h"
#include "SolidElement.h"
#include "FluidElement.h"

extern void initializeSolverStatic(int maxNr, bool disableWisdomFFTW) {
    // fftw
    SolverFFTW::importWisdom(disableWisdomFFTW);
    SolverFFTW_1::initialize(maxNr);
    SolverFFTW_3::initialize(maxNr); 
    SolverFFTW_N3::initialize(maxNr);
    SolverFFTW_N6::initialize(maxNr);
    SolverFFTW_N9::initialize(maxNr);
    SolverFFTW::exportWisdom();
    // PreloopFFTW::initialize(maxNr);
    // element
    SolidElement::initWorkspace(maxNr / 2);
    FluidElement::initWorkspace(maxNr / 2);
};

extern void finalizeSolverStatic() {
    // fftw
    SolverFFTW_1::finalize();
    SolverFFTW_3::finalize(); 
    SolverFFTW_N3::finalize();
    SolverFFTW_N6::finalize();
    SolverFFTW_N9::finalize();
    PreloopFFTW::finalize();
};
