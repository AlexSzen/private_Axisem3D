//KernerIO.cpp

#include "KernerIO.h"
#include "NetCDF_Writer.h"
#include "NetCDF_Reader.h"
#include "Parameters.h"
#include "XMPI.h"
#include <iostream>
#include <numeric>

KernerIO::KernerIO(bool dumpTimeKernels, int temp_startElem, int temp_countElem, int totSteps, int bufferSize):
mDumpTimeKernels(dumpTimeKernels), mStartElem(temp_startElem), mCountElem(temp_countElem), mTotSteps(totSteps), mBufferSize(bufferSize) 
 {
	
	mNetCDF_w = new NetCDF_Writer();
	mNetCDF_r = new NetCDF_Reader();
	mNumLoads = 0;
	
	
	// open forward 
	std::string fname_wvf;
	
	#ifdef _USE_PARALLEL_NETCDF
		fname_wvf = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd.nc4";
		mNetCDF_r->openParallel(fname_wvf);
	#else 
		fname_wvf = Parameters::sOutputDirectory + "/wavefields/wavefield_db_fwd_" + std::to_string(XMPI::rank()) + ".nc4";
		mNetCDF_r->open(fname_wvf);	
	#endif
}

void KernerIO::initialize(int totNuProc, int totNu, int startElemNu, int countElemNu, int totElemProc, int totElem, const std::vector<int> &nusKer,const std::vector<int> &nrsKer) {

	std::string fname_ker ;
	
	#ifdef _USE_PARALLEL_NETCDF
		// init kernel variables
		fname_ker = Parameters::sOutputDirectory + "/kernels/kernels_db.nc4";
		mStartElemNuKernels = startElemNu;
		mCountElemNuKernels = countElemNu;
		
		if (XMPI::root()) {
		
				//create dims and vars in kernel file 
				mNetCDF_w->open(fname_ker, true);
				
				std::vector<size_t> dims_kernels; 
				std::vector<size_t> dims_nus;
				
				if (mDumpTimeKernels) {
					dims_kernels.push_back( mTotSteps );
				} else {
					dims_kernels.push_back( 1 );
				}		
				dims_kernels.push_back( totNu );
				dims_kernels.push_back( 6 ); //real and imag parts of elastic radial ani kernels. 
				dims_kernels.push_back( nPntEdge );
				dims_kernels.push_back( nPntEdge );
				
				dims_nus.push_back(totElem);
				
				mNetCDF_w->defModeOn();
				mNetCDF_w->defineVariable<Real>("Kernels", dims_kernels);
				mNetCDF_w->defineVariable<int>("Nus", dims_nus);
				mNetCDF_w->defineVariable<int>("Nrs", dims_nus);

				
				mNetCDF_w->defModeOff();
				mNetCDF_w->close();
				
			}
			
			mNetCDF_w->openParallel(fname_ker);

			std::vector<size_t> startNus, countNus;
			startNus.push_back(mStartElem);
			countNus.push_back(mCountElem);	
			mNetCDF_w->writeVariableChunk("Nus", nusKer, startNus, countNus);
			mNetCDF_w->writeVariableChunk("Nrs", nrsKer, startNus, countNus);
		
	#else //serial netcdf
		fname_ker = Parameters::sOutputDirectory + "/kernels/kernels_db_" + std::to_string(XMPI::rank()) + ".nc4";
		mStartElemNuKernels = 0;
		mCountElemNuKernels = totNuProc;
		
		//create dims and vars in kernel file 
		mNetCDF_w->open(fname_ker, true);
		
		std::vector<size_t> dims_kernels; 
		std::vector<size_t> dims_nus;
		
		if (mDumpTimeKernels) {
			dims_kernels.push_back( mTotSteps );
		} else {
			dims_kernels.push_back( 1 );
		}		
		dims_kernels.push_back( totNuProc );
		dims_kernels.push_back( 6 ); //real and imag parts of elastic radial ani kernels. 
		dims_kernels.push_back( nPntEdge );
		dims_kernels.push_back( nPntEdge );
		
		dims_nus.push_back(totElemProc);
		
		mNetCDF_w->defModeOn();
		mNetCDF_w->defineVariable<Real>("Kernels", dims_kernels);
		mNetCDF_w->defineVariable<int>("Nus", dims_nus);
		mNetCDF_w->defineVariable<int>("Nrs", dims_nus);

		
		mNetCDF_w->defModeOff();
		
		mNetCDF_w->writeVariableWhole("Nus", nusKer);
		mNetCDF_w->writeVariableWhole("Nrs", nrsKer);
		
	#endif
	


}

void KernerIO::finalize() {
	
	mNetCDF_r->close();
	mNetCDF_w->close();
	delete mNetCDF_w;
	delete mNetCDF_r;
	
}

void KernerIO::dumpToFile(const vec_vec_ar6_RMatPP &kernels, int bufferSize) {
	
	std::vector<size_t> startKernels, countKernels;
	startKernels.push_back(mTimeLine);
	startKernels.push_back(mStartElemNuKernels);
	startKernels.push_back(0);
	startKernels.push_back(0);
	startKernels.push_back(0);
	
	countKernels.push_back(1); //have to write one by one because of netcdf issue
	countKernels.push_back(mCountElemNuKernels);
	countKernels.push_back(6);
	countKernels.push_back(nPntEdge);
	countKernels.push_back(nPntEdge); 
	
	int dumpSize = 1;
	if (mDumpTimeKernels) {
		dumpSize = bufferSize;
	}
	for (int it = 0; it < dumpSize ; it ++) {
		
		mNetCDF_w->writeVariableChunk("Kernels", kernels[it], startKernels, countKernels);
		startKernels[0]++;
		
	}
	mTimeLine += dumpSize;

}

void KernerIO::loadNus(std::vector<int> &Nus) {
	
	//read nus 
	std::vector<size_t> startElem, countElem;
	startElem.push_back(mStartElem);
	countElem.push_back(mCountElem);
	Nus.assign(mCountElem, 0);
	
	#ifdef _USE_PARALLEL_NETCDF
		mNetCDF_r->readVariableChunk("Nus", Nus, startElem, countElem);
	#else 
		mNetCDF_r->read1D("Nus", Nus);
	#endif
		
}

void KernerIO::loadNrs(std::vector<int> &Nrs) {
	
	//read nrs 
	std::vector<size_t> startElem, countElem;
	startElem.push_back(mStartElem);
	countElem.push_back(mCountElem);
	Nrs.assign(mCountElem, 0);
	
	#ifdef _USE_PARALLEL_NETCDF
		mNetCDF_r->readVariableChunk("Nrs", Nrs, startElem, countElem);
	#else 
		mNetCDF_r->read1D("Nrs", Nrs);
	#endif

		
}

void KernerIO::loadWavefield(vec_vec_ar6_RMatPP &disp, int bufferSize) {
	
	std::vector<size_t> startElemNu, countElemNu;
	
	int startChunk = mTotSteps - bufferSize * (1 + mNumLoads);
	if (bufferSize != mBufferSize) startChunk = 0; //at last step we change bufferSize, i.e. read last chunk

	startElemNu.push_back(startChunk); //load starting from end chunks
	startElemNu.push_back(mStartElemNuFwd);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	startElemNu.push_back(0);

	countElemNu.push_back(1); // have to load one by one because of netcdf unresolved issue 
	countElemNu.push_back(mCountElemNuFwd);
	countElemNu.push_back(6);
	countElemNu.push_back(nPntEdge);
	countElemNu.push_back(nPntEdge);
	
	// fill disp with 0 
	vec_ar6_RMatPP initBuf(mCountElemNuFwd, zero_ar6_RMatPP);
	for (int it = 0; it < bufferSize; it++) {//have to read one by one 
		mNetCDF_r->readVariableChunk("displacement_wavefield", initBuf, startElemNu, countElemNu);
		disp.push_back(initBuf);
		startElemNu[0]++;
	}
	
	mNumLoads ++;
	
	
	
	
			
}

void KernerIO::loadMaterial(vec_ar12_RMatPP &materials) {

	//read nus 
	std::vector<size_t> startElem, countElem;
	startElem.push_back(mStartElem);
	countElem.push_back(mCountElem);
	std::vector<int> Nus(mCountElem, 0);
	std::vector<int> Nrs(mCountElem, 0);
	
	#ifdef _USE_PARALLEL_NETCDF
		mNetCDF_r->readVariableChunk("Nus", Nus, startElem, countElem);
	#else 
		mNetCDF_r->read1D("Nus", Nus);
	#endif

	// create start and count for elemNu
	int totNuProc = 0;
	for (int i = 0; i<Nus.size(); i++) totNuProc+=Nus[i];
	int temp_startElemNu;
	std::vector<int> temp_countElemNu(XMPI::nproc(),0);
	XMPI::gather(totNuProc, temp_countElemNu, true);
	temp_startElemNu = std::accumulate(temp_countElemNu.begin(), temp_countElemNu.begin() + XMPI::rank(), 0);
	
	#ifdef _USE_PARALLEL_NETCDF
		mCountElemNuFwd = temp_countElemNu[XMPI::rank()];
		mStartElemNuFwd = temp_startElemNu;
	#else 
		mCountElemNuFwd = totNuProc;
		mStartElemNuFwd = 0;
	#endif
	std::vector<size_t> startElemNu, countElemNu;
	
	startElemNu.push_back(mStartElemNuFwd);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	startElemNu.push_back(0);
	
	countElemNu.push_back(mCountElemNuFwd);
	countElemNu.push_back(12);
	countElemNu.push_back(nPntEdge);
	countElemNu.push_back(nPntEdge);

	// fill with 0
	materials.assign(totNuProc, zero_ar12_RMatPP);
	
	mNetCDF_r->readVariableChunk("material_fields", materials, startElemNu, countElemNu);
		
}