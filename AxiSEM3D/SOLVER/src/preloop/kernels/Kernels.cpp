// kernels.cpp
//build params for kernel computation 

#include "Kernels.h"
#include "Domain.h"
#include "Mesh.h"
#include "Quad.h"
#include "Element.h"
#include "Parameters.h"
#include "XMPI.h"
#include <sstream>
#include <boost/algorithm/string.hpp>
#include "Kerner.h"
#include "KernerElement.h"

Kernels::Kernels(bool computeKer): mComputeKernels(computeKer) {
	
}


Kernels::~Kernels() {
	
}

void Kernels::buildInparam(Kernels *&kernels, const Parameters &par, int totalStepsSTF,  int verbose) {
	
	bool computeKer = par.getValue<bool>("COMPUTE_KERNELS");

	kernels = new Kernels(computeKer);
	if (!computeKer) return;
	
	kernels->mMaxStep = totalStepsSTF;
	kernels->mDumpTimeKernels = par.getValue<bool>("OUTPUT_TIME_KERNELS");
	
	int recInterval = par.getValue<int>("WAVEFIELD_RECORD_INTERVAL");
	kernels->mRecordInterval = recInterval;
	if (recInterval <= 0) {
		recInterval = 1;
	}
	
	kernels->mTotalStepsKernels = totalStepsSTF / recInterval;
	if (totalStepsSTF % recInterval > 0) {
		kernels->mTotalStepsKernels += 1;
	}
	
	kernels->mBufferSize = par.getValue<int>("WAVEFIELD_DUMP_INTERVAL");
	if (kernels->mBufferSize <= 0) {
		kernels->mBufferSize = 1;
	}
	
	if (verbose) XMPI::cout << kernels->verbose();

}

void Kernels::release(Domain &domain, Mesh &mesh, const Parameters &par, double deltaT) {
	
	if (!mComputeKernels) return;
        
        // ---------- <GET PARAMETERS> ----------
        double rmin = par.getValue<double>("RMIN")*1e3;
        double rmax = par.getValue<double>("RMAX")*1e3;
        double tmin = par.getValue<double>("THETA_MIN")*degree;
        double tmax = par.getValue<double>("THETA_MAX")*degree;
        int gll_ani = par.getValue<int>("GLL_ANIMATION");
        // ---------- </GET PARAMETERS> ----------
	
	Kerner *kerner = new Kerner(mDumpTimeKernels, mTotalStepsKernels, mBufferSize, mRecordInterval, mMaxStep); 
     
	for (int ielem = 0; ielem < domain.getNumElements(); ielem ++) {
        if (domain.getElement(ielem)->needDumping(rmin,rmax,tmin,tmax)) {
		    KernerElement *kerElem = new KernerElement(domain.getElement(ielem), mesh.getQuad(ielem)->getIntegralFactor(), mBufferSize, (Real) deltaT * mRecordInterval);
		    kerner->addKernerElement(kerElem);
        }
	}
	kerner->setDomainRecorder(domain.getDomainRecorder());
	domain.setKerner(kerner);
}

std::string Kernels::verbose() {
	
	std::stringstream ss;
	ss << "\n========================= Kernels ========================" << std::endl;
	ss << "Kernels will be computed" <<std::endl; 
	
	ss << "========================= Kernels ========================\n" << std::endl;
	return ss.str();
}


