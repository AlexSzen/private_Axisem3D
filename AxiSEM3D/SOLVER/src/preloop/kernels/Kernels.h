// builds params for kerner and signal proc 

#pragma once 

#include <vector>
#include <string>
#include "global.h"
#include "eigenc.h"

class Parameters;
class Domain;
class Mesh;

class Kernels {

public:
	
	Kernels(bool computeKer); // if computeKer = 0 we call this constructor 
	~Kernels();
	
	std::string verbose();
	
	static void buildInparam(Kernels *&kernels, const Parameters &par, int totalStepsSTF, int verbose);
	void release(Domain &domain, Mesh &mesh, const Parameters &par, double deltaT); //actually need params for partial dumping
	
private:
	
	bool mComputeKernels;
	bool mDumpTimeKernels;

	double mRmin = 0., mRmax = 7.e6;
	double mThetaMin = 0., mThetaMax = 180.;
	
	int mTotalStepsKernels, mBufferSize, mRecordInterval;
	int mMaxStep; // max step for time loop
	
	
};