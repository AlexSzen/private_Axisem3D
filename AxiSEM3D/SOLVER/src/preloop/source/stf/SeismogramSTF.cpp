// SeismogramSTF.cpp
// STF produced by seismogram
#include "SeismogramSTF.h"
#include "Processor.h"
#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>      // std::setprecision

SeismogramSTF::SeismogramSTF(const RMatX3 trace, double dt_fwd, double duration_fwd, double hdur_fwd, double decay_fwd,
	const RDMatXX &adjoint_params, const IColX &filter_types):
mHalfDuration(hdur_fwd), mDecay(decay_fwd) {

	mDeltaT = dt_fwd;
    int nStepBeforeZero = ceil(1.5 * mHalfDuration / mDeltaT);
    int nStepAfterZero = ceil(duration_fwd / mDeltaT);
    mShift = nStepBeforeZero * mDeltaT;
    int nStep = nStepBeforeZero + nStepAfterZero;

	if (nStep + 1 != trace.rows()) throw std::runtime_error("SeismogramSTF::SeismogramSTF : length of seismogram and STF inconsistent.");

	mSTFs.assign(nStep + 1, (double) 0.);
	mSTFp.assign(nStep + 1, (double) 0.);
	mSTFz.assign(nStep + 1, (double) 0.);

	// apply each measurement and add to the STF
	for (int i_measurement = 0; i_measurement < adjoint_params.rows(); i_measurement++) {
		RDCol2 window = adjoint_params.block(i_measurement, 0, 1, 2);
		double measurement = adjoint_params(i_measurement, 2);

		RMatX3 trace_measurement_T_disp = trace;
		RMatX3 trace_measurement_T_vel = trace;
		RMatX3 trace_measurement_T_accel = trace;
		CMatX3 trace_measurement_F_vel;
		CMatX3 trace_measurement_F_accel;

		Processor::transformT2F(trace_measurement_T_vel, trace_measurement_F_vel);
		Processor::transformT2F(trace_measurement_T_accel, trace_measurement_F_accel);
		Processor::derivate(trace_measurement_F_vel); //disable for misfit test
		Processor::derivate(trace_measurement_F_accel); //disable for misfit test
		Processor::derivate(trace_measurement_F_accel); //disable for misfit test
		Processor::filter(trace_measurement_F_vel, filter_types[i_measurement]);
		Processor::filter(trace_measurement_F_accel, filter_types[i_measurement]);
		Processor::transformF2T(trace_measurement_F_vel, trace_measurement_T_vel);
		Processor::transformF2T(trace_measurement_F_accel, trace_measurement_T_accel);
		Processor::taper(trace_measurement_T_vel, (Real) window(0), (Real) window(1));
		Processor::taper(trace_measurement_T_accel, (Real) window(0), (Real) window(1));
		Processor::taper(trace_measurement_T_disp, (Real) window(0), (Real) window(1));

		Real norm_s = 0.;
		Real norm_p = 0.;
		Real norm_z = 0.;
		Real max = 0.;
		std::ofstream myfile;
		myfile.open("/home/alex/Desktop/phd/axisem3d_alex/AxiSEM3D/build/output/stations/test.txt");
		for (int it = 0; it <= nStep; it++) { //normalization factor : for traveltime tomo it's time integrated squared velocity
											 // for amplitude it's displacement
			norm_s += mDeltaT * trace_measurement_T_vel(it, 0) * trace_measurement_T_vel(it, 0);
			norm_p += mDeltaT * trace_measurement_T_vel(it, 1) * trace_measurement_T_vel(it, 1);
			norm_z += mDeltaT * trace_measurement_T_vel(it, 2) * trace_measurement_T_vel(it, 2);
		}
		
	//	std::cout<<std::setprecision(15)<<norm_s/(Real)2.<<std::endl; //this is value of misfit in misfit test
		if (norm_s == 0.) norm_s = 1.;
		if (norm_p == 0.) norm_p = 1.;
		if (norm_z == 0.) norm_z = 1.;

		for (int i = 0; i <= nStep; i++) { //time reversed seismogram.
				
	   	mSTFs[i] = trace_measurement_T_vel(nStep - i, 0) / norm_s;
	   	// mSTFp[i] = measurement * trace_measurement_T(nStep - i, 1) / norm_p;
	   	// mSTFz[i] = measurement * trace_measurement_T(nStep - i, 2) / norm_z;
		
		//mSTFs[i] = trace_measurement_T(nStep - i, 0); // stf for misfit test
		myfile<<mSTFs[i]<<"\n";

	    }
		myfile.close();

	}


}

std::string SeismogramSTF::verbose() const {
	std::stringstream ss;
	ss << "\n=================== Source Time Function ===================" << std::endl;
	ss << "  Time Step               =   " << mDeltaT << std::endl;
	ss << "  Number of Steps         =   " << mSTFs.size() << std::endl;
	ss << "  Total Duration          =   " << mDeltaT * mSTFs.size() << std::endl;
	ss << "  Duration after Origin   =   " << mDeltaT * mSTFs.size() - mShift << std::endl;
	ss << "  Shift before Origin     =   " << mShift << std::endl;
	ss << "  Time Series Type        =   Seismogram" << std::endl;
	ss << "  Half Duration           =   N/A"  << std::endl;
	ss << "  Decay Factor            =   N/A" << std::endl;
	ss << "=================== Source Time Function ===================\n" << std::endl;
	return ss.str();

}
