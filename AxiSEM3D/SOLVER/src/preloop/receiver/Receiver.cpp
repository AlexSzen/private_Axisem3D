// Receiver.cpp
// created by Kuangdai on 13-May-2016
// receiver

#include "Receiver.h"
#include "XMath.h"
#include "Geodesy.h"
#include "Domain.h"
#include "Mesh.h"
#include "Quad.h"
#include "Element.h"
#include "SpectralConstants.h"
#include "PointwiseRecorder.h"
#include "Parameters.h"

#include <sstream>
#include <iomanip>
#include "XMPI.h"

Receiver::Receiver(const std::string &name, const std::string &network,
    double theta_lat, double phi_lon, bool geographic,
    double depth, double srcLat, double srcLon, double srcDep):
mName(name), mNetwork(network), mDepth(depth) {
    RDCol3 rtpG, rtpS;
    if (geographic) {
        rtpG(0) = 1.;
        rtpG(1) = Geodesy::lat2Theta_d(theta_lat, mDepth);
        rtpG(2) = Geodesy::lon2Phi(phi_lon);
        rtpS = Geodesy::rotateGlob2Src(rtpG, srcLat, srcLon, srcDep);

        // For off axis point force tests, during normal (axial) simulation,
        // get coordinates of equivalent receiver for the off axis simulation with
        // off axis source at given source centered theta and phi.

        // convert to cartesian coords
        RDCol3 xyzG = Geodesy::toCartesian(rtpG);

        // theta and phi of off axis source (here on the equator), in source-centered coords.
        double off_axis_theta = pi/2;
        double off_axis_phi = 0.;

        // rotation matrix
        RDMat33 rot_matrix = Geodesy::rotationMatrix(off_axis_theta, off_axis_phi);

        // rotate cartesian coords
        RDCol3 off_axis_xyzG = rot_matrix * xyzG;

        // convert to spherical
        bool defined;
        RDCol3 off_axis_rtpG = Geodesy::toSpherical(off_axis_xyzG, defined);

        // convert to geographic
        double off_axis_theta_lat = Geodesy::theta2Lat_d(off_axis_rtpG(1), mDepth);
        double off_axis_phi_lon = Geodesy::phi2Lon(off_axis_rtpG(2));
        double off_axis_depth = mDepth;
/*
        if (XMPI::rank()==0) {
        //write rotated receivers for point force test
            std::ofstream myfile(Parameters::sInputDirectory + "/rotated_stations.txt", std::ios_base::app);
        
            if (myfile.is_open()) {
                myfile << mName <<"\t" << mNetwork << "\t" << off_axis_theta_lat
                << "\t" << off_axis_phi_lon << "\t" << 0. << "\t" <<mDepth<<"\n";
                myfile.close();
            }
        }
        XMPI::cout<<"New receiver "<< mName <<" is latitude: "<< off_axis_theta_lat<< " longitude "<<off_axis_phi_lon<< " and depth "<<mDepth<<XMPI::endl;

        //TODO:save this to file STATIONS_ADJ_PFTEST
*/

    } else {
        rtpS(0) = 1.;
        rtpS(1) = theta_lat * degree;
        rtpS(2) = phi_lon * degree;
        rtpG = Geodesy::rotateSrc2Glob(rtpS, srcLat, srcLon, srcDep);
    }
    mTheta = rtpS(1);
    mPhi = rtpS(2);
    mLat = Geodesy::theta2Lat_d(rtpG(1), mDepth);
    mLon = Geodesy::phi2Lon(rtpG(2));
    mBackAzimuth = Geodesy::backAzimuth(srcLat, srcLon, srcDep, mLat, mLon, mDepth);
    // // test
    // XMPI::cout << name << " " << network << " ";
    // XMPI::cout << mLat << " " << mLon << " " << " 0.0 " << mDepth << XMPI::endl;
}

void Receiver::release(PointwiseRecorder &recorderPW, const Domain &domain,
    int elemTag, const RDMatPP &interpFact) {
    Element *myElem = domain.getElement(elemTag);
    recorderPW.addReceiver(mName, mNetwork, mPhi, interpFact, myElem, mTheta, mBackAzimuth,
        mLat, mLon, mDepth);
}

bool Receiver::locate(const Mesh &mesh, int &elemTag, RDMatPP &interpFact) const {
    RDCol2 recCrds, srcXiEta;
    double r = mesh.computeRadiusRef(mDepth, mLat, mLon);
    recCrds(0) = r * sin(mTheta);
    recCrds(1) = r * cos(mTheta);
    if (recCrds(0) > mesh.sMax() + tinySingle || recCrds(0) < mesh.sMin() - tinySingle) {
        return false;
    }
    if (recCrds(1) > mesh.zMax() + tinySingle || recCrds(1) < mesh.zMin() - tinySingle) {
        return false;
    }
    for (int iloc = 0; iloc < mesh.getNumQuads(); iloc++) {
        const Quad *quad = mesh.getQuad(iloc);
        if (!quad->nearMe(recCrds(0), recCrds(1))) {
            continue;
        }
        if (quad->invMapping(recCrds, srcXiEta)) {
            if (std::abs(srcXiEta(0)) <= 1.000001 && std::abs(srcXiEta(1)) <= 1.000001) {
                elemTag = quad->getElementTag();
                RDColP interpXi, interpEta;
                XMath::interpLagrange(srcXiEta(0), nPntEdge,
                    quad->isAxial() ? SpectralConstants::getP_GLJ().data():
                    SpectralConstants::getP_GLL().data(), interpXi.data());
                XMath::interpLagrange(srcXiEta(1), nPntEdge,
                    SpectralConstants::getP_GLL().data(), interpEta.data());
                interpFact = interpXi * interpEta.transpose();
                return true;
            }
        }
    }
    return false;
}

std::string Receiver::verbose(bool geographic, int wname, int wnet) const {
    std::stringstream ss;
    ss << std::setw(wname) << mName << "   ";
    ss << std::setw(wnet) << mNetwork << "   ";
    if (geographic) {
        ss << std::setw(8) << mLat << "   ";
        ss << std::setw(8) << mLon << "   ";
    } else {
        ss << std::setw(8) << mTheta / degree << "   ";
        ss << std::setw(8) << mPhi / degree << "   ";
    }
    ss << mDepth;
    return ss.str();
}
