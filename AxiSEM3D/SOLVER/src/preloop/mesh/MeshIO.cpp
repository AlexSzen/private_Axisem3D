#include "MeshIO.h"
#include "Parameters.h"
#include "NetCDF_Writer.h"
#include "NetCDF_Reader.h"
#include "eigenc.h"
#include "Domain.h"
#include "Element.h"
#include "Mesh.h"
#include "Point.h"
#include "XMPI.h"
#include "XMath.h"
#include "SlicePlot.h"
#include "Quad.h"
#include <iostream>
#include <numeric>

MeshIO::MeshIO(const Mesh *mesh, const std::string &fname): mMesh(mesh), mFileName(fname) {
	//nothing
}


void MeshIO::dumpFields(const Domain &domain, const Parameters &par) {
	
	// ---------- <GET PARAMETERS> ----------
	double rmin = par.getValue<double>("RMIN")*1e3;
	double rmax = par.getValue<double>("RMAX")*1e3;
	double tmin = par.getValue<double>("THETA_MIN")*degree;
	double tmax = par.getValue<double>("THETA_MAX")*degree;
	int gll_ani = par.getValue<int>("GLL_ANIMATION");
	// ---------- </GET PARAMETERS> ----------

	// ---------- <GENERAL> ----------
	
	
	int gllpoints_proc = domain.getNumPoints();
	int gllpoints_all = XMPI::sum(gllpoints_proc);

	int elems_proc = domain.getNumElements();
	int elems_all = XMPI::sum(elems_proc);

	int elems_proc_ani = 0, elemNus_proc = 0, elemNus_proc_ani = 0;

	for (int i = 0; i < elems_proc; i++) {
		int maxNu = mMesh->getQuad(i)->getNu() + 1;
		elemNus_proc += maxNu;
	//	if (domain.getElement(i)->needDumping(rmin,rmax,tmin,tmax)) {
			elems_proc_ani++;
			elemNus_proc_ani += maxNu;	
	//	}
	}
	
	
	int elems_all_ani = XMPI::sum(elems_proc_ani);
	int elemNus_all = XMPI::sum(elemNus_proc);
	int elemNus_all_ani = XMPI::sum(elemNus_proc_ani);
	// ---------- </GENERAL> ----------
		
	// ---------- <GET FIELDS> ----------	
	std::vector<int> NusAni(elems_proc,0); 
	std::vector<int> NrsAni(elems_proc,0);
	std::vector<int> Nus(elems_proc,0); 
	std::vector<int> Nrs(elems_proc,0);
	vec_RMatPP mesh_S(elems_proc, RMatPP::Zero());
	vec_RMatPP mesh_Z(elems_proc, RMatPP::Zero());
	vec_RMatPP integral_factor(elems_proc, RMatPP::Zero());
	vec_ar2_RMatPP vp(elemNus_proc, zero_ar2_RMatPP);
	vec_ar2_RMatPP vp1D(elemNus_proc, zero_ar2_RMatPP);
	vec_ar2_RMatPP rho(elemNus_proc, zero_ar2_RMatPP);
	vec_ar2_RMatPP rho1D(elemNus_proc, zero_ar2_RMatPP);

	RMatXX_RM zero_matAni(gll_ani,gll_ani); //animations 
	zero_matAni.setZero();	
	vec_RMatXX_RM zero_ar2_matAni(2,zero_matAni);

	vec_ar12_RMatPP materials;
	int nuline = 0;
	for (int i=0; i<elems_proc; i++) {
		
		Element* elem = domain.getElement(i);
		Quad* quad = mMesh->mQuads[i];
		 
	//	if (elem->needDumping(rmin,rmax,tmin,tmax)) { //for now keep the flexible dumping even for wavefields for kernels
      	NusAni[i] = elem->getMaxNu()+1;
      	NrsAni[i] = elem->getMaxNr();
		Nus[i] = elem->getMaxNu()+1;
		Nrs[i] = elem->getMaxNr();
		
		vec_CMatPP vph_elem = quad->getMaterialFourier("vph", SlicePlot::PropertyRefTypes::Property3D);
		vec_CMatPP vpv_elem = quad->getMaterialFourier("vpv", SlicePlot::PropertyRefTypes::Property3D);		
		vec_CMatPP vsh_elem = quad->getMaterialFourier("vsh", SlicePlot::PropertyRefTypes::Property3D);
		vec_CMatPP vsv_elem = quad->getMaterialFourier("vsv", SlicePlot::PropertyRefTypes::Property3D);
		vec_CMatPP rho_elem = quad->getMaterialFourier("rho", SlicePlot::PropertyRefTypes::Property3D);
		vec_CMatPP rho_elem1D = quad->getMaterialFourier("rho", SlicePlot::PropertyRefTypes::Property1D); //for misfit tests 
		vec_CMatPP eta_elem = quad->getMaterialFourier("eta", SlicePlot::PropertyRefTypes::Property3D);
		vec_CMatPP vp_elem = quad->getMaterialFourier("vp", SlicePlot::PropertyRefTypes::Property3D);
		vec_CMatPP vp_elem1D = quad->getMaterialFourier("vp", SlicePlot::PropertyRefTypes::Property1D); //for misfit tests 
		vec_CMatPP vs_elem = quad->getMaterialFourier("vs", SlicePlot::PropertyRefTypes::Property3D);

		for (int ialpha=0; ialpha<Nus[i]; ialpha++) { //i assume materials have this expansion order..? they are padded to it anw for kernel computation 


			materials.push_back(zero_ar12_RMatPP);
			
			if (ialpha<rho_elem.size()) {
				materials.back()[0] = rho_elem[ialpha].real();
				materials.back()[1] = rho_elem[ialpha].imag();
				rho[nuline + ialpha][0] = rho_elem[ialpha].real();
				rho[nuline + ialpha][1] = rho_elem[ialpha].imag();
				rho1D[nuline + ialpha][0] = rho_elem1D[ialpha].real();
				rho1D[nuline + ialpha][1] = rho_elem1D[ialpha].imag();
			}
			if (ialpha<vph_elem.size()) {
				materials.back()[2] = vp_elem[ialpha].real();
				materials.back()[3] = vp_elem[ialpha].imag();
				vp[nuline + ialpha][0] = vp_elem[ialpha].real();
				vp[nuline + ialpha][1] = vp_elem[ialpha].imag();
				vp1D[nuline + ialpha][0] = vp_elem1D[ialpha].real();
				vp1D[nuline + ialpha][1] = vp_elem1D[ialpha].imag();
			}
			if (ialpha<vpv_elem.size()) {
				materials.back()[4] = vpv_elem[ialpha].real();
				materials.back()[5] = vpv_elem[ialpha].imag();
			}
			if (ialpha<vsh_elem.size()) {
				materials.back()[6] = vsh_elem[ialpha].real();
				materials.back()[7] = vsh_elem[ialpha].imag();
			}
			if (ialpha<vsv_elem.size()) {
				materials.back()[8] = vsv_elem[ialpha].real();
				materials.back()[9] = vsv_elem[ialpha].imag();
			}
			if (ialpha<eta_elem.size()) {
				materials.back()[10] = eta_elem[ialpha].real();
				materials.back()[11] = eta_elem[ialpha].imag();
			}

			
		}
		nuline += Nus[i];
		
		// coords 
		mesh_S[i] = (elem->getCoordsPoints())[0];
		mesh_Z[i] = (elem->getCoordsPoints())[1];
		XMath::structuredUseFirstRow( (quad->getIntegralFactor()).cast <Real> (), integral_factor[i] ); //integral factor 


		
    }

	// ---------- </GET FIELDS> ----------	


	// ---------- <DEFINE WRITING OFFSETS> ----------
	
	int temp_startElem, temp_startElemNu;
	std::vector<int> temp_countElem(XMPI::nproc(),0);
	std::vector<int> temp_countElemNu(XMPI::nproc(),0);

	XMPI::gather(elems_proc, temp_countElem, true); // non fourier fields 
	temp_startElem = std::accumulate(temp_countElem.begin(), temp_countElem.begin()+XMPI::rank(), 0);
	
	XMPI::gather(elemNus_proc, temp_countElemNu, true); //fourier fields 
	temp_startElemNu = std::accumulate(temp_countElemNu.begin(), temp_countElemNu.begin()+XMPI::rank(),0);

	std::vector<size_t> startElem, countElem;
	startElem.push_back(temp_startElem);
	countElem.push_back(temp_countElem[XMPI::rank()]);
	
	std::vector<size_t> startElemGll, countElemGll;
	startElemGll.push_back(temp_startElem);
	countElemGll.push_back(temp_countElem[XMPI::rank()]);
	startElemGll.push_back(0);
	countElemGll.push_back(nPntEdge);
	startElemGll.push_back(0);
	countElemGll.push_back(nPntEdge);
	
	std::vector<size_t> startElemNuGll, countElemNuGll;
	startElemNuGll.push_back(temp_startElemNu);
	countElemNuGll.push_back(temp_countElemNu[XMPI::rank()]);
	startElemNuGll.push_back(0);
	countElemNuGll.push_back(12); //real/imag for each of the 6 material fields 
	startElemNuGll.push_back(0);
	countElemNuGll.push_back(nPntEdge);
	startElemNuGll.push_back(0);
	countElemNuGll.push_back(nPntEdge);
	

	int temp_startElem_ani, temp_startElemNu_ani;
	std::vector<int> temp_countElem_ani(XMPI::nproc(),0);
	std::vector<int> temp_countElemNu_ani(XMPI::nproc(),0);

	XMPI::gather(elems_proc_ani, temp_countElem_ani, true); //non fourier fields
	temp_startElem_ani = std::accumulate(temp_countElem_ani.begin(), temp_countElem_ani.begin()+XMPI::rank(), 0);

	XMPI::gather(elemNus_proc_ani, temp_countElemNu_ani, true); //fourier fields 
	temp_startElemNu_ani = std::accumulate(temp_countElemNu_ani.begin(), temp_countElemNu_ani.begin()+XMPI::rank(),0);
	
	std::vector<size_t> startElem_ani, countElem_ani;
	startElem_ani.push_back(temp_startElem_ani);
	countElem_ani.push_back(temp_countElem_ani[XMPI::rank()]);
	
	std::vector<size_t> startElemGll_ani, countElemGll_ani;
	startElemGll_ani.push_back(temp_startElem_ani);
	countElemGll_ani.push_back(temp_countElem_ani[XMPI::rank()]);
	startElemGll_ani.push_back(0);
	countElemGll_ani.push_back(gll_ani);
	startElemGll_ani.push_back(0);
	countElemGll_ani.push_back(gll_ani);
	
	std::vector<size_t> startElemNuGll_ani, countElemNuGll_ani;
	startElemNuGll_ani.push_back(temp_startElemNu);
	countElemNuGll_ani.push_back(temp_countElemNu[XMPI::rank()]);
	startElemNuGll_ani.push_back(0);
	countElemNuGll_ani.push_back(2);
	startElemNuGll_ani.push_back(0);
	countElemNuGll_ani.push_back(nPntEdge);
	startElemNuGll_ani.push_back(0);
	countElemNuGll_ani.push_back(nPntEdge);

	// ---------- </DEFINE WRITING OFFSETS> ----------

	// ---------- <DEFINE NETCDF FILE> ----------
	NetCDF_Writer nc_writer = NetCDF_Writer();	
	
	
	#ifdef _USE_PARALLEL_NETCDF // use parallel netcdf
		if (XMPI::root()) { 
			
			std::vector<size_t> dimsElem, dimsElemGll, dimsElemNuGll;
			
			dimsElem.push_back(elems_all);	
			dimsElemGll.push_back(elems_all);
			dimsElemGll.push_back(nPntEdge);
			dimsElemGll.push_back(nPntEdge);	
			dimsElemNuGll.push_back(elemNus_all);
			dimsElemNuGll.push_back(12); // real/imag for each of the 6 material fields 
			dimsElemNuGll.push_back(nPntEdge);
			dimsElemNuGll.push_back(nPntEdge);
			
			std::vector<size_t> dimsElem_ani, dimsElemGll_ani, dimsElemNuGll_ani;
			
			dimsElem_ani.push_back(elems_all_ani);	
			dimsElemGll_ani.push_back(elems_all_ani);
			dimsElemGll_ani.push_back(gll_ani);
			dimsElemGll_ani.push_back(gll_ani);
			dimsElemNuGll_ani.push_back(elemNus_all);
			dimsElemNuGll_ani.push_back(2); //just use this to test vp
			dimsElemNuGll_ani.push_back(nPntEdge);
			dimsElemNuGll_ani.push_back(nPntEdge);


			nc_writer.open(mFileName + ".nc4",true);	
			nc_writer.defModeOn();
			nc_writer.defineVariable<int>("Nus", dimsElem);
			nc_writer.defineVariable<int>("Nrs", dimsElem);
			nc_writer.defineVariable<int>("domain_decomposition", dimsElem);
			nc_writer.defineVariable<int>("element_mesh", dimsElem);
			nc_writer.defineVariable<int>("sem_mesh", dimsElemGll);
			nc_writer.defineVariable<Real>("mesh_S", dimsElemGll);
			nc_writer.defineVariable<Real>("mesh_Z", dimsElemGll);
			nc_writer.defineVariable<Real>("integral_factor", dimsElemGll);
			nc_writer.defineVariable<Real>("material_fields", dimsElemNuGll);
			nc_writer.defineVariable<Real>("vp", dimsElemNuGll_ani);
			nc_writer.defineVariable<Real>("vp1D", dimsElemNuGll_ani);
			nc_writer.defineVariable<Real>("rho", dimsElemNuGll_ani);
			nc_writer.defineVariable<Real>("rho1D", dimsElemNuGll_ani);




			nc_writer.defModeOff();
			
			nc_writer.writeVariableWhole("domain_decomposition", mMesh->mMsgInfo->mElemToProc);
			
			nc_writer.close();
		}
		// ---------- </DEFINE NETCDF FILE> ----------
		
					
		// ---------- <WRITE FIELDS> ----------	
		nc_writer.openParallel(mFileName + ".nc4");

		nc_writer.writeVariableChunk("Nus", Nus, startElem, countElem);
		nc_writer.writeVariableChunk("Nrs", Nrs, startElem, countElem);
		nc_writer.writeVariableChunk("element_mesh", mMesh->mMsgInfo->mLocElemToGlobElem, startElem, countElem);
		nc_writer.writeVariableChunk("sem_mesh", mMesh->mMsgInfo->mLocElemToGlobPoints, startElemGll, countElemGll);
		nc_writer.writeVariableChunk("mesh_S", mesh_S, startElemGll, countElemGll);
		nc_writer.writeVariableChunk("mesh_Z", mesh_Z, startElemGll, countElemGll);
		nc_writer.writeVariableChunk("integral_factor", integral_factor, startElemGll, countElemGll);
		nc_writer.writeVariableChunk("material_fields", materials, startElemNuGll, countElemNuGll);
		nc_writer.writeVariableChunk("vp", vp, startElemNuGll_ani, countElemNuGll_ani);
		nc_writer.writeVariableChunk("vp1D", vp1D, startElemNuGll_ani, countElemNuGll_ani);
		nc_writer.writeVariableChunk("rho", rho, startElemNuGll_ani, countElemNuGll_ani);
		nc_writer.writeVariableChunk("rho1D", rho1D, startElemNuGll_ani, countElemNuGll_ani);
	
	#else /// use serial netcdf 
		std::vector<size_t> dimsElem, dimsElemGll, dimsElemNuGll;
		std::vector<size_t> dimsElemAll; //this is just for domain decomposition, where each proc has the whole domain.
		
		dimsElemAll.push_back(elems_all);
		dimsElem.push_back(elems_proc);	
		dimsElemGll.push_back(elems_proc);
		dimsElemGll.push_back(nPntEdge);
		dimsElemGll.push_back(nPntEdge);	
		dimsElemNuGll.push_back(elemNus_proc);
		dimsElemNuGll.push_back(12); // real/imag for each of the 6 material fields 
		dimsElemNuGll.push_back(nPntEdge);
		dimsElemNuGll.push_back(nPntEdge);



		nc_writer.open(mFileName + "_" + std::to_string(XMPI::rank())+".nc4",true);	// create filenames with proc ranks as indices. 
		
		nc_writer.defModeOn();
		
		nc_writer.defineVariable<int>("Nus", dimsElem);
		nc_writer.defineVariable<int>("Nrs", dimsElem);
		nc_writer.defineVariable<int>("element_mesh", dimsElem);
		nc_writer.defineVariable<int>("sem_mesh", dimsElemGll);
		nc_writer.defineVariable<Real>("mesh_S", dimsElemGll);
		nc_writer.defineVariable<Real>("mesh_Z", dimsElemGll);
		nc_writer.defineVariable<Real>("integral_factor", dimsElemGll);
		nc_writer.defineVariable<Real>("material_fields", dimsElemNuGll);
		nc_writer.defineVariable<Real>("vp", dimsElemNuGll);
		nc_writer.defineVariable<Real>("vp1D", dimsElemNuGll);
		nc_writer.defineVariable<Real>("rho", dimsElemNuGll);
		nc_writer.defineVariable<Real>("rho1D", dimsElemNuGll);
		
		if (XMPI::root()) // we only need the domain decomposition on root 
			nc_writer.defineVariable<int>("domain_decomposition", dimsElemAll);

				
		nc_writer.defModeOff();
		
		if (XMPI::root()) {// we only need the domain decomposition on root 
			nc_writer.writeVariableWhole("domain_decomposition", mMesh->mMsgInfo->mElemToProc);
		}
		startElemNuGll[0] = 0;
		countElemNuGll[0] = elemNus_proc;
		startElemNuGll_ani[0] = 0;
		countElemNuGll_ani[0] = elemNus_proc;
			
		nc_writer.writeVariableWhole("Nus", Nus);
		nc_writer.writeVariableWhole("Nrs", Nrs);
		nc_writer.writeVariableWhole("element_mesh", mMesh->mMsgInfo->mLocElemToGlobElem);
		nc_writer.writeVariableWhole("sem_mesh", mMesh->mMsgInfo->mLocElemToGlobPoints);
		nc_writer.writeVariableWhole("mesh_S", mesh_S);
		nc_writer.writeVariableWhole("mesh_Z", mesh_Z);
		nc_writer.writeVariableWhole("integral_factor", integral_factor);

		// Another interesting netcdf behaviour. I should be able 
		// to use writeVariableWhole for the variables below. 
		// However, if they're 'too big', it crashes, but works with writeVariableChunk.
		nc_writer.writeVariableChunk("material_fields", materials, startElemNuGll, countElemNuGll);
		nc_writer.writeVariableChunk("vp", vp, startElemNuGll_ani, countElemNuGll_ani);
		nc_writer.writeVariableChunk("vp1D", vp1D, startElemNuGll_ani, countElemNuGll_ani);
		nc_writer.writeVariableChunk("rho", rho, startElemNuGll_ani, countElemNuGll_ani);
		nc_writer.writeVariableChunk("rho1D", rho1D, startElemNuGll_ani, countElemNuGll_ani);
		


		
	#endif



	nc_writer.close();
	// ---------- </WRITE FIELDS> ----------	


	
}