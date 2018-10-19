// Alex
// Dumps various quantities necessary for animations and/or kernels

#pragma once


#include <string>
class Mesh;
class Parameters;
class Domain;
class Source;

class MeshIO {
	
public:	
	MeshIO(const Mesh *mesh, const std::string &fname); 
	void dumpFields(const Domain &domain, const Parameters &par);

private:
	const Mesh *mMesh;
	std::string mFileName;
	
};