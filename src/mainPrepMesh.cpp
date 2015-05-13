#include "Reader.h"
#include "VDBWrapper.h"
#include "PostProcessor.h"

#include <tclap/CmdLine.h>

using namespace FractureSim;
using namespace std;

#include <fstream>
#include <sstream>
int assignRegions(id_map& regions, node_map& regionDefs, node_map& nodes, elem_map& elems);
int readRegionDefinitions(node_map& regionDefs, string filename); //node_map maps id to vector<double>

int main(int argc, char* argv[]){
	try {
		printf("%% parsing command line ...\n");
		TCLAP::CmdLine cmd("Preprocessing mesh", ' ', "using OpenVDB 2.2 backend");
		TCLAP::UnlabeledValueArg<string> inpFile("filename", 
			"Elmer format input file containing the BEM mesh", true, "default", "filename",
			cmd
		);
		TCLAP::ValueArg<string> outFile("o","out","VTK/VBD/mesh output file/dir name", false, "out","out-file", cmd);
		TCLAP::ValueArg<double> resMesh("r","res-mesh","resolution of BEM mesh, default is 0.1", false, 0.1,"res-mesh", cmd);
		TCLAP::ValueArg<double> resGrid("R","res-grid","resolution of VDB grid, default is 0.01", false, 0.01,"res-grid", cmd);
		cmd.parse( argc, argv );

		BEMReader reader(inpFile.getValue());
		// reading model
		id_set region_ids, //leave empty to read all regions
			   cracks;      //read these line-elements as crack tips
		//crack_tip_ids.insert(crackTipIDs.getValue().begin(), crackTipIDs.getValue().end()); //fill from CLI
		elem_map elems, ctElems, ctParents;
		node_map nodes;
		id_map regions;
		
		reader.readModel(
			elems, regions, ctElems, ctParents, nodes,
			FractureSim::TRI, region_ids, FractureSim::LINE, cracks
		);
		printf("%% read input file, have %d elems and %d nodes.\n",elems.size(),nodes.size());
		cracks.clear(); //cracks is now treated as set of regions representing fracture surfaces

		//pretend we didn't know regions from the input mesh and assign them based on a region-spec file
		node_map regionDefs;
		int nRegions=readRegionDefinitions(regionDefs,inpFile.getValue()+".regions");
		assignRegions(regions, regionDefs,nodes,elems);

		VDBWrapper ls(resGrid.getValue());
		ls.init(nodes,elems,regions,cracks);

		printf("%% level-set initialized.\n");

		ls.writeAll(outFile.getValue()+".vdb");

		state_map empty;
		PostProcessor* toVTK = new PostProcessor(nodes, elems, regions, cracks, ctElems, ctParents, empty);
		vector_type regionVect(elems.size());
		for(int i=0; i<elems.size(); ++i){
			regionVect[i]=regions[i+ELEM_BASE_INDEX];
		}
		toVTK->vtkAddData("region",1,true,regionVect);
		toVTK->vtkWrite(outFile.getValue()+".vtk",VTK_TRIANGLE);

		printf("%% output written.\n");

		ls.mesh(nodes,elems,resMesh.getValue());
		regions.clear(); cracks.clear(); ctElems.clear(); ctParents.clear();
		delete toVTK;
		toVTK = new PostProcessor(nodes, elems, regions, cracks, ctElems, ctParents, empty);

		assignRegions(regions, regionDefs,nodes,elems);
		regionVect.resize(elems.size());
		for(int i=0; i<elems.size(); ++i){
			regionVect[i]=regions[i+ELEM_BASE_INDEX];
		}
		toVTK->vtkAddData("region",1,true,regionVect);
		toVTK->vtkWrite(outFile.getValue()+"_remeshed.vtk",VTK_TRIANGLE);


	}catch (TCLAP::ArgException &e)  {// catch any exceptions
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return -1;
	}
	return 0;
}

int assignRegions(id_map& regions, node_map& regionDefs, node_map& nodes, elem_map& elems){
	regions.clear();
	Eigen::Vector3d n; // tmp storage for a region-constraint
	double d; // such that the constraint is evaluated on p as n.dot(p)<=d
	bool found, match;

	unsigned int i=0;
	for(elem_map::const_iterator it=elems.begin(); it!=elems.end(); ++it,++i){
		// check which region fits this tri
		Eigen::Vector3d // a,b,c are node coordinates
			a (nodes[it->second[0]][0], nodes[it->second[0]][1], nodes[it->second[0]][2]),
			b (nodes[it->second[1]][0], nodes[it->second[1]][1], nodes[it->second[1]][2]),
			c (nodes[it->second[2]][0], nodes[it->second[2]][1], nodes[it->second[2]][2]);
		found=false;
		for(node_map::iterator rd=regionDefs.begin(); rd!=regionDefs.end() && !found; ++rd){ //iterate region defs
			match=true;
			for(int j=0; j<rd->second.size() && match; j+=4){
				n[0]=rd->second[j  ]; n[1]=rd->second[j+1]; n[2]=rd->second[j+2];
				d   =rd->second[j+3];
				if( n.dot(a) > d && // || // using && means tri is added to region if at least 1 vertex matches
					n.dot(b) > d && // || // using || means tri is added to region if all of its vertices match
					n.dot(c) > d ) match=false;
			}
			if(match){
				found=true;
				regions[it->first]=rd->first;
			}
		}
	}
	return 0;
}

int readRegionDefinitions(node_map& regionDefs, string filename){ //node_map maps id to vector<double>
	int nRegions=0; // number of regionDefs read
	bool regKwdFound=false; // look for the keyword "regions" followed by the number of regions to read
	string line, token;
	int ret, done=0;
	char check;
	double a,b,c,d;
	unsigned int nextId;

	std::istringstream strstream;
	ifstream in(filename.c_str());
	if(!in.is_open()) return -1;

	while(in.good()){
		getline(in, line);
		if(line.empty() || line.at(0)=='#') continue; // comments have # in the first column
		strstream.clear();
		strstream.str(line);

		getline(strstream, token, ' ');
		if(!regKwdFound){
			if(token.compare("regions")==0){
				regKwdFound=true;
				getline(strstream, token, ' '); //next token is number of regions to read
				ret=sscanf(token.c_str(), "%u", &nRegions);
				if(ret!=1){
					in.close();
					//printf("invalid ret=%d should be 1 on token %s\n",ret, token.c_str());
					return -1;
				}
			}
		}else if(done<nRegions){ //reading next region until we have plenty
			//printf("processing '%s' done %d\n",line.c_str(),done);
			ret=sscanf(token.c_str(), "%u", &nextId);
			if(ret==1) while(getline(strstream, token, ' ')){ //token is now one condition of the region definition
				//printf("parsing token '%s'",token.c_str());
				ret=sscanf(token.c_str(), "(%lf,%lf,%lf,%lf%c",&a,&b,&c,&d,&check);
				if(ret==5 && check==')'){ // correct format
					regionDefs[nextId].push_back(a);
					regionDefs[nextId].push_back(b);
					regionDefs[nextId].push_back(c);
					regionDefs[nextId].push_back(d);
					//printf(" ... ok\n");
				}
				//else printf(" ... reject ret=%d, check='%c'\n",ret,check);
			}
			++done;
		}
	}
	//printf("\nregionDefs:\n");
	//printMap(regionDefs);
	// finally add an empty region def which will be the default
	if(regionDefs.empty()) nextId=ELEM_BASE_INDEX;
	else nextId = regionDefs.rbegin()->first +1;
	regionDefs[nextId].assign(0,0.0);

	return nRegions;
}