#include "PostProcessor.h"

#include <fstream>
#include <cmath>
#include <cfloat>

#include <iostream>

using namespace std;

namespace FractureSim{

	int PostProcessor::computeNodeSIFs(
            const vector_type& displacements, double E, double nu,
			vect3d_map& sifs, vect3d_map& faceNormals, vect3d_map& tangents
    ){
		// compute local coordinate systems and SIFs at all nodes on a crack-tip
		// uses information from the crack-tip parent elements (not necessarily all adjacent nodes)
		// - run along the crack-tip, for each segment, take the parent triangle
		// - compute the averaged face normal of parent-tris at each endpoint
		// - compute the averaged SIFs by displacement correlation from the interior nodes
		sifs.clear();
		faceNormals.clear();
		tangents.clear();
		id_map count; //printf("\n%% computing SIFs...");
		for(elem_map::iterator it=crackTips.begin(); it!=crackTips.end(); ++it){
            unsigned int nd_a, nd_b, nd_c; // nodes of the current triangle
            nd_a = it->second[0]; nd_b = it->second[1];
            nd_c = findInteriorNode(nd_a, nd_b, elems[parents[it->first][0]]);
            Eigen::Vector3d a,b,c, n1,n2,n3;
            copyNode(nd_a, a);
            copyNode(nd_b, b);
            copyNode(nd_c, c);
            getLocalCoordFrame(a,b,c, n1,n2,n3);

            if(tangents.count(nd_a)==0) tangents[nd_a]=Eigen::Vector3d(0.0,0.0,0.0);
            if(tangents.count(nd_b)==0) tangents[nd_b]=Eigen::Vector3d(0.0,0.0,0.0);
            if(faceNormals.count(nd_a)==0) faceNormals[nd_a]=Eigen::Vector3d(0.0,0.0,0.0);
            if(faceNormals.count(nd_b)==0) faceNormals[nd_b]=Eigen::Vector3d(0.0,0.0,0.0);
            faceNormals[nd_a]+=n1;
            faceNormals[nd_b]+=n1;
            tangents[nd_a]+=n3;
            tangents[nd_b]+=n3;
		}
		for(vect3d_map::iterator it=faceNormals.begin();
			it!=faceNormals.end(); ++it){ // normalize face normals
			it->second.normalize();
			//printf("\n .. normal  at node %d:\t(%.3lf, %.3lf, %.3lf)", it->first, it->second[0], it->second[1], it->second[2]);
		}
		for(vect3d_map::iterator it=tangents.begin();
			it!=tangents.end(); ++it){ // normalize tangents
			it->second.normalize();
			//printf("\n .. tangent at node %d:\t(%.3lf, %.3lf, %.3lf)", it->first, it->second[0], it->second[1], it->second[2]);
		}
		// now we have surface normals and crack-tip tangents at all crack-tip nodes
		// next: use displacement correlation to compute SIFs
		for(elem_map::iterator it=crackTips.begin(); it!=crackTips.end(); ++it){
            unsigned int nd_a, nd_b, nd_c; // nodes of the current triangle
            nd_a = it->second[0]; nd_b = it->second[1];
            nd_c = findInteriorNode(nd_a, nd_b, elems[parents[it->first][0]]);
			Eigen::Vector3d a,b,c, n1,n2,n3, cod; // crack opening displacement at correlation point
			copyNode(nd_a, a);
            copyNode(nd_b, b);
            copyNode(nd_c, c);
            getLocalCoordFrame(a,b,c, n1,n2,n3);
            // correlation point is the interior node (nd_c)
            cod[0] = displacements[3*(nd_c-NODE_BASE_INDEX)  ];
            cod[1] = displacements[3*(nd_c-NODE_BASE_INDEX)+1];
            cod[2] = displacements[3*(nd_c-NODE_BASE_INDEX)+2];
			double r; // distance from the crack tip to the correlation point
            r=((0.5*(a+b))-c).dot(n2); //project midpoint-interior to edge normal (shortest distance from correlation pt to crack tip)
            double mu = E/(2*(1+nu));
            double u1,u2,u3; // projected cod in n1,n2,n3 directions
            u1=cod.dot(faceNormals[nd_a]);
			u2=cod.dot(tangents[nd_a].cross(faceNormals[nd_a]));
            u3=cod.dot(tangents[nd_a]);
			if(sifs.count(nd_a)==0) sifs[nd_a].setZero();
            sifs[nd_a][0] += -mu*sqrt(2*M_PI)*u1/(sqrt(r)*2*(1-nu)); //flip signs as opening displacement is opposite normals
            sifs[nd_a][1] += -mu*sqrt(2*M_PI)*u2/(sqrt(r)*2*(1-nu));
            sifs[nd_a][2] += -mu*sqrt(M_PI)*u3/sqrt(2*r);
			//printf("\n .. intermediate K1 at node %d is %.1le", nd_a, -mu*sqrt(2*M_PI)*u1/(sqrt(r)*2*(1-nu)) );
			++count[nd_a];
			
            u1=cod.dot(faceNormals[nd_b]);
			u2=cod.dot(tangents[nd_b].cross(faceNormals[nd_b]));
            u3=cod.dot(tangents[nd_b]);
			if(sifs.count(nd_b)==0) sifs[nd_b].setZero();
            sifs[nd_b][0] += -mu*sqrt(2*M_PI)*u1/(sqrt(r)*2*(1-nu)); //flip signs as opening displacement is opposite normals
            sifs[nd_b][1] += -mu*sqrt(2*M_PI)*u2/(sqrt(r)*2*(1-nu));
            sifs[nd_b][2] += -mu*sqrt(M_PI)*u3/sqrt(2*r);
			//printf("\n .. intermediate K1 at node %d is %.1le", nd_b, -mu*sqrt(2*M_PI)*u1/(sqrt(r)*2*(1-nu)) );
			++count[nd_b];
		}
		for(vect3d_map::iterator it=sifs.begin(); it!=sifs.end(); ++it)
			it->second /= count[it->first]; // divide SIFs by number of contributions

		return 0;
	}
    
    int PostProcessor::vtkWrite(std::string fileName, VTK_CELL_TYPE cellType){
        int nNodes = nodes.size();
        int nCoordsPerNode = nodes.begin()->second.size();
        int nElems = elems.size();
        int nNodesPerElem = elems.begin()->second.size();
		ofstream out(fileName.c_str());
        
        if(!out.is_open()) return -1;
        out.precision(12);
        //out.setf(std::ios::scientific);
        //write header
        out << "# vtk DataFile Version 2.0" << endl;
        out << "VTK exported mesh" << endl;
        out << "ASCII" << endl;
        out << "DATASET UNSTRUCTURED_GRID" << endl;
        //node coordinates
        out << "POINTS " << nNodes << " double" << endl;
//        for(int i=0;i<nNodes; ++i){
        for(node_map::iterator i=nodes.begin(); i!=nodes.end(); ++i){
//            out << nodes[3*i] << " " << nodes[3*i+1] << " " << nodes[3*i+2] << endl;
            for(int j=0; j<nCoordsPerNode; ++j) out << i->second[j] << " ";
            for(int j=nCoordsPerNode; j<3; ++j) out << "0.0 "; // fill with zeros if less than 3 coords per node given
            out << endl;
        }
        //cells
        out << "CELLS " << nElems << " " <<  (1+nNodesPerElem)*nElems << endl;
//        for(int i=0;i<nElems; ++i){
        for(elem_map::iterator i=elems.begin(); i!=elems.end(); ++i){
            out << nNodesPerElem << " ";
            // VTK Unstructured Grid identifies nodes by 0-based index into POINTS list
            for(int j=0; j<nNodesPerElem; ++j) out << i->second[j]-NODE_BASE_INDEX << " ";
            out << endl;
        }
        //cell types
        out << "CELL_TYPES " << nElems << endl;
        for(int i=0;i<nElems; ++i) out << cellType << endl;
        // mesh information is done now, next up point and cell data ...
        
        out << "POINT_DATA " << nNodes << endl;
        for(int i=0;i<vtkDataNames.size(); ++i){
            if(vtkDataIsCellData[i]==false) // this is point data
                vtkWriteData(out, i, nNodes);
        }

        out << "CELL_DATA " << nElems << endl;
        for(int i=0;i<vtkDataNames.size(); ++i){
            if(vtkDataIsCellData[i]==true) // this is cell data
                vtkWriteData(out, i, nElems);
        }

        return 0;
    }
    
    void PostProcessor::vtkWriteData(ofstream& out, int i, int n){
        if(vtkDataDimension[i]==1){ //scalars
            out << "SCALARS " << vtkDataNames[i] << " double" << endl;
            out << "LOOKUP_TABLE default" << endl; //could specify color lookup
            for(int j=0;j<n; ++j) out << (*vtkData[i])[j] << endl;
        }else if(vtkDataDimension[i]==3){ //vectors
            out << "VECTORS " << vtkDataNames[i] << " double" << endl;
            for(int j=0;j<n; ++j)
                out << (*vtkData[i])[3*j  ] << " "
                    << (*vtkData[i])[3*j+1] << " "
                    << (*vtkData[i])[3*j+2] << endl;
        }
    }
    
    int PostProcessor::vtkAddData(
            std::string name, int dimension, bool isCellData, vector_type& data
    ){
        if( //check that number of data entries matches expected value
            (isCellData==true  && data.rows()==dimension*elems.size()) ||
            (isCellData==false && data.rows()==dimension*nodes.size())
        ){
            vtkDataNames.push_back(name);
            vtkDataDimension.push_back(dimension);
            vtkDataIsCellData.push_back(isCellData);
            vtkData.push_back(&data);
            return 0;
        }else{
//            fprintf(stderr, "can't add %s data because data size is %d but expected %d (dim=%d)\n",
//                isCellData?"cell":"point",data.rows(),
//                isCellData?(dimension*elems.size()):(dimension*nodes.size()), dimension);
            return -1;
        }
    }

	/** // node based version -- OUTDATED - USE THE ELEMENT BASED VERSION!
	int PostProcessor::computeSurfaceStresses(
        node_map& maxPrincipalStress, const vector_type& displacements,
        double E, double nu, bool addToVTK
    ){
		// compute principal and cartesian stresses based on the linear FEM formulation found in
		// Bargteil et al. 2007 "A Finite Element Method for Animating Large Viscoplastic Flow"
		// by extending every triangle to a tet with no out-of-plane deformation
		// ...
		// average deformation gradient from elements to nodes, area weighted
		unsigned int n = nodes.size();
		s_xx.resize(n); s_xx.setZero();
		s_yy.resize(n); s_yy.setZero();
		s_zz.resize(n); s_zz.setZero();
		s_xy.resize(n); s_xy.setZero();
		s_xz.resize(n); s_xz.setZero();
		s_yz.resize(n); s_yz.setZero();
        maxPrincipalStress.clear();
		std::map<unsigned int, double> nd_area;
		std::map<unsigned int, Eigen::Matrix3d> nd_F;
		id_set ctNodes = nodeSet(crackTips); // we want to ignore elements that contain a crack-tip node here, since they should be handled by stress intensity.

		//matrices are: B inverse of edge vector matrix in material space, X edge vector matrix in world space, F deformation gradient
		// U*S*V' SVD of F, P stress (Piola-Kirchhoff)
		Eigen::Matrix3d B,F,U,S,Vt,P,X, I=Eigen::Matrix3d::Identity();
		double mu=E/(2*(1+nu)), lambda=E*nu/((1+nu)*(1-2*nu)); // convert (E, nu) to Lame parameters (mu = shear modulus, lambda) ; as usual expect trouble for nu=0.5

		for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it){ // loop over elements			
			bool flag=false;
			if( ctNodes.count(it->second[0]) !=0 ||
				ctNodes.count(it->second[1]) !=0 ||
				ctNodes.count(it->second[2]) !=0 ) continue; // skip crack-tip adjacent elements
			Eigen::Vector3d // a,b,c are node coordinates in material space; ua,ub,uc are nodal displacements, so world coordinates are a+ua etc.
				a (nodes[it->second[0]][0], nodes[it->second[0]][1], nodes[it->second[0]][2]),
				b (nodes[it->second[1]][0], nodes[it->second[1]][1], nodes[it->second[1]][2]),
				c (nodes[it->second[2]][0], nodes[it->second[2]][1], nodes[it->second[2]][2]),
				ua(u[3*(it->second[0]-NODE_BASE_INDEX)], u[3*(it->second[0]-NODE_BASE_INDEX)+1], u[3*(it->second[0]-NODE_BASE_INDEX)+2]),
				ub(u[3*(it->second[1]-NODE_BASE_INDEX)], u[3*(it->second[1]-NODE_BASE_INDEX)+1], u[3*(it->second[1]-NODE_BASE_INDEX)+2]),
				uc(u[3*(it->second[2]-NODE_BASE_INDEX)], u[3*(it->second[2]-NODE_BASE_INDEX)+1], u[3*(it->second[2]-NODE_BASE_INDEX)+2]);
            if(cracks.count(regions[it->first])!=0){ // this is a crack element
                ua*=0.5; ub*=0.5; uc*=0.5; // use half of crack opening displacement (because COD should be distributed to 2 surfaces)
            }
			X.col(0) = a - c; // using world space matrix as temporary storage before inversion
			X.col(1) = b - c;
			X.col(2) = X.col(0).cross(X.col(1)).normalized();
			X.computeInverseWithCheck(B,flag);
			if(!flag) return -1; // matrix not invertible ~> possibly degenerate element
			// now build the actual world space matrix
			X.col(0) = a+ua -c-uc;
			X.col(1) = b+ub -c-uc;
			X.col(2) = X.col(0).cross(X.col(1)).normalized();
			F=X*B; //F= 0.5*( F+F.transpose());

			// add 1/3 of F to each node of the triangle
			double tri_area_3 = ((a-c).cross(b-c).norm())/6.0; // 1/3 of triangle area
			for(int k=0; k<3; ++k){
				if(nd_area.count(it->second[k])==0){
					nd_area[it->second[k]]=0.0;
					nd_F[it->second[k]].setZero();
				}
				nd_area[it->second[k]] += tri_area_3;
				nd_F[it->second[k]] += (F *tri_area_3);

				if(!(F(0,0)==F(0,0) && F(0,1)==F(0,1) && F(0,2)==F(0,2) &&
					 F(1,0)==F(1,0) && F(1,1)==F(1,1) && F(1,2)==F(1,2) &&
					 F(2,0)==F(2,0) && F(2,1)==F(2,1) && F(2,2)==F(2,2))||
					 std::abs(F(0,0))>DBL_MAX || std::abs(F(0,1))>DBL_MAX || std::abs(F(0,2))>DBL_MAX ||  
					 std::abs(F(1,0))>DBL_MAX || std::abs(F(1,1))>DBL_MAX || std::abs(F(2,2))>DBL_MAX ||  
					 std::abs(F(2,0))>DBL_MAX || std::abs(F(2,1))>DBL_MAX || std::abs(F(2,2))>DBL_MAX )
					 std::cout << "\nerror in F: " << F << "; element " << it->first << "; a=" << a << " b=" << b << " c=" << c << "; ua=" << ua << " ub=" << ub << " uc=" << uc << std::endl;   
			}
		}

		// compute nodal stresses from interpolated deformation gradient
		unsigned int i=0;
		for(node_map::iterator it= nodes.begin(); it!=nodes.end(); ++it, ++i){
			if(nd_F.count(it->first)!=0){	
				F=nd_F[it->first];
				if(!(F(0,0)==F(0,0) && F(0,1)==F(0,1) && F(0,2)==F(0,2) &&
					 F(1,0)==F(1,0) && F(1,1)==F(1,1) && F(1,2)==F(1,2) &&
					 F(2,0)==F(2,0) && F(2,1)==F(2,1) && F(2,2)==F(2,2))||
					 std::abs(F(0,0))>DBL_MAX || std::abs(F(0,1))>DBL_MAX || std::abs(F(0,2))>DBL_MAX ||  
					 std::abs(F(1,0))>DBL_MAX || std::abs(F(1,1))>DBL_MAX || std::abs(F(2,2))>DBL_MAX ||  
					 std::abs(F(2,0))>DBL_MAX || std::abs(F(2,1))>DBL_MAX || std::abs(F(2,2))>DBL_MAX)
					 std::cout << "\npre-div error in F: " << F << "; node " << it->first << " nd_area=" << nd_area[it->first] <<  std::endl;   
			}

			// since we skip crack-tip adjacent elements, some nodes may not have any associated area
			if(nd_area.count(it->first)==0) nd_F[it->first].setIdentity();
			else if(std::abs(1.0/nd_area[it->first])>DBL_MAX) nd_F[it->first].setIdentity();
			else F = ((1.0/nd_area[it->first])* nd_F[it->first]); //nd_F[it->first] /= nd_area[it->first];

			//F=nd_F[it->first];
			if(!(F(0,0)==F(0,0) && F(0,1)==F(0,1) && F(0,2)==F(0,2) &&
				 F(1,0)==F(1,0) && F(1,1)==F(1,1) && F(1,2)==F(1,2) &&
				 F(2,0)==F(2,0) && F(2,1)==F(2,1) && F(2,2)==F(2,2))||
				 std::abs(F(0,0))>DBL_MAX || std::abs(F(0,1))>DBL_MAX || std::abs(F(0,2))>DBL_MAX ||  
				 std::abs(F(1,0))>DBL_MAX || std::abs(F(1,1))>DBL_MAX || std::abs(F(2,2))>DBL_MAX ||  
				 std::abs(F(2,0))>DBL_MAX || std::abs(F(2,1))>DBL_MAX || std::abs(F(2,2))>DBL_MAX )
				 std::cout << "\nerror in F: " << F << "; node " << it->first << " nd_area=" << nd_area[it->first] <<  std::endl;   

			Eigen::JacobiSVD<Eigen::Matrix3d> svd(F,Eigen::ComputeFullU | Eigen::ComputeFullV);
			S = svd.singularValues().asDiagonal(); // these will be returned in order (highest to lowest)
			U = svd.matrixU();
			Vt= svd.matrixV().transpose();
			P = 2*mu*(S-I) + lambda*(S-I).trace()*I;
			if(!(P(0,0)==P(0,0))){ printf("\nerror in P(0,0)=%.3le\n",P(0,0)); std::cout << "F: " << nd_F[it->first] << std::endl;}
			if(!(Vt(0,0)==Vt(0,0))){ printf("\nerror in Vt(0,0)=%.3le -- node %d\n",Vt(0,0), it->first); }
            maxPrincipalStress[it->first].assign(8,0.0);
            maxPrincipalStress[it->first][0]=P(0,0); // first entry is the max principal stress value
			maxPrincipalStress[it->first][1]=Vt(0,0); // next 3 entries
            maxPrincipalStress[it->first][2]=Vt(0,1); // are the plane-
            maxPrincipalStress[it->first][3]=Vt(0,2); // normal vector
			maxPrincipalStress[it->first][4]=P(2,2);
			maxPrincipalStress[it->first][5]=Vt(2,0); // next 3 entries
            maxPrincipalStress[it->first][6]=Vt(2,1); // are the plane-
            maxPrincipalStress[it->first][7]=Vt(2,2); // normal vector
			s_xx[i]=P(0,0); s_yy[i]=P(1,1); s_zz[i]=P(2,2); // testing only, store principal stress values
			
			P = U*P*Vt;
			if(addToVTK){
				//s_xx[i]=P(0,0); s_yy[i]=P(1,1); s_zz[i]=P(2,2);
				s_xy[i]=0.5*( P(0,1)+P(1,0) ); s_xz[i]=0.5*( P(0,2)+P(2,0) ); s_yz[i]=0.5*( P(2,1)+P(1,2) );
			}
		}
        if(addToVTK){
            vtkAddData("stress_xx",1,false,s_xx);
            vtkAddData("stress_yy",1,false,s_yy);
            vtkAddData("stress_zz",1,false,s_zz);
            vtkAddData("stress_xy",1,false,s_xy);
            vtkAddData("stress_xz",1,false,s_xz);
            vtkAddData("stress_yz",1,false,s_yz);
        }
		return 0;
	}
	/*/ // element based version
	int PostProcessor::computeSurfaceStresses(
        node_map& maxPrincipalStress, const vector_type& u, const vector_type& u_c,
        double E, double nu, bool addToVTK
    ){
		// compute principal and cartesian stresses based on the linear FEM formulation found in
		// Bargteil et al. 2007 "A Finite Element Method for Animating Large Viscoplastic Flow"
		// by extending every triangle to a tet with no out-of-plane deformation
		// ...
		unsigned int n = elems.size(), i=0;
		s_xx.resize(n);
		s_yy.resize(n);
		s_zz.resize(n);
		s_xy.resize(n);
		s_xz.resize(n);
		s_yz.resize(n);
        maxPrincipalStress.clear();

		//matrices are: B inverse of edge vector matrix in material space, X edge vector matrix in world space, F deformation gradient
		// U*S*V' SVD of F, P stress (Piola-Kirchhoff)
		Eigen::Matrix3d B,F,U,S,Vt,P,X, I=Eigen::Matrix3d::Identity();

		for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it, ++i){ // loop over elements
			double mu=E/(2*(1+nu)), lambda=E*nu/((1+nu)*(1-2*nu)); // convert (E, nu) to Lame parameters (mu = shear modulus, lambda) ; as usual expect trouble for nu=0.5
			Eigen::Vector3d // a,b,c are node coordinates in material space; ua,ub,uc are nodal displacements, so world coordinates are a+ua etc.
				a (nodes[it->second[0]][0], nodes[it->second[0]][1], nodes[it->second[0]][2]),
				b (nodes[it->second[1]][0], nodes[it->second[1]][1], nodes[it->second[1]][2]),
				c (nodes[it->second[2]][0], nodes[it->second[2]][1], nodes[it->second[2]][2]),
				ua(u[3*(it->second[0]-NODE_BASE_INDEX)], u[3*(it->second[0]-NODE_BASE_INDEX)+1], u[3*(it->second[0]-NODE_BASE_INDEX)+2]),
				ub(u[3*(it->second[1]-NODE_BASE_INDEX)], u[3*(it->second[1]-NODE_BASE_INDEX)+1], u[3*(it->second[1]-NODE_BASE_INDEX)+2]),
				uc(u[3*(it->second[2]-NODE_BASE_INDEX)], u[3*(it->second[2]-NODE_BASE_INDEX)+1], u[3*(it->second[2]-NODE_BASE_INDEX)+2]),
				uan,ubn,ucn; // for negative COD of fracture elements
            if(cracks.count(regions[it->first])!=0){ // this is a crack element
				Eigen::Vector3d // uba,ubb,ubc are nodal crack base displacements
					uba(u_c[3*(it->second[0]-NODE_BASE_INDEX)], u_c[3*(it->second[0]-NODE_BASE_INDEX)+1], u_c[3*(it->second[0]-NODE_BASE_INDEX)+2]),
					ubb(u_c[3*(it->second[1]-NODE_BASE_INDEX)], u_c[3*(it->second[1]-NODE_BASE_INDEX)+1], u_c[3*(it->second[1]-NODE_BASE_INDEX)+2]),
					ubc(u_c[3*(it->second[2]-NODE_BASE_INDEX)], u_c[3*(it->second[2]-NODE_BASE_INDEX)+1], u_c[3*(it->second[2]-NODE_BASE_INDEX)+2]);
                ua*=0.5; ub*=0.5; uc*=0.5; // use half of crack opening displacement (because COD is distributed over 2 surfaces)
				uan=uba-ua; ubn=ubb-ub; ucn=ubc-uc; // negative COD + base
				ua+=uba; ub+=ubb; uc+=ubc;          // positive COD + base displacement
            }

			if( computeElementPrincipalStresses(U,S,Vt,P, a,b,c, ua,ub,uc, mu,lambda) !=0) return -1;

            maxPrincipalStress[it->first].assign(9,0.0);
			// store max. principal stress and its plane-normal vector for each element (also min. for compressive fracture?)
            // in the SVD, V is the pre-rotation and U is the post-rotation
            // so the direction of max. principal stress is the first column of V
            // which becomes the first row of V.transpose()
            maxPrincipalStress[it->first][0]=P(0,0);  // first entry is the stress value
            maxPrincipalStress[it->first][1]=Vt(0,0); // next 3 entries
            maxPrincipalStress[it->first][2]=Vt(0,1); // are the plane-
            maxPrincipalStress[it->first][3]=Vt(0,2); // normal vector
            // same for min. principal stress used for compressive fractures
			maxPrincipalStress[it->first][4]=P(2,2);
			maxPrincipalStress[it->first][5]=Vt(2,0); // next 3 entries
            maxPrincipalStress[it->first][6]=Vt(2,1); // are the plane-
            maxPrincipalStress[it->first][7]=Vt(2,2); // normal vector
            //s_xx[i]=P(0,0); s_yy[i]=P(1,1); s_zz[i]=P(2,2); // for testing

			// flag fracture elements based on sign of COD
			// [8] --> 1: both positive; 2: max. negative, min. positive; 3: max. pos., min. neg.; 4: both neg.
			if(cracks.count(regions[it->first])!=0){
				maxPrincipalStress[it->first][8]=1.0;
				bool tmp=false;
				// repeat stress computation with crack-base displacement MINUS half crack-opening displacement and see which sign of the COD produces more principal stress
				if( computeElementPrincipalStresses(U,S,Vt,P, a,b,c, uan,ubn,ucn, mu,lambda) !=0) return -2;
				if( std::abs(maxPrincipalStress[it->first][0]) < std::abs(P(0,0)) ){
					maxPrincipalStress[it->first][0]=P(0,0);  // first entry is the stress value
					maxPrincipalStress[it->first][1]=Vt(0,0); // next 3 entries
					maxPrincipalStress[it->first][2]=Vt(0,1); // are the plane-
					maxPrincipalStress[it->first][3]=Vt(0,2); // normal vector
					maxPrincipalStress[it->first][8]=2.0; tmp=true;
					//printf("\n%% *** have neg. ten.");
				}
				if( std::abs(maxPrincipalStress[it->first][4]) < std::abs(P(2,2)) ){
					maxPrincipalStress[it->first][4]=P(2,2);
					maxPrincipalStress[it->first][5]=Vt(2,0); // next 3 entries
					maxPrincipalStress[it->first][6]=Vt(2,1); // are the plane-
					maxPrincipalStress[it->first][7]=Vt(2,2); // normal vector
					maxPrincipalStress[it->first][8]=tmp?(4.0):(3.0);
					//printf("\n%% *** have neg. comp.");
				}
			}else{ //not a fracture element
				// for testing: only consider branching by setting non-fracture stress to 0
				//printf("\n%% warning: we're executing experimental code here ...");
				//maxPrincipalStress[it->first][0]=0.0;
				//maxPrincipalStress[it->first][4]=0.0;
			}

            P = U*P*Vt; // P is now a matrix of cartesian stresses
			//printf("\n%% -- cartesian: %.3le, %.3le, %.3le", P(0,0), P(1,1), P(2,2));
			//printf("\n%% --  w/ shear: %.3le, %.3le, %.3le", P(0,1), P(0,2), P(1,2));
			if(addToVTK){
                s_xx[i]=P(0,0); s_yy[i]=P(1,1); s_zz[i]=P(2,2);
                s_xy[i]=0.5*( P(0,1)+P(1,0) ); s_xz[i]=0.5*( P(0,2)+P(2,0) ); s_yz[i]=0.5*( P(2,1)+P(1,2) );
            }
		}
        if(addToVTK){
            vtkAddData("stress_xx",1,true,s_xx);
            vtkAddData("stress_yy",1,true,s_yy);
            vtkAddData("stress_zz",1,true,s_zz);
            vtkAddData("stress_xy",1,true,s_xy);
            vtkAddData("stress_xz",1,true,s_xz);
            vtkAddData("stress_yz",1,true,s_yz);
        }
		return 0;
	}
	/**/

    
	int PostProcessor::computeElementPrincipalStresses(
		Eigen::Matrix3d& U, Eigen::Matrix3d& S, Eigen::Matrix3d& Vt, Eigen::Matrix3d& P,
		const Eigen::Vector3d&  a, const Eigen::Vector3d&  b, const Eigen::Vector3d&  c,
		const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc,
		double mu, double lambda
	){
		Eigen::Matrix3d B,F,X, I=Eigen::Matrix3d::Identity();
		bool flag=true;
		X.col(0) = a - c; // using world space matrix as temporary storage before inversion
		X.col(1) = b - c;
		X.col(2) = X.col(0).cross(X.col(1)).normalized();
		X.computeInverseWithCheck(B,flag);
		if(!flag) return -1; // matrix not invertible ~> possibly degenerate element
		// now build the actual world space matrix
		X.col(0) = a+ua -c-uc;
		X.col(1) = b+ub -c-uc;
		X.col(2) = X.col(0).cross(X.col(1)).normalized();
		F=X*B;
		// compute singular value decomposition of F --> U*S*V'
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(F,Eigen::ComputeFullU | Eigen::ComputeFullV);
		S = svd.singularValues().asDiagonal(); // these will be returned in order (highest to lowest)
		U = svd.matrixU();
		Vt= svd.matrixV().transpose();
		P = 2*mu*(S-I) + lambda*(S-I).trace()*I; // P is now a diagonal matrix of principal stresses
		//printf("\n%% diag. deform: %.3le, %.3le, %.3le", S(0,0), S(1,1), S(2,2));
        //printf("\n%% stresses for element %d:",it->first);
        //printf("\n%% -- principal: %.3le, %.3le, %.3le", P(0,0), P(1,1), P(2,2));
		return 0;
	}




//******************************************************************************
//******************************************************************************
//******************************************************************************
// old stuff ...

    int PostProcessor::computeSIFs(
        const vector_type& displacements, vector_type& sifs,
        double E, double nu
    ){
        printf("\nOLD VERSION - DON'T USE (PostProcessor::computeSIFs)\n");
        sifs.resize(3*crackTips.size()); sifs.setZero();
        // loop through all crack-tip line segments and compute SIFs
        // using displacement correlation at the interior node
        // of the first parent triangle (boundary lines can have only 1 parent tri)
        int i=0;
        for(elem_map::const_iterator it = crackTips.begin();
            it!= crackTips.end(); ++it,++i
        ){
            unsigned int nd_a, nd_b, nd_c; // nodes of the current triangle
            Eigen::Vector3d a,b,c, n1,n2,n3;
            // nd_a and nd_b span the crack tip line segment
            nd_a = it->second[0]; nd_b = it->second[1];
            // nd_c is the interior node of the parent triangle
            nd_c = findInteriorNode(nd_a, nd_b, elems[parents[it->first][0]]);
            
            copyNode(nd_a, a);
            copyNode(nd_b, b);
            copyNode(nd_c, c);
            
            getLocalCoordFrame(a,b,c, n1,n2,n3);
            
            Eigen::Vector3d cod; // crack opening displacement at correlation point
            double r; // distance from the crack tip to the correlation point
            // correlation point is the interior node (nd_c)
            cod[0] = displacements[3*(nd_c-NODE_BASE_INDEX)  ];
            cod[1] = displacements[3*(nd_c-NODE_BASE_INDEX)+1];
            cod[2] = displacements[3*(nd_c-NODE_BASE_INDEX)+2];
            r = (((a+b)/2) -c).norm();
            double u1,u2,u3; // projected cod in n1,n2,n3 directions
            u1=cod.dot(n1);
            u2=cod.dot(n2);
            u3=cod.dot(n3);
			const double cf=1.0;//0.45; // u is linear over r, we can move the correlation point closer to the edge and u will fall linearly ~~> correction factor determined from penny shaped crack test case vs. theory in K_I
            //K1 = mu*sqrt(2*pi)*u1/(sqrt(r)*2*(1-nu));
            //K2 = mu*sqrt(2*pi)*u2/(sqrt(r)*2*(1-nu));
            //K3 = mu*sqrt(2*pi)*u3/sqrt(2*r);
            double mu = E/(2*(1+nu));
            sifs[3*i  ] =  -mu*sqrt(2*M_PI)*u1*cf/(sqrt(r*cf)*2*(1-nu)); //flip signs as opening displacement is opposite normals
            sifs[3*i+1] =  -mu*sqrt(2*M_PI)*u2*cf/(sqrt(r*cf)*2*(1-nu));
            sifs[3*i+2] =  -mu*sqrt(M_PI)*u3*cf/sqrt(2*r*cf);

        } // end for ( crackTip elements ), now i is the number of crack tip elements
		return 0;
    }
}
