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

	/** // node based version
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
			bool flag=false;
			double mu=E/(2*(1+nu)), lambda=E*nu/((1+nu)*(1-2*nu)); // convert (E, nu) to Lame parameters (mu = shear modulus, lambda) ; as usual expect trouble for nu=0.5
			Eigen::Vector3d // a,b,c are node coordinates in material space; ua,ub,uc are nodal displacements, so world coordinates are a+ua etc.
				a (nodes[it->second[0]][0], nodes[it->second[0]][1], nodes[it->second[0]][2]),
				b (nodes[it->second[1]][0], nodes[it->second[1]][1], nodes[it->second[1]][2]),
				c (nodes[it->second[2]][0], nodes[it->second[2]][1], nodes[it->second[2]][2]),
				ua(u[3*(it->second[0]-NODE_BASE_INDEX)], u[3*(it->second[0]-NODE_BASE_INDEX)+1], u[3*(it->second[0]-NODE_BASE_INDEX)+2]),
				ub(u[3*(it->second[1]-NODE_BASE_INDEX)], u[3*(it->second[1]-NODE_BASE_INDEX)+1], u[3*(it->second[1]-NODE_BASE_INDEX)+2]),
				uc(u[3*(it->second[2]-NODE_BASE_INDEX)], u[3*(it->second[2]-NODE_BASE_INDEX)+1], u[3*(it->second[2]-NODE_BASE_INDEX)+2]);
            if(cracks.count(regions[it->first])!=0){ // this is a crack element
				Eigen::Vector3d // uba,ubb,ubc are nodal crack base displacements
					uba(u_c[3*(it->second[0]-NODE_BASE_INDEX)], u_c[3*(it->second[0]-NODE_BASE_INDEX)+1], u_c[3*(it->second[0]-NODE_BASE_INDEX)+2]),
					ubb(u_c[3*(it->second[1]-NODE_BASE_INDEX)], u_c[3*(it->second[1]-NODE_BASE_INDEX)+1], u_c[3*(it->second[1]-NODE_BASE_INDEX)+2]),
					ubc(u_c[3*(it->second[2]-NODE_BASE_INDEX)], u_c[3*(it->second[2]-NODE_BASE_INDEX)+1], u_c[3*(it->second[2]-NODE_BASE_INDEX)+2]);
                ua*=0.5; ub*=0.5; uc*=0.5; // use half of crack opening displacement (because COD is distributed over 2 surfaces)
				ua+=uba; ub+=ubb; uc+=ubc; // add the base displacement
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
            maxPrincipalStress[it->first].assign(8,0.0);
            maxPrincipalStress[it->first][0]=P(0,0); // first entry is the stress value
			// store max. principal stress and its plane-normal vector for each element (also min. for compressive fracture?)
            // in the SVD, V is the pre-rotation and U is the post-rotation
            // so the direction of max. principal stress is the first column of V
            // which becomes the first row of V.transpose()
            maxPrincipalStress[it->first][1]=Vt(0,0); // next 3 entries
            maxPrincipalStress[it->first][2]=Vt(0,1); // are the plane-
            maxPrincipalStress[it->first][3]=Vt(0,2); // normal vector
            // same for min. principal stress used for compressive fractures
			maxPrincipalStress[it->first][4]=P(2,2);
			maxPrincipalStress[it->first][5]=Vt(2,0); // next 3 entries
            maxPrincipalStress[it->first][6]=Vt(2,1); // are the plane-
            maxPrincipalStress[it->first][7]=Vt(2,2); // normal vector
            //s_xx[i]=P(0,0); s_yy[i]=P(1,1); s_zz[i]=P(2,2); // for testing
			
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
    
    //no longer supported: PostProcessor::propagateCrackOLD -- use SubsampledCrackTip class instead
//    int PostProcessor::propagateCrackOLD(
//        node_map& newNodes, elem_map& newElems, id_map& newRegions,
//        elem_map& newCrackTips, elem_map& newParents,
//        state_map& newCrackTipStates, const vector_type& sifs,
//		const FractureModel& fractureModel, bool isSubsteppedFractureModel
//    ){
//		printf("\nOLD VERSION - DON'T USE (PostProcessor::propagateCrackOLD)\n");
//        newNodes.clear(); newElems.clear(); newRegions.clear();
//        newCrackTips.clear(); newParents.clear();
//        newCrackTipStates.clear();
//        // 1. for each crack tip element compute new positions of both endpoints
//        // 2. store old node ID of endpoints
//        // 3. average positions of nodes belonging to multiple elements
//        // 4. assign new node IDs
//        // 5. triangulate area between old and new crack tip (make sure ordering is correct!)
//        // 5. build output maps (regions are copied from first parent of old crack tip)
//        
//        // storage for new positions (both endpoints of all crack tip elements)
//        unsigned int i; // counter var
//		int n = crackTips.size();
//        vector<Eigen::Vector3d> newPos(2*n,Eigen::Vector3d());
//        vector<bool> hasNewPos(2*n,false);
//        vector<int> oldNodeIDs(2*n,0);
//        Eigen::Vector3d dn; // store result of fracture criterion
//
//		if(isSubsteppedFractureModel){
//			//discarded, use SubsampledCrackTip class!		
//		}else{ // older version for simple fracture models
//			i=0;
//			for(elem_map::const_iterator it = crackTips.begin();
//				it!= crackTips.end(); ++it,++i
//			){
//				// sifs on current edge are in sifs[3*i..]
//				// node IDs of endpoints are it->second[0] and it->second[1]
//				// write new positions to newPos[2*i..]
//				unsigned int nd_a, nd_b, nd_c; // nodes of the current triangle
//				nd_a = it->second[0]; nd_b = it->second[1];
//				// store old node IDs
//				oldNodeIDs[2*i  ] = nd_a;
//				oldNodeIDs[2*i+1] = nd_b;
//				// copy old positions to newPos
//				copyNode(nd_a, newPos[2*i  ]);
//				copyNode(nd_b, newPos[2*i+1]);
//				// check propagation criterion
//				Eigen::Vector3d edgeSIFs(sifs[3*i], sifs[3*i+1], sifs[3*i+2]);
//				if( crackTipStates[it->first]!=INACTIVE &&
//					fractureModel.fractureCriterion(dn, edgeSIFs)
//				){ // propagation criterion met
//					hasNewPos[2*i  ] = true;
//					hasNewPos[2*i+1] = true;
//                                
//					// build the local coordinate frame
//					Eigen::Vector3d a,b,c, n1,n2,n3;
//					// nd_a and nd_b span the crack tip line segment
//					nd_a = it->second[0]; nd_b = it->second[1];
//					// nd_c is the interior node of the parent triangle
//					nd_c = findInteriorNode(nd_a, nd_b, elems[parents[it->first][0]]);
//					copyNode(nd_a, a);
//					copyNode(nd_b, b);
//					copyNode(nd_c, c);
//					getLocalCoordFrame(a,b,c, n1,n2,n3);
//					Eigen::Vector3d dx1, dx2; // transform dn to world coordinate displacements of endpoints
//					// dn[0] is radius of displacement, dn[1]=theta is angle of rotation around the segment's axis
//					// allow the crack segment to "twist" around the forward edge normal (n2) by angle phi=dn[2]
//					dx1 = n2*dn[0]*cos(dn[1]) + n1*dn[0]*sin(dn[1]);
//					double halfL=0.5*(b-a).norm();
//					dx2 = n3*halfL*cos(dn[2]) - n1*halfL*sin(dn[2]) ;
//					// update positions of both endpoints
//					newPos[2*i  ] = 0.5*(a+b)+dx1-dx2;
//					newPos[2*i+1] = 0.5*(a+b)+dx1+dx2;
//				}else{ // propagation criterion not met
//					hasNewPos[2*i  ] = false;
//					hasNewPos[2*i+1] = false;
//				}            
//			} // end for ( crack tip elements )
//		}
//        // now we have a bunch of new positions for the element endpoints
//        // average them together to form the new nodes
//        // and reconnect to form the new crack tip elements
//
//        id_map newNodeIDs, nodeCount; //newNodeIDs[oldNodeID] = newNodeID
//		id_map edgeCount;//, newNodeIDsBckwd; // newNodeIDsBckwd[newNodeID] = oldNodeID
//        unsigned int nextID=(nodes.rbegin()->first)+1; // take the last (highest) node ID and add 1
//        // this is a "smoothing" mechanism on the crack tip (enables averaging over not-propagated nodes)
//        for(int j=0; j<2*n; ++j){
//            if(hasNewPos[j]){ // count how many edges adjacent to a node have propagated
//                if(edgeCount.count(oldNodeIDs[j])==0) // no entry exists yet
//                    edgeCount[oldNodeIDs[j]]=1;       // create it here
//                else
//                    edgeCount[oldNodeIDs[j]]++;       // add to existing entry
//            }
//        }
//        // use hasNewPos[j] instead of edgeCount[...] for less smooth results
//        for(int j=0; j<2*n; ++j){ // iterate through newPos
//            if(/**hasNewPos[j] /*/edgeCount[oldNodeIDs[j]]>0/**/){ // node has propagated from at least one side
//                // assign new node ID as newNodeIDs[oldNodeIDs[j]] = nextNodeID
//                // add newPos[j] to newNodes[newNodeIDs[oldNodeIDs[j]]] and count
//                vector<double> pos(newPos[j].data(),newPos[j].data()+3);
//                if(newNodeIDs.count(oldNodeIDs[j])==0){ // assign new number at first occurence of node
//                    i=nextID++; //post-increment!
//                    newNodeIDs[oldNodeIDs[j]]= i;
//					//newNodeIDsBckwd[i] = oldNodeIDs[j];
//                    newNodes[i] = pos;
//                    nodeCount[i] = 1;
//                }else{ // node already known, add up
//                    i=newNodeIDs[oldNodeIDs[j]];
//                    newNodes[i][0] += pos[0];
//                    newNodes[i][1] += pos[1];
//                    newNodes[i][2] += pos[2];
//                    ++nodeCount[i];
//                }
//			}//else{} // nothing to do here ...
//        }
//        
//		//vector<int> reject;
//        for(node_map::iterator it=newNodes.begin(); it!=newNodes.end(); ++it){
//			// divide by count
//            it->second[0] /= nodeCount[it->first];
//            it->second[1] /= nodeCount[it->first];
//            it->second[2] /= nodeCount[it->first];
//			// averaging done
//			//printf("nd %d count %d\n",it->first, nodeCount[it->first]);
//
//			// check positions
//			//if(levelSet.isInside(it->second[0],it->second[1],it->second[2])){
//			//	// position is inside ...
//			//
//			//}else{ // position is outside of the object
//			//	// old position (original node)
//			//	vector<double>& oldPos = nodes[newNodeIDsBckwd[it->first]];
//			//	Eigen::Vector3d q(oldPos[0], oldPos[1], oldPos[2]);
//			//	// new position (which is invalid) 
//			//	Eigen::Vector3d p(it->second[0],it->second[1],it->second[2]);
//
//			//	int tri = levelSet.getClosestTriangleToLine(p,q); // find triangle closest to line pq
//			//	printf("new position is outside (%.2lf,%.2lf,%.2lf), closest triangle %d\n", it->second[0],it->second[1],it->second[2], tri);
//
//			//	Eigen::Vector3d v = (q-p); v.normalize(); // v is unit vector along the line
//			//	Eigen::Vector3d a,b,c, n;
//			//	copyNode(elems[tri][0], a);
//			//	copyNode(elems[tri][1], b);
//			//	copyNode(elems[tri][2], c);
//			//	n=(b-a).cross(c-a); n.normalize(); // n is triangle normal unit vector
//
//			//	// intersect line with plane of triangle
//			//	double d = n.dot(v), s=0.0;
//			//	if(std::abs(d) > DBL_EPSILON) s = n.dot(a-p) / d;
//
//			//	if(s<DBL_EPSILON){ // consider rejecting propagation?
//			//	}else{
//			//		// new position
//			//		q = p+s*v;
//			//		it->second[0]=q[0];
//			//		it->second[1]=q[1];
//			//		it->second[2]=q[2];
//			//	}
//			//}
//        }
//		//for(i=0; i<reject.size(); ++i){
//		//	newNodes.erase(i);
//		//	newNodeIDs.erase(i);
//		//}
//        
////        //debug output
////        printf("%%new nodes from crack propagation\nnew_nodes=[\n");
////        for(node_map::iterator it=newNodes.begin(); it!=newNodes.end(); ++it){
////            printf(" %.6le %.6le %.6le %% %d\n",
////                it->second[0], it->second[1], it->second[2], it->first);
////        } printf("];\n");
//        
//        // now newNodes has all the updated positions of the crack tip nodes
//        // they are already known by their new index
//        // mapping of old to new index is in newNodeIDs
//        
//		//printf("** building new mesh\n");
//        // here we check the "state" of a node, based on it's adjacent elements
//        // the point is to connect the new crack tip to an ACTIVE node if possible
//        state_map nodeStates; //sum old crackTipStates to nodes, using ACTIVE, ACTIVE_A, INACTIVE
//        for(elem_map::iterator it=crackTips.begin();
//            it!=crackTips.end(); ++it){
//            for(int k=0; k<=1; ++k){ // for both nodes of the segment
//                if(nodeStates.count(it->second[k])==0){ //assign
//                    nodeStates[it->second[k]]=
//                            (crackTipStates[it->first]==INACTIVE)?ACTIVE_A:ACTIVE;
//                    //printf("\n%% ... init node %d to state %d",it->second[k],nodeStates[it->second[k]]);
//                }else{ //merge
//                    if(crackTipStates[it->first]==INACTIVE)
//                        nodeStates[it->second[k]]=
//                            (nodeStates[it->second[k]]==ACTIVE)?ACTIVE_A:INACTIVE;
//                    //printf("\n%% ... merged node %d to state %d",it->second[k],nodeStates[it->second[k]]);
//                }
//            }
//        }//now nodes connecting active segments are ACTIVE
//        // and nodes connecting one inactive segment are ACTIVE_A
//        // while nodes connecting more than one (ie. two) inactive segments are INACTIVE
////        for(state_map::iterator it=nodeStates.begin(); it!=nodeStates.end();++it)
////            printf("\n%% node %d has state %d ... ",it->first,it->second);
//        
//        // loop over the old crack tip, check whether it has propagated and triangulate
//        nextID=(elems.rbegin()->first)+1; // take the last (highest) element ID and add 1
//        i=0;
//        for(elem_map::const_iterator it = crackTips.begin();
//            it!= crackTips.end(); ++it,++i
//        ){
//            unsigned int nd_a, nd_b;
//            nd_a = it->second[0]; nd_b = it->second[1];
//            int newA=-1, newB=-1;
//            if( newNodeIDs.count(nd_a) >0 )
//                    newA=newNodeIDs[nd_a]; // only access map if element exists, otherwise it default-constructs!
//            if( newNodeIDs.count(nd_b) >0 )
//                    newB=newNodeIDs[nd_b];
//            
//            // build new crack tip element
//            newCrackTips[i+ELEM_BASE_INDEX].push_back( (newA>0)?newA:(nd_a) );
//            newCrackTips[i+ELEM_BASE_INDEX].push_back( (newB>0)?newB:(nd_b) );
//            newCrackTips[i+ELEM_BASE_INDEX].push_back(nd_a);
//            newCrackTips[i+ELEM_BASE_INDEX].push_back(nd_b);
//            
//            // build triangulation behind crack tip
//            if(newA <0 && newB >0){ // a not propagated, only "right" triangle
//                newElems[nextID].push_back(nd_a);
//                newElems[nextID].push_back(newB);
//                newElems[nextID].push_back(nd_b); // b is interior node in tri
//                newRegions[nextID]=regions[parents[it->first][0]];
//                newParents[i+ELEM_BASE_INDEX].push_back(nextID);
//                newParents[i+ELEM_BASE_INDEX].push_back(0); // keep Elmer format intact (each edge can have two parents)
//                if( crackTipStates[it->first]==ACTIVE )
//					newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE;
//				else if( crackTipStates[it->first]==ACTIVE_A )
//                    newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE_A;
//				else if( crackTipStates[it->first]==ACTIVE_B )
//                    newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE_B;
//                else
//					newCrackTipStates[i+ELEM_BASE_INDEX]=crackTipStates[it->first];
//                ++nextID;
//            }else if(newA >0 && newB <0){ // b not propagated, only "left" triangle
//                newElems[nextID].push_back(newA);
//                newElems[nextID].push_back(nd_b);
//                newElems[nextID].push_back(nd_a); // a is interior node in tri
//                newRegions[nextID]=regions[parents[it->first][0]];
//                newParents[i+ELEM_BASE_INDEX].push_back(nextID);
//                newParents[i+ELEM_BASE_INDEX].push_back(0); // keep Elmer format intact (each edge can have two parents)
//                if( crackTipStates[it->first]==ACTIVE )
//					newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE;
//				else if( crackTipStates[it->first]==ACTIVE_A )
//                    newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE_A;
//				else if( crackTipStates[it->first]==ACTIVE_B )
//                    newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE_B;
//                else
//					newCrackTipStates[i+ELEM_BASE_INDEX]=crackTipStates[it->first];
//                ++nextID;
//            }else if(newA >0 && newB >0){ // both triangles
//                unsigned int back_node=nd_a, front_node=newB;
//                if(nodeStates[back_node]!=ACTIVE){ //flip diagonal of quad
//                    back_node=nd_b; front_node=newA;
//                }
//                newElems[nextID].push_back(newA);
//                newElems[nextID].push_back(newB);
//                newElems[nextID].push_back(back_node);
//                newRegions[nextID]=regions[parents[it->first][0]];
//                newParents[i+ELEM_BASE_INDEX].push_back(nextID);
//                newParents[i+ELEM_BASE_INDEX].push_back(0); // keep Elmer format intact (each edge can have two parents)
//                if( crackTipStates[it->first]==ACTIVE )
//					newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE;
//				else if( crackTipStates[it->first]==ACTIVE_A )
//                    newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE_A;
//				else if( crackTipStates[it->first]==ACTIVE_B )
//                    newCrackTipStates[i+ELEM_BASE_INDEX]=UPDATE_B;
//                else
//					newCrackTipStates[i+ELEM_BASE_INDEX]=crackTipStates[it->first];
//                ++nextID;
//                newElems[nextID].push_back(nd_a);
//                newElems[nextID].push_back(front_node);
//                newElems[nextID].push_back(nd_b);
//                newRegions[nextID]=regions[parents[it->first][0]];
//                ++nextID;
//			}else{
//                // neither node propagated, keep geometry unchanged
//                // important: copy parents and state
//                newParents[i+ELEM_BASE_INDEX]=parents[it->first];
//                newCrackTipStates[i+ELEM_BASE_INDEX]=crackTipStates[it->first];
//            }
//        } // all done
//        //printf("** done\n");
////        printf("new_elems=[\n");
////        for(elem_map::iterator it = newElems.begin(); it!= newElems.end(); ++it){
////            printf(" %d %d %d %% %d\n", it->second[0], it->second[1], it->second[2], it->first);
////        } printf("];\n");
////        printf("new_ct=[\n");
////        for(elem_map::iterator it = newCrackTip.begin(); it!= newCrackTip.end(); ++it){
////            printf(" %d %d\n", it->second[0], it->second[1]);
////        } printf("];\n");
////        printf("new_ct_parents=[\n");
////        for(elem_map::iterator it = newParents.begin(); it!= newParents.end(); ++it){
////            printf(" %d %d (%d)\n", it->first, it->second[0], it->second[1]);
////        } printf("];\n");
//        
//        //delete[] newPos;
//        //delete[] oldNodeIDs;
//        //delete[] hasNewPos;
//        
//		return 0;
//    }
}
