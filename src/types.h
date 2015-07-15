#ifndef TYPES_H
#define	TYPES_H

#include <map>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include <cstdio>

namespace FractureSim{
	enum CONSTANTS{
        NODE_BASE_INDEX = 1, // index origin of nodes
        ELEM_BASE_INDEX = 1
    };
    enum CRACK_STATE{ // state of crack-tips
        /*0*/ACTIVE,   // can propagate, if fracture criterion matches
		/*1*/ACTIVE_A, // only node a is valid, b has intersected boundary
		/*2*/ACTIVE_B, // only node b is valid, a has intersected boundary
        /*3*/INACTIVE, // not allowed to propagate for any reason
        /*4*/UPDATE,   // modified ("dirty" state)
        /*5*/UPDATE_A, // modified, was ACTIVE_A before
        /*6*/UPDATE_B  // modified, was ACTIVE_B before
    };

	typedef std::map<unsigned int, std::vector<double> > node_map;
	typedef std::map<unsigned int, std::vector<unsigned int> > elem_map;
	typedef std::map<unsigned int, unsigned int> id_map;
    typedef std::set<unsigned int> id_set;
    typedef std::map<unsigned int, CRACK_STATE> state_map;
	typedef std::map<unsigned int, Eigen::Vector3d> vect3d_map;
    
	typedef Eigen::MatrixXd matrix_type;
	typedef Eigen::VectorXd vector_type;
    
    inline id_set nodeSet(elem_map elems){
        id_set nodes;
        for(elem_map::iterator i = elems.begin(); i != elems.end(); ++i)
            for(elem_map::mapped_type::iterator j = i->second.begin();
                j != i->second.end(); ++j)
                nodes.insert(*j);
        return nodes;
    }

	// NEW FOR FractureRB:
	typedef std::pair<unsigned int, unsigned int> edge;
	typedef std::map<edge, unsigned int> edge_imap;

    inline id_set nodeSet(edge_imap edges){
        id_set nodes;
        for(edge_imap::iterator i = edges.begin(); i != edges.end(); ++i){
			nodes.insert(i->first.first);
			nodes.insert(i->first.second);
		}
        return nodes;
    }
}

#endif
