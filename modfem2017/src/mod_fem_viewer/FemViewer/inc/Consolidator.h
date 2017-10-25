#ifndef _CONSOLIDATROR_H__
#define _CONSOLIDATROR_H__

/* for mesh module */
#include "mmh_intf.h"

/* for approximation module */
#include "aph_intf.h"

class Consolidator {

private:
	Consolidator() {}
	friend Consolidator& ConsolInstance();

public:
	/* returns maximal element id */
	inline int get_elems_count(int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_get_max_elem_id(Mesh_id);
	}

	inline int get_next_act_elem(int Mesh_id = MMC_CUR_MESH_ID, int Nel = 0) {
		return mmr_get_next_act_elem(Mesh_id,Nel);
	}

	inline int get_elem_status(int Mesh_id = MMC_CUR_MESH_ID, int Nel = 0) {
		return mmr_el_status(Mesh_id,Nel);
	}

	inline int get_elem_faces(int Nel, int *Faces, int *Orient, int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_el_faces(Mesh_id,Nel,Faces,Orient);
	}

	inline int get_elem_struct(int *El_str,int El,int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_elem_structure(Mesh_id,El,El_str);
	}

	inline int get_elem_node_coor(double *Coords,int *Nodes,int El,int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_el_node_coor(Mesh_id,El,Nodes,Coords);
	}

	inline int get_face_struct(int Fa,int *Face_str, int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_face_structure(Mesh_id,Fa,Face_str);
	}

	inline int get_face_node_coor(int Fa,int *Nodes,double *Coords,int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_fa_node_coor(Mesh_id,Fa,Nodes,Coords);
	}

	inline int get_max_edge_id(int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_get_max_edge_id(Mesh_id);
	}

	inline int get_edge_status(int Ed,int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_edge_status(Mesh_id,Ed);
	}

	inline int get_edge_nodes(int *Nodes,int Ed,int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_edge_nodes(Mesh_id,Ed,Nodes);
	}

	inline int get_node_coor(double* Coords,int No, int Mesh_id = MMC_CUR_MESH_ID) {
		return mmr_node_coor(Mesh_id,No,Coords);
	}

	inline int get_elem_pdeg() {
		return(0);
	}

};

extern Consolidator* pConsolidator;
Consolidator& ConsolInstance();


#endif /* _CONSOLIDATROR_H__ */
