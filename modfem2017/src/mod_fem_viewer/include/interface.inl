
	// Common mesh module interfaces
    int init_mesh(int Control, char * FileName)
	{
		return mmr_init_mesh(Control,FileName);
	}

	int free_mesh(int Mesh_id)
	{
		return mmr_free_mesh(Mesh_id);
	}

	int  get_nmno(int Mesh_id)
	{
		return mmr_get_max_node_id(Mesh_id);
	}

	int  get_nmed(int Mesh_id)
	{
		return mmr_get_max_edge_id(Mesh_id);
	}

	int  get_nmfa(int Mesh_id)
	{
		return mmr_get_max_face_id(Mesh_id);
	}

	int  get_nmel(int Mesh_id)
	{
		return mmr_get_max_elem_id(Mesh_id);
	}
	
	int  get_next_act_elem(int Mesh_id,int El)
	{
		return mmr_get_next_act_elem(Mesh_id,El);
	}

	int  get_next_elem(int Mesh_id,int El)
	{
		return mmr_get_next_elem_all(Mesh_id,El);
	}


	int  get_next_face(int Mesh_id,int Fa)
	{
		return mmr_get_next_face_all(Mesh_id,Fa);
	}

	int  get_next_edge(int Mesh_id,int Ed)
	{
		return mmr_get_next_edge_all(Mesh_id,Ed);
	}

	int  get_next_node(int Mesh_id,int Node)
	{
		return mmr_get_next_node_all(Mesh_id,Node);
	}

	int  get_el_status(int Mesh_id,int El)
	{
		return mmr_el_status(Mesh_id,El);
	}

	int  get_el_type(int Mesh_id,int El)
	{
		return mmr_el_type(Mesh_id,El);
	}

	int  get_el_mate(int Mesh_id,int El)
	{
		return mmr_el_mate(Mesh_id,El);
	}

	int  get_el_faces(int Mesh_id,int El,int* Faces, int* Orient)
	{
		return mmr_el_faces(Mesh_id,El,Faces,Orient);
	}

	int  get_el_fam(int Mesh_id,int El,int* Elsons,int* Type)
	{
		return mmr_el_fam(Mesh_id,El,Elsons,Type);
	}

	int  get_el_struct(int Mesh_id,int El,int *El_struct)
	{
		return mmr_elem_structure(Mesh_id,El,El_struct);
	}

	int  get_el_node_coor(int Mesh_id,int El,int* Nodes, double *Xcoor)
	{
		return mmr_el_node_coor(Mesh_id,El,Nodes,Xcoor);
	}
	
	int  get_fa_status(int Mesh_id,int Fa)
	{
		return mmr_fa_status(Mesh_id,Fa);
	}

	int  get_fa_type(int Mesh_id,int Fa)
	{
		return mmr_fa_type(Mesh_id,Fa);
	}

	int  get_fa_sub_bnd(int Mesh_id,int Fa)
	{
		return mmr_fa_sub_bnd(Mesh_id,Fa);
	}

	int  get_fa_node_coor(int Mesh_id,int Fa,int* Nodes,double* Coords)
	{
		return mmr_fa_node_coor(Mesh_id,Fa,Nodes,Coords);
	}

	void get_fa_neig(int Mesh_id, int Fa, int* Fa_neig,
			int* Neig_sides, int* Node_shift, int* Diff_gen,
			double* Acceff, double* Bcceff)
	{
		mmr_fa_neig(Mesh_id,Fa,Fa_neig,Neig_sides,Node_shift,Diff_gen,Acceff,Bcceff);
	}

	int  get_face_edges(int Mesh_id,int Fa,int *Edges,int *Orient)
	{
		return mmr_fa_edges(Mesh_id,Fa,Edges,Orient);
	}

	int  get_face_struct(int Mesh_id,int Fa,int *Fa_struct)
	{
		return mmr_face_structure(Mesh_id,Fa,Fa_struct);
	}

	int  get_edge_status(int Mesh_id,int Ed)
	{
		return mmr_edge_status(Mesh_id,Ed);
	}

	int  get_edge_struct(int Mesh_id,int Ed,int *Ed_struct)
	{
		return mmr_edge_structure(Mesh_id,Ed,Ed_struct);
	}

	int  get_edge_nodes(int Mesh_id,int Ed,int *Ed_nodes)
	{
		return mmr_edge_nodes(Mesh_id,Ed,Ed_nodes);
	}

	int  get_node_status(int Mesh_id,int Node)
	{
		 return mmr_node_status(Mesh_id,Node);
	}

	int  get_node_coor(int Mesh_id,int Node,double *Coor)
	{
		return mmr_node_coor(Mesh_id,Node,Coor);
	}

	// Common approximation module interfaces
	int  init_field(char Field_type,int Mesh_id,const char *Filename)
	{
		return apr_init_field(Field_type,APC_READ,Mesh_id,-1,-1,101,const_cast<char*>(Filename),NULL);
	}

	int free_field(int Field_id)
	{
		return apr_free_field(Field_id);
	}

	int  get_nreq(int Field_id)
	{
		return apr_get_nreq(Field_id);
	}

	int  get_nr_sol(int Field_id)
	{
		return apr_get_nr_sol(Field_id);
	}

	int  get_base_type(int Field_id)
	{
		return apr_get_base_type(Field_id);
	}
	
	int  get_el_pdeg(int Field_id,int El,int* Pdeg_vec)
	{
		return apr_get_el_pdeg(Field_id,El,Pdeg_vec);
	}

	int  get_element_dofs(int Field_id,int El_id,int Vec_id,double *Dofs)
	{
		return apr_get_el_dofs(Field_id,El_id,Vec_id,Dofs);
	}

