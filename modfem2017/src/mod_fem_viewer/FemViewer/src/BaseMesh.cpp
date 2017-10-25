/*
 * BaseMesh.cpp

 *
 *  Created on: 2011-04-22
 *      Author: Paweł Macioł
 */
#include <stdio.h>
#include "Log.h"
#include "Enums.h"
#include "fv_config.h"
#include "mmh_intf.h"
#include "BaseMesh.h"
#include "fv_dictstr.h"
#ifdef _USE_FV_EXT_MOD
#include "ApproxModule.h"
#endif

#ifndef _USE_FV_LIB
const int MMC_MOD_FEM_MESH_DATA = 0;
const int MMC_GRADMESH_DATA = 10;
const int MMC_MOD_FEM_PRISM_DATA = 1;
const int MMC_CUR_MESH_ID = 0;
#endif
#include "Enums.h"

// system
#include<string.h>
#include<sstream>
#include<iostream>

const char APPROXIMATION_DG  = 'D';
const char APPROXIMATION_STD = 'S';

namespace FemViewer {

#ifdef _USE_FV_LIB
intfint BaseMesh::Get_nmel = mmr_get_max_elem_id;
intfint BaseMesh::Get_nmno = mmr_get_max_node_id;
intf2int BaseMesh::Get_next_act_fa = mmr_get_next_act_face;
intf2int BaseMesh::Get_next_act_el = mmr_get_next_act_elem;
intf2int BaseMesh::Get_face_bc = mmr_fa_bc;
voidf2int4intp2doublep BaseMesh::Get_face_neig = mmr_fa_neig;
intf2int BaseMesh::Get_el_status = mmr_el_status;
intf2int2intp BaseMesh::Get_el_faces = mmr_el_faces;
intf2intintp BaseMesh::Get_el_struct = mmr_elem_structure;
intf2intintpdoublep BaseMesh::Get_el_node_coor = mmr_el_node_coor;
intf2intintp BaseMesh::Get_face_struct= mmr_face_structure;
intf2intintpdoublep BaseMesh::Get_fa_node_coor = mmr_fa_node_coor;
intfint BaseMesh::Get_nmed = mmr_get_max_edge_id;
intfint BaseMesh::Get_nmfa = mmr_get_max_face_id;
intf2intintp BaseMesh::Get_edge_nodes = mmr_edge_nodes;
intf2int BaseMesh::Get_edge_status = mmr_edge_status;
intf2intdoublep BaseMesh::Get_node_coor = mmr_node_coor;
intf2int BaseMesh::Get_node_status = mmr_node_status;
intfcharp BaseMesh::module_introduce = mmr_module_introduce;
intfintcharpfilep BaseMesh::init_mesh = mmr_init_mesh;
intfint BaseMesh::free_mesh = mmr_free_mesh;
#else
intfint BaseMesh::Get_nmel = NULL;
intfint BaseMesh::Get_nmno = NULL;
intf2int BaseMesh::Get_next_act_fa = NULL;
intf2int BaseMesh::Get_next_act_el = NULL;
intf2int BaseMesh::Get_face_bc = NULL;
voidf2int4intp2doublep BaseMesh::Get_face_neig = NULL;
intf2int BaseMesh::Get_el_status = NULL;
intf2int2intp BaseMesh::Get_el_faces = NULL;
intf2intintp BaseMesh::Get_el_struct = NULL;
intf2intintpdoublep BaseMesh::Get_el_node_coor = NULL;
intf2intintp BaseMesh::Get_face_struct= NULL;
intf2intintpdoublep BaseMesh::Get_fa_node_coor = NULL;
intfint BaseMesh::Get_nmed = NULL;
intfint BaseMesh::Get_nmfa = NULL;
intf2intintp BaseMesh::Get_edge_nodes = NULL;
intf2int BaseMesh::Get_edge_status = NULL;
intf2intdoublep BaseMesh::Get_node_coor = NULL;
intf2int BaseMesh::Get_node_status = NULL;
intfcharp BaseMesh::module_introduce = NULL;
intfintcharpfilep BaseMesh::init_mesh = NULL;
intfint BaseMesh::free_mesh = NULL;
#endif

int BaseMesh::_meshCounter = 0;
int BaseMesh::GetMeshModuleType()
{
	//mfp_debug("BaseMesh::GetMeshModuleType\n");
	char type[32];
	HandleType result = Unknown;
	module_introduce(type);
	if (strcmp(type,"3D_PRISM")) result = MeshPrizm;
	else if (strcmp(type,"3D_Hybrid")) result = MeshHybrid;
	else if (strcmp(type,"3D_remesh")) result = MeshRemesh;
	return result;
}

int BaseMesh::Init(const char* name_)
{
	#ifdef _USE_FV_LIB
	this->idx() = MMC_CUR_MESH_ID;
	#else
	char name[12];
	module_introduce(name);
	int control = can_read(name,name_.c_str());
	this->idx() = init_mesh(control,const_cast<char*>(name_);
	if (this->idx() < 0) return(-1);
	#endif
	this->type() = static_cast<HandleType>(GetMeshModuleType());
	_meshCounter++;
	_name = name_ ? std::string(name_) : "";
	return this->idx();
}

int BaseMesh::Free()
{
	#ifndef _USE_FV_LIB
	free_mesh(MMC_CUR_MESH_ID);
	this->idx() = -1;
	#endif
	_meshCounter--;
	_name.clear();
	return(0);
}


} // end namespace FemViewer

