#ifndef GHOST_BLOCK_EXCHANGER_HPP
#define GHOST_BLOCK_EXCHANGER_HPP

#include<stdlib.h>
#include<stdio.h>

#ifdef PARALLEL
#include "mmh_intf.h"
#include "mmph_intf.h"
#include "aph_intf.h"
#include <mpi.h>

/* interface for parallel communication modules */
#include "pch_intf.h"
#endif

void init_posglobs(int* posglobs, int posglob_offset);

#ifdef PARALLEL

int exchange_global_row_indices(
                   /* returns: >=0 -success code, <0 -error code */
  int Nr_dof_ent,  /* in: number of DOF entities in the level */
  int* L_dof_ent_id,  /* in: list of DOF entities associated with DOF structs */
  int* L_struct_nrdof,    /* in: list of nrdofs for each dof struct */
  int* L_struct_posg,     /* in: list of positions within the global */
                      /*     vector of dofs for each dof struct */
  int* L_thing_to_struct, /*in: list of DOF structs associated with DOF entities */
  int field_id
  );

#endif
#endif
