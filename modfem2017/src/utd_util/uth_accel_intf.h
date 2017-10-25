/***********************************************************************
File uth_accel_intf - utility routines supporting streaming processing 
                      on accelerators (common to all problem modules)

Contains declarations of routines:

local:
  utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
                           and assemble them to the global SM using ACCELERATOR
------------------------------
History:
	08.2015 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
*************************************************************************/

#ifndef _uth_accel_intf_
#define _uth_accel_intf_


/*------------------------------------------------------------
 utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
                                 and assemble them to the global SM using ACCELERATOR
------------------------------------------------------------*/
int utr_create_assemble_stiff_mat_elem_accel(
  int Problem_id, 
  int Level_id, 
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int* Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
  int Nr_int_ent,
  int* L_int_ent_type,
  int* L_int_ent_id,
  int Nr_colors_elems,
  int* L_color_index_elems, 
  int Nr_colors_accel,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int Max_dofs_int_ent
					 );


#endif
