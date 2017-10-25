/************************************************************************
File pds_conv_diff_dg_util.c - approximation dependent utility procedures
                               for convection-diffusion equations

Contains routines:
  pdr_limit_slope - perform slope limiting
 ------------------------------  			
History:    
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

/* header files for the problem dependent module for Conv_Diff's equation */
#include "../../include/pdh_conv_diff.h"	

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include "pdh_control_intf.h"

#include "pdh_intf.h"

/* interface of the mesh manipulation module */
#include "mmh_intf.h"	

/* interface for all approximation modules */
#include "aph_intf.h"	

/* utilities - including simple time measurement library */
#include "uth_intf.h"

/* magical parameters for slope limiting */
#define MAGIC1 1.0e-12
#define SMALL  1.0e-10


/*---------------------------------------------------------
pdr_slope_limit - to limit the slope of linear solution
---------------------------------------------------------*/
int pdr_slope_limit( /* returns: >0 - success, <=0 - failure */
	int Problem_id	/* in: data structure to be used  */
	)
{


  int nmel;		/* number of elements */
  int el_type;		/* element type */
  int nreq;		/* number of equations */
  double time;          /* time instant */
  int num_shap;         /* number of element shape functions */
  int el_faces[MMC_MAXELFAC+1]; /* element faces */
  int el_neig[MMC_MAXELFAC+1]; /* equal or larger element's neighbors */
  double node_coor[3*MMC_MAXELVNO];  /* coord of nodes */
  double dofs_loc[APC_MAXELSD]; /* element solution dofs */
  double node_coor1[3*MMC_MAXELVNO];  /* coord of nodes */
  double dofs_loc1[APC_MAXELSD]; /* element solution dofs */

  double hsize, zsize, xysize;          /* element sizes */
  double uc0,uc1,uc2,uc3,xc0[3],xc1[3],xc2[3],xc3[3];
  double daux,ur1,ur2,ur3;
  double utilda_x,utilda_y, utilda_z;

/* auxiliary variables */
  int mesh_id, field_id, el_mate;
  int i, iaux, iel, boundary, neig, num_nodes;
  int bc_type1, bc_type2, bc_type3, slope;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* get formulation parameters */
  //i=1; name=pdr_ctrl_i_params(Problem_id,i); changed for string
  i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  i=3; field_id=pdr_ctrl_i_params(Problem_id,i);
  i=5; nreq=pdr_ctrl_i_params(Problem_id,i);
  i=4; time=pdr_time_d_params(Problem_id,i);
  i=11; slope=pdr_ctrl_i_params(Problem_id,i);
  double coeff3=1.0; // limiting strength parameter

  if (slope == 1) { /* initial technique based on reference planes */

    if(nreq>1) return(0);

/* for each active element */
    iel=0;
    while((iel=mmr_get_next_act_elem(mesh_id,iel))!=0){

/* find element type */
	el_type = mmr_el_type(mesh_id,iel);
	el_mate = mmr_el_groupID(mesh_id,iel);

/* find degree of polynomial and number of element scalar dofs */
	num_shap = apr_get_ent_numshap(field_id, APC_ELEMENT, iel);

/* find faces and equal size neighbors */
	mmr_el_faces(mesh_id, iel, el_faces, NULL);
	mmr_el_eq_neig(mesh_id, iel, el_neig, NULL);

/* get the coordinates of element nodes in the right order */
	num_nodes = mmr_el_node_coor(mesh_id, iel, NULL, node_coor);

/* get  element size */
	hsize=mmr_el_hsize(mesh_id, iel, NULL,NULL,NULL);

/* get current solution degrees of freedom */
	i=1;
	apr_read_ent_dofs(field_id,APC_ELEMENT,iel,num_shap,i,dofs_loc);

/*kbw
if(iel>0&&put_vec3_length(&dofs_loc[1])*fabs(dofs_loc[0])>SMALL){
printf("In SLOPE_LIMIT: element %d, type %d, hsize %lf\n",
			iel,el_type,num_shap,hsize);
printf("Neighbors:");
for(i=1;i<=el_neig[0];i++) printf(" %d",el_neig[i]);
printf("\nNodes' coordinates and solution:\n");
for(i=0;i<num_nodes;i++) printf("%.9e %.9e %.9e \n",
		 node_coor[3*i],node_coor[3*i+1],node_coor[3*i+2]);
for(i=0;i<4;i++) printf("%20.15lf",dofs_loc[i]);
printf("\n");
getchar();
}
/*kew*/

/* skip elements on inter-generation boundary */
	boundary=0;
	for(i=1;i<=el_neig[0];i++){
/* if larger neighbor or inactive neighbor */
	  if(el_neig[i]<0||
	     (el_neig[i]>0 && mmr_el_status(mesh_id,el_neig[i])<=0)) {
	    boundary=1;
	    break;
	  }
	}

/* skip if piecewise constant elements or inter-generation boundary */
	if(num_shap==1 || boundary==1){
	  continue;
	}

/* for prismatic element */
	if(el_type == MMC_PRISM) {

/* ad hoc sizes */
	  zsize=fabs(node_coor[2]-node_coor[11]);
	  xysize=sqrt(hsize*hsize*hsize/zsize);

/* center coordinates */
	  xc0[0]=(node_coor[0]+node_coor[3]+node_coor[6]+
		  node_coor[9]+node_coor[12]+node_coor[15])/6.0;
	  xc0[1]=(node_coor[1]+node_coor[4]+node_coor[7]+
		  node_coor[10]+node_coor[13]+node_coor[16])/6.0;
	  xc0[2]=(node_coor[2]+node_coor[5]+node_coor[8]+
		  node_coor[11]+node_coor[14]+node_coor[17])/6.0;

/* average value (for p=1 value at geometrical center) */
	  if(num_shap==3||num_shap==4){
	    uc0=dofs_loc[0]+dofs_loc[1]/3.0+dofs_loc[2]/3.0;
	  }
	  else{
	    uc0=pdr_average_sol_el(Problem_id,iel);
	  }

#ifdef DEBUG
	  if(num_shap==4 && fabs(uc0-pdr_average_sol_el(Problem_id,iel)) > SMALL){
	    printf("Error in element %d average solution %lf != %lf\n",
		   iel,uc0,pdr_average_sol_el(Problem_id,iel));
	  }
#endif

/*kbw
daux=pdr_average_sol_el(Problem_id,iel);
printf("Element %d - center %lf %lf %lf, value %lf (average %lf)\n",
iel,xc0[0],xc0[1],xc0[2],uc0,daux);
/*kew*/

/* limit in z direction if and only if */
	  if(num_shap>3 && fabs(dofs_loc[3])>SMALL) {

/* perform limiting in z-direction */

/* get neighbor across first base */
	    neig=el_neig[1];

/* get bc type of first (bottom) face */
	    bc_type1=pdr_get_bc_type(mmr_fa_bc(mesh_id,el_faces[1]));

/* if not on the boundary */
	    if(neig>0){

/* get the coordinates of element nodes in the right order */
	      mmr_el_node_coor(mesh_id, neig, NULL, node_coor1);

/* center coordinates */
	      xc1[0]=(node_coor1[0]+node_coor1[3]+node_coor1[6]+
		      node_coor1[9]+node_coor1[12]+node_coor1[15])/6.0;
	      xc1[1]=(node_coor1[1]+node_coor1[4]+node_coor1[7]+
		      node_coor1[10]+node_coor1[13]+node_coor1[16])/6.0;
	      xc1[2]=(node_coor1[2]+node_coor1[5]+node_coor1[8]+
		      node_coor1[11]+node_coor1[14]+node_coor1[17])/6.0;

/* average value (for p=1 value at geometrical center) */
	      iaux = apr_get_ent_pdeg(field_id, APC_ELEMENT, neig);
	      if(iaux==1){
/* get solution degrees of freedom */
		i=1;
		apr_read_ent_dofs(field_id,APC_ELEMENT,neig,num_shap,i,
				  dofs_loc1);

		uc1=dofs_loc1[0]+dofs_loc1[1]/3.0+dofs_loc1[2]/3.0;
	      }
	      else{
		uc1=pdr_average_sol_el(Problem_id,neig);
	      }

	    } /* end if not on the boundary */
	    else{
	      if (bc_type1==PDC_BC_DIRI) {
	        pdr_bc_diri_coeff(el_faces[1], el_mate, xc1, NULL,
				  time, NULL,NULL,NULL,NULL,&uc1);
	      }	      
	      else if (bc_type1==PDC_BC_MIXED) {
	        pdr_bc_mixed_coeff(el_faces[1], el_mate, xc1, NULL,
				     time,NULL,NULL,NULL,NULL,&uc1,&daux);
	      }	      
	    }

/* get neighbor across second base */
	    neig=el_neig[2];

/* get bc type of first (bottom) face */
	    bc_type2=pdr_get_bc_type(mmr_fa_bc(mesh_id,el_faces[2]));

/* if not on the boundary */
	    if(neig>0){

/* get the coordinates of element nodes in the right order */
	      mmr_el_node_coor(mesh_id, neig, NULL, node_coor1);

/* center coordinates */
	      xc2[0]=(node_coor1[0]+node_coor1[3]+node_coor1[6]+
		      node_coor1[9]+node_coor1[12]+node_coor1[15])/6.0;
	      xc2[1]=(node_coor1[1]+node_coor1[4]+node_coor1[7]+
		      node_coor1[10]+node_coor1[13]+node_coor1[16])/6.0;
	      xc2[2]=(node_coor1[2]+node_coor1[5]+node_coor1[8]+
		      node_coor1[11]+node_coor1[14]+node_coor1[17])/6.0;

/* average value (for p=1 value at geometrical center) */
	      iaux=apr_get_ent_pdeg(field_id, APC_ELEMENT, neig);
	      if(iaux==1){
/* get solution degrees of freedom */
		i=1;
		apr_read_ent_dofs(field_id,APC_ELEMENT,neig,num_shap,i,
				  dofs_loc1);
		uc2=dofs_loc1[0]+dofs_loc1[1]/3.0+dofs_loc1[2]/3.0;
	      }
	      else{
		uc2=pdr_average_sol_el(Problem_id,neig);
	      }
	      
	    } /* end if not on the boundary */
	    else{
	      if (bc_type2==PDC_BC_DIRI) {
	        pdr_bc_diri_coeff(el_faces[2], el_mate, xc2, NULL,
				     time,NULL,NULL,NULL,NULL,&uc2);
	      }	      
	      else if (bc_type2==PDC_BC_MIXED) {
	        pdr_bc_mixed_coeff(el_faces[2], el_mate, xc2, NULL,
				     time,NULL,NULL,NULL,NULL,&uc2,&daux);
	      }	      
	    }

/* check whether elements form structured layers */
/* we assume neighbors are linear prisms as well */
#ifdef DEBUG
	    if((el_neig[1]>0&&
		(fabs(xc1[0]-xc0[0])>SMALL||fabs(xc1[1]-xc0[1])>SMALL))||
	       (el_neig[2]>0&&
		(fabs(xc2[0]-xc0[0])>SMALL||fabs(xc2[1]-xc0[1])>SMALL))){

	      printf("Elements not in layers!\n");
	      printf("el0 - center: %lf, %lf, %lf\n", xc0[0], xc0[1],xc0[2]);
	      printf("el1 - center: %lf, %lf, %lf\n", xc1[0], xc1[1],xc1[2]);
	      printf("el2 - center: %lf, %lf, %lf\n", xc2[0], xc2[1],xc2[2]);

              return(0);
	    }
#endif


/* limit only if gradient in z greater than treshold value */
	    if(fabs(dofs_loc[3])>MAGIC1*zsize*zsize){

/* limited value - gradient in z direction */
	      utilda_z = dofs_loc[3];

	      if(el_neig[1]>0) ur1=uc0-uc1;
	      else if(bc_type1==PDC_BC_DIRI||bc_type1==PDC_BC_MIXED) 
		ur1=2.0*(uc0-uc1);
	      else ur1 = utilda_z;
	      if(el_neig[2]>0) ur2=uc2-uc0;
	      else if(bc_type2==PDC_BC_DIRI||bc_type2==PDC_BC_MIXED) 
		ur2=2.0*(uc2-uc0);
	      else ur2 = utilda_z;

/*kbw
printf("Before limiting in z (zsize %lf):\n",zsize);
printf("u,z = %lf, ur1 = %lf, ur2 = %lf\n",
utilda_z,ur1,ur2);
kew*/

/* compute minmod limited value */
	      if((utilda_z>0&&ur1>0&&ur2>0)||(utilda_z<0&&ur1<0&&ur2<0)){
		daux=utm_min(fabs(ur2),utm_min(fabs(ur1),fabs(utilda_z)));
		utilda_z = daux*utilda_z/fabs(utilda_z);
	      }
	      else utilda_z=0;

/*kbw
if(fabs(dofs_loc[3]-utilda_z)>SMALL){
  printf("ELEMENT %d LIMITED in Z\n",iel);
  printf("After limiting in z: old value %lf, new value %lf\n",
	 dofs_loc[3],utilda_z);
}
/*kew*/

/* write back corrected degrees of freedom */
	      if(fabs(dofs_loc[3]-utilda_z)>SMALL) {


		if(num_shap>4){
/* zero all higher order degrees of freedom */
		  for(i=4;i<num_shap;i++) {
		    dofs_loc[i]=0.0;
		  }

/* modify the first dof to get the proper average value */
		  dofs_loc[0]=uc0-dofs_loc[1]/3.0-dofs_loc[2]/3.0;
		}

		dofs_loc[3] = utilda_z;
		i=1;
		apr_write_ent_dofs(field_id, APC_ELEMENT, iel, num_shap, i, 
				   dofs_loc);

	      } /* end if limiting took place */

	    } /* end if gradient in z big enough */

	  } /* end if z gradient exist */

	  if(fabs(dofs_loc[1])+fabs(dofs_loc[2])>SMALL){

/* perform limiting in xy plane */

/* take values at centers of three neighbors */

/* get neighbor across first side */
	    neig=el_neig[3];

/* get bc type of first (bottom) face */
	    bc_type1=pdr_get_bc_type(mmr_fa_bc(mesh_id,el_faces[3]));

/* if not on the boundary */
	    if(neig>0){

/* get the coordinates of element nodes in the right order */
	      mmr_el_node_coor(mesh_id, neig, NULL, node_coor1);

/* center coordinates */
	      xc1[0]=(node_coor1[0]+node_coor1[3]+node_coor1[6]+
		      node_coor1[9]+node_coor1[12]+node_coor1[15])/6.0;
	      xc1[1]=(node_coor1[1]+node_coor1[4]+node_coor1[7]+
		      node_coor1[10]+node_coor1[13]+node_coor1[16])/6.0;
	      xc1[2]=(node_coor1[2]+node_coor1[5]+node_coor1[8]+
		      node_coor1[11]+node_coor1[14]+node_coor1[17])/6.0;

/* average value (for p=1 value at geometrical center) */
	      iaux=apr_get_ent_pdeg(field_id, APC_ELEMENT, neig);
	      if(iaux==1){
/* get solution degrees of freedom */
		i=1;
		apr_read_ent_dofs(field_id,APC_ELEMENT,neig,num_shap,i,
				  dofs_loc1);
		uc1=dofs_loc1[0]+dofs_loc1[1]/3.0+dofs_loc1[2]/3.0;
	      }
	      else{
		uc1=pdr_average_sol_el(Problem_id,neig);
	      }

	    } /* end if not on the boundary */
	    else{
	      if (bc_type1==PDC_BC_DIRI||bc_type1==PDC_BC_MIXED) {
		
/* center of face coordinates */
		xc1[0]=(node_coor[0]+node_coor[3]+
			node_coor[9]+node_coor[12])/4.0;
		xc1[1]=(node_coor[1]+node_coor[4]+
			node_coor[10]+node_coor[13])/4.0;
		xc1[2]=(node_coor[2]+node_coor[5]+
			node_coor[11]+node_coor[14])/4.0;

		if (bc_type1==PDC_BC_DIRI) {
		  pdr_bc_diri_coeff(el_faces[3], el_mate, xc1, NULL,
				     time,NULL,NULL,NULL,NULL,&uc1);
		}
		else if (bc_type1==PDC_BC_MIXED) {
		  pdr_bc_mixed_coeff(el_faces[3], el_mate, xc1, NULL,
				     time,NULL,NULL,NULL,NULL,&uc1,&daux);
		}

	      }
	    } /* end if on the boundary */

/*kbw
printf("Element1 %d - center %lf %lf %lf, value %lf\n",
neig,xc1[0],xc1[1],xc1[2],uc1);
/*kew*/

/* get neighbor across second side */
	    neig=el_neig[4];

/* get bc type of first (bottom) face */
	    bc_type2=pdr_get_bc_type(mmr_fa_bc(mesh_id,el_faces[4]));

/* if not on the boundary */
	    if(neig>0){

/* get the coordinates of element nodes in the right order */
	      mmr_el_node_coor(mesh_id, neig, NULL, node_coor1);

/* center coordinates */
	      xc2[0]=(node_coor1[0]+node_coor1[3]+node_coor1[6]+
		      node_coor1[9]+node_coor1[12]+node_coor1[15])/6.0;
	      xc2[1]=(node_coor1[1]+node_coor1[4]+node_coor1[7]+
		      node_coor1[10]+node_coor1[13]+node_coor1[16])/6.0;
	      xc2[2]=(node_coor1[2]+node_coor1[5]+node_coor1[8]+
		      node_coor1[11]+node_coor1[14]+node_coor1[17])/6.0;

/* average value (for p=1 value at geometrical center) */
	      iaux=apr_get_ent_pdeg(field_id, APC_ELEMENT, neig);
	      if(iaux==1){
/* get solution degrees of freedom */
		i=1;
		apr_read_ent_dofs(field_id,APC_ELEMENT,neig,num_shap,i,
				  dofs_loc1);
		uc2=dofs_loc1[0]+dofs_loc1[1]/3.0+dofs_loc1[2]/3.0;
	      }
	      else{
		uc2=pdr_average_sol_el(Problem_id,neig);
	      }
	    } /* end if not on the boundary */
	    else{
	      if (bc_type2==PDC_BC_DIRI||bc_type2==PDC_BC_MIXED) {

		xc2[0]=(node_coor[3]+node_coor[6]+
			node_coor[12]+node_coor[15])/4.0;
		xc2[1]=(node_coor[4]+node_coor[7]+
			node_coor[13]+node_coor[16])/4.0;
		xc2[2]=(node_coor[5]+node_coor[8]+
			node_coor[14]+node_coor[17])/4.0;

		if (bc_type2==PDC_BC_DIRI) {
		  pdr_bc_diri_coeff(el_faces[4], el_mate, xc2, NULL,
				     time,NULL,NULL,NULL,NULL,&uc2);
		}
		if (bc_type2==PDC_BC_MIXED) {
		  pdr_bc_mixed_coeff(el_faces[4], el_mate, xc2, NULL,
				     time,NULL,NULL,NULL,NULL,&uc2,&daux);
		}

	      }
	    } /* end if on the boundary */

/*kbw
printf("Element2 %d - center %lf %lf %lf, value %lf\n",
neig,xc2[0],xc2[1],xc2[2],uc2);
/*kew*/


/* get neighbor across third side */
	    neig=el_neig[5];

/* get bc type first (bottom) face */
	    bc_type3=pdr_get_bc_type(mmr_fa_bc(mesh_id,el_faces[5]));

/* if not on the boundary */
	    if(neig>0){

/* get the coordinates of element nodes in the right order */
	      mmr_el_node_coor(mesh_id, neig, NULL, node_coor1);

/* center coordinates */
	      xc3[0]=(node_coor1[0]+node_coor1[3]+node_coor1[6]+
		      node_coor1[9]+node_coor1[12]+node_coor1[15])/6.0;
	      xc3[1]=(node_coor1[1]+node_coor1[4]+node_coor1[7]+
		      node_coor1[10]+node_coor1[13]+node_coor1[16])/6.0;
	      xc3[2]=(node_coor1[2]+node_coor1[5]+node_coor1[8]+
		      node_coor1[11]+node_coor1[14]+node_coor1[17])/6.0;

/* average value (for p=1 value at geometrical center) */
	      iaux=apr_get_ent_pdeg(field_id, APC_ELEMENT, neig);
	      if(iaux==1){
/* get solution degrees of freedom */
		i=1;
		apr_read_ent_dofs(field_id,APC_ELEMENT,neig,num_shap,i,
				  dofs_loc1);
		uc3=dofs_loc1[0]+dofs_loc1[1]/3.0+dofs_loc1[2]/3.0;
	      }
	      else{
		uc3=pdr_average_sol_el(Problem_id,neig);
	      }
	    } /* end if not on the boundary */
	    else{
	      if (bc_type3==PDC_BC_DIRI||bc_type3==PDC_BC_MIXED) {

		xc3[0]=(node_coor[0]+node_coor[6]+
			node_coor[9]+node_coor[15])/4.0;
		xc3[1]=(node_coor[1]+node_coor[7]+
			node_coor[10]+node_coor[16])/4.0;
		xc3[2]=(node_coor[2]+node_coor[8]+
			node_coor[11]+node_coor[17])/4.0;

		if (bc_type3==PDC_BC_DIRI) {
		  pdr_bc_diri_coeff(el_faces[5], el_mate, xc3, NULL,
				     time,NULL,NULL,NULL,NULL,&uc3);
		}
		else if (bc_type3==PDC_BC_MIXED) {
		  pdr_bc_mixed_coeff(el_faces[5], el_mate, xc3, NULL,
				     time,NULL,NULL,NULL,NULL,&uc3,&daux);
		}

	      }
	    } /* end if on the boundary */

/*kbw
printf("Element3 %d - center %lf %lf %lf, value %lf\n",
       neig,xc3[0],xc3[1],xc3[2],uc3);
/*kew*/

#ifdef DEBUG
/* check whether elements form structured layers */
	    if((el_neig[3]>0&&fabs(xc1[2]-xc0[2])>SMALL)||
	       (el_neig[4]>0&&fabs(xc2[2]-xc0[2])>SMALL)||
	       (el_neig[5]>0&&fabs(xc3[2]-xc0[2])>SMALL)){
	      
	      printf("Elements not in layers!\n");
	      printf("el0 %d - center: %lf, %lf, %lf\n",
		     iel,xc0[0],xc0[1],xc0[2]);
	      printf("el1 %d - center: %lf, %lf, %lf\n",
		     el_neig[3],xc1[0],xc1[1],xc1[2]);
	      printf("el2 %d - center: %lf, %lf, %lf\n",
		     el_neig[4],xc2[0],xc2[1],xc2[2]);
	      printf("el3 %d - center: %lf, %lf, %lf\n",
		     el_neig[5],xc3[0],xc3[1],xc3[2]);
	      
	      return(0);
	    }
#endif

/* limited value: gradient in x direction */
	    daux = (node_coor[3]-node_coor[0])*(node_coor[7]-node_coor[1]) -
	           (node_coor[4]-node_coor[1])*(node_coor[6]-node_coor[0]);
	    utilda_x = -( (node_coor[4]-node_coor[1])*dofs_loc[2] -
                          (node_coor[7]-node_coor[1])*dofs_loc[1] ) / daux ;

/* limit only if gradient in x greater than treshold value */
	    if(fabs(utilda_x)>MAGIC1*xysize*xysize){


/* project three gradients on the x direction to get reference values */
	      if((el_neig[3]>0||bc_type1==PDC_BC_DIRI||bc_type1==PDC_BC_MIXED)
	      &&(el_neig[4]>0||bc_type2==PDC_BC_DIRI||bc_type2==PDC_BC_MIXED))
		{
		daux = (xc1[0]-xc0[0])*(xc2[1]-xc0[1]) -
		       (xc1[1]-xc0[1])*(xc2[0]-xc0[0]);
		if(fabs(daux)>0.01*xysize){
		  ur1 = -( (xc1[1]-xc0[1])*(uc2-uc0) -
			   (xc2[1]-xc0[1])*(uc1-uc0) ) * coeff3 / daux ;
		}
		else{
		  ur1 = utilda_x;
		}
	      }
	      else ur1 = utilda_x;
	      if((el_neig[4]>0||bc_type2==PDC_BC_DIRI||bc_type2==PDC_BC_MIXED)
	      &&(el_neig[5]>0||bc_type3==PDC_BC_DIRI||bc_type3==PDC_BC_MIXED))
		{
		daux = (xc2[0]-xc0[0])*(xc3[1]-xc0[1]) -
		       (xc2[1]-xc0[1])*(xc3[0]-xc0[0]);
		if(fabs(daux)>0.01*xysize){
		  ur2 = -( (xc2[1]-xc0[1])*(uc3-uc0) -
			   (xc3[1]-xc0[1])*(uc2-uc0) ) * coeff3 / daux ;
		}
		else{
		  ur2 = utilda_x;
		}
	      }
	      else ur2 = utilda_x;
	      if((el_neig[5]>0||bc_type3==PDC_BC_DIRI||bc_type3==PDC_BC_MIXED)
	      &&(el_neig[3]>0||bc_type1==PDC_BC_DIRI||bc_type1==PDC_BC_MIXED))
		{
		daux = (xc3[0]-xc0[0])*(xc1[1]-xc0[1]) -
		       (xc3[1]-xc0[1])*(xc1[0]-xc0[0]);
		if(fabs(daux)>0.01*xysize){
		  ur3 = -( (xc3[1]-xc0[1])*(uc1-uc0) -
			   (xc1[1]-xc0[1])*(uc3-uc0) ) * coeff3 / daux ;
		}
		else{
		  ur3 = utilda_x;
		}
	      }
	      else ur3 = utilda_x;

/*kbw
printf("Before limiting in x, el %d xysize %lf:\n",iel,xysize);
printf("u,x = %lf, ur1 = %lf, ur2 = %lf, ur3 = %lf\n",
       utilda_x,ur1,ur2,ur3);
/*kew*/

/* compute minmod limited value */
	      if((utilda_x>0&&ur1>0&&ur2>0&&ur3>0)||
		 (utilda_x<0&&ur1<0&&ur2<0&&ur3<0)){
		daux=utm_min(fabs(ur1),
			    utm_min(fabs(ur2),
				   utm_min(fabs(ur3),fabs(utilda_x))));
		utilda_x = daux*utilda_x/fabs(utilda_x);
	      }
	      else utilda_x=0;


/*kbw
printf("After limiting in x: new value u,x = %lf\n", utilda_x);
/*kew*/

	    }

/* limited value: gradient in y direction */
	    daux = (node_coor[3]-node_coor[0])*(node_coor[7]-node_coor[1]) -
	           (node_coor[4]-node_coor[1])*(node_coor[6]-node_coor[0]);
	    utilda_y = -( (node_coor[6]-node_coor[0])*dofs_loc[1] -
			  (node_coor[3]-node_coor[0])*dofs_loc[2] ) / daux ;

/* limit only if gradient in y greater than treshold value */
	    if(fabs(utilda_y)>MAGIC1*xysize*xysize){

/* project three gradients on the y direction to get reference values */
	      if((el_neig[3]>0||bc_type1==PDC_BC_DIRI||bc_type1==PDC_BC_MIXED)
	      &&(el_neig[4]>0||bc_type2==PDC_BC_DIRI||bc_type2==PDC_BC_MIXED))
		{
		daux = (xc1[0]-xc0[0])*(xc2[1]-xc0[1]) -
		       (xc1[1]-xc0[1])*(xc2[0]-xc0[0]);
		if(fabs(daux)>0.01*xysize){
		  ur1 = -( (xc2[0]-xc0[0])*(uc1-uc0) -
			   (xc1[0]-xc0[0])*(uc2-uc0) ) / daux ;
		}
		else{
		  ur1 = utilda_y;
		}
	      }
	      else ur1 = utilda_y;
	      if((el_neig[4]>0||bc_type2==PDC_BC_DIRI||bc_type2==PDC_BC_MIXED)
	      &&(el_neig[5]>0||bc_type3==PDC_BC_DIRI||bc_type3==PDC_BC_MIXED))
		{
		daux = (xc2[0]-xc0[0])*(xc3[1]-xc0[1]) -
		       (xc2[1]-xc0[1])*(xc3[0]-xc0[0]);
		if(fabs(daux)>0.01*xysize){
		  ur2 = -( (xc3[0]-xc0[0])*(uc2-uc0) -
			   (xc2[0]-xc0[0])*(uc3-uc0) ) / daux ;
		}
		else{
		  ur2 = utilda_y;
		}
	      }
	      else ur2 = utilda_y;
	      if((el_neig[5]>0||bc_type3==PDC_BC_DIRI||bc_type3==PDC_BC_MIXED)
	      &&(el_neig[3]>0||bc_type1==PDC_BC_DIRI||bc_type1==PDC_BC_MIXED))
		{
		daux = (xc3[0]-xc0[0])*(xc1[1]-xc0[1]) -
                       (xc3[1]-xc0[1])*(xc1[0]-xc0[0]);
		if(fabs(daux)>0.01*xysize){
		  ur3 = -( (xc1[0]-xc0[0])*(uc3-uc0) -
			   (xc3[0]-xc0[0])*(uc1-uc0) ) / daux ;
		}
		else{
		  ur3 = utilda_y;
		}
	      }
	      else ur3 = utilda_y;

/*kbw
printf("Before limiting in y (xysize %lf):\n",xysize);
printf("u,y = %lf, ur1 = %lf, ur2 = %lf, ur3 = %lf\n",
utilda_y,ur1,ur2,ur3);
/*kew*/

/* compute minmod limited value */
	      if((utilda_y>0&&ur1>0&&ur2>0&&ur3>0)||
		 (utilda_y<0&&ur1<0&&ur2<0&&ur3<0)){
		daux=utm_min(fabs(ur1),
			    utm_min(fabs(ur2),
				   utm_min(fabs(ur3),fabs(utilda_y))));
		utilda_y = daux*utilda_y/fabs(utilda_y);
	      }
	      else utilda_y=0;

/*kbw
printf("After limiting in y: new value u,y = %lf\n", utilda_y);
/*kew*/

	    }

/* compute back element dofs */
	    dofs_loc[1] = utilda_x*(node_coor[3]-node_coor[0]) +
	                  utilda_y*(node_coor[4]-node_coor[1]);
	    dofs_loc[2] = utilda_x*(node_coor[6]-node_coor[0]) +
                          utilda_y*(node_coor[7]-node_coor[1]);

/* correct first dof so that the average remains the same */
	    daux=dofs_loc[0];
	    dofs_loc[0]=uc0-dofs_loc[1]/3.0-dofs_loc[2]/3.0;

/* if there where any changes reset remaining dofs */
	    if(fabs(daux-dofs_loc[0])>SMALL){
	      if(num_shap>4){
/* zero all higher order degrees of freedom */
		for(i=4;i<num_shap;i++) {
		  dofs_loc[i]=0.0;
		}
	      }
	    }
/*kbw
if(fabs(daux-dofs_loc[0])>SMALL){
printf("After SLOPE_LIMIT - dofs:\n");
for(i=0;i<4;i++) printf("%20.15lf",dofs_loc[i]);
printf("\n");
getchar();
}
/*kew*/

/* write back corrected degrees of freedom */
	    i=1;
	    apr_write_ent_dofs(field_id, APC_ELEMENT, iel, num_shap, i, 
			       dofs_loc);
	  
	  } /* end if gradients in x and y big enough */

	} /* end if linear prism */


    }/* end for each element */


  }

  return(1);

}

