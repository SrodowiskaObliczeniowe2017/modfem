#ifndef _DEFS_H_
#define _DEFS_H_

/* Types of mesh entities */
//const int MMC_QUAD          = 4;   /* quadrilateral element or face */
//const int MMC_TRIA          = 3;   /* triangular element or face */
//const int MMC_TETRA         = 7;   /* tetrahedral element */
//const int MMC_PRISM         = 5;   /* prismatic element */
//const int MMC_BRICK         = 6;   /* hexahedral element */
//
///* Status indicators */
//const int MMC_ACTIVE        = 1;   /* active mesh entity */
//const int MMC_INACTIVE      = -1;   /* inactive (refined) mesh entity */
//const int MMC_FREE          = 0;   /* free space in data structure */

//const int APC_NO_DOFS       = -1;

/* Refinement types */
//const int MMC_NOT_REF       = 0;   /* not refined */
//const int MMC_REF_ISO       = 1;   /* isotropic refinement */

/* Basis functions types */
#define APC_TENSOR        0
#define APC_COMPLETE      1

/* Initialization options */
//const int APC_ZERO          = 0;
//const int APC_READ          = 1;
//const int APC_INIT          = 2;
//
///* standard macro for max and min and abs */
#define ut_max(x,y) ((x)>(y)?(x):(y))
#define ut_min(x,y) ((x)<(y)?(x):(y))
#define ut_abs(x)   ((x)<0?-(x):(x))
//
///* Constants */
////const int MMC_MAXELFAC       = 5;  /* maximal number of faces of an element */
////const int MMC_MAXFAVNO        = 4;  /* maximal number of vertices of a face */
////const int MMC_MAXELVNO        = 6;  /* maximal number of vertices of an element */
////const int MMC_MAXELSONS       = 8;  /* maximal number of faces of an element */
//
//const int APC_MAXEQ          = 1;
//
////#define APC_MAXELP_COMP 9
//const int APC_MAXELP_COMP          = 9;
////#define APC_MAXELP_TENS 909
//const int APC_MAXELP_TENS          = 909;
//
#define APC_MAXELVD		600
//const int APC_MAXELVD          = 600;
////#define APC_MAXELSD APC_MAXEQ*APC_MAXELVD
////const int APC_MAXELP_TENS          = 909;
//
///* Other */
//const int MMC_NO_FATH       = 0;   /* no father indicator */
//const int MMC_SAME_ORIENT   = 1;   /* indicator for the same orientation */
//const int MMC_OPP_ORIENT   = -1;   /* indicator for the opposite orientation */
////extern const int MMC_NO_FATH;		/* no father indicator */
////extern const int MMC_SAME_ORIEN;		 /* indicator for the same orientation */
////extern const int MMC_OPP_ORIENT;   /* indicator for the opposite orientation */
//
//
////#define APC_MAXELSD APC_MAXEQ*APC_MAXELVD
////const int APC_MAXELSD   = APC_MAXEQ*APC_MAXELVD;
//
//#define APC_FALSE	0
//#define APC_TRUE	1


#endif /* _DEFS_H_
*/