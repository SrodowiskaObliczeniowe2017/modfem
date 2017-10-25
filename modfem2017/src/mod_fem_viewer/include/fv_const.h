#ifndef _FV_CONST_H_
#define _FV_CONST_H_ 

#ifndef APC_FALSE
#define APC_FALSE 0
#endif

#ifndef APC_TRUE
#define APC_TRUE 1
#endif

#ifndef FV_USE_EXTERNALS

/* Types of mesh entities */
const int FV_MMC_QUAD          = 4;   /* quadrilateral element or face */
const int FV_MMC_TRIA          = 3;   /* triangular element or face */
const int FV_MMC_TETRA         = 7;   /* tetrahedral element */
const int FV_MMC_PRISM         = 5;   /* prismatic element */
const int FV_MMC_BRICK         = 6;   /* hexahedral element */

///* Status indicators */
const int FV_MMC_ACTIVE        = 1;   /* active mesh entity */
const int FV_MMC_INACTIVE      = -1;   /* inactive (refined) mesh entity */
const int FV_MMC_FREE          = 0;   /* free space in data structure */

//const int APC_NO_DOFS       = -1;

/* Refinement types */
const int FV_MMC_NOT_REF       = 0;   /* not refined */
const int FV_MMC_REF_ISO       = 1;   /* isotropic refinement */

/* Basis functions types */
const int APC_TENSOR        = 0;
const int APC_COMPLETE      = 1;

/* Initialization options */
//const int APC_ZERO          = 0;
//const int APC_READ          = 1;
//const int APC_INIT          = 2;

/* standard macro for max and min and abs */
#define ut_max(x,y) ((x)>(y)?(x):(y))
#define ut_min(x,y) ((x)<(y)?(x):(y))
#define ut_abs(x)   ((x)<0?-(x):(x))

/* Constants */
//const int MMC_MAXELFAC       = 5;  /* maximal number of faces of an element */
//const int MMC_MAXFAVNO        = 4;  /* maximal number of vertices of a face */
//const int MMC_MAXELVNO        = 6;  /* maximal number of vertices of an element */
//const int MMC_MAXELSONS       = 8;  /* maximal number of faces of an element */

const int APC_MAXEQ          = 1;

//#define APC_MAXELP_COMP 9
//const int APC_MAXELP_COMP          = 9;
//#define APC_MAXELP_TENS 909
//const int APC_MAXELP_TENS          = 909;

//#define APC_MAXELVD 600
//const int APC_MAXELVD          = 600;
//#define APC_MAXELSD APC_MAXEQ*APC_MAXELVD
//const int APC_MAXELP_TENS          = 909;

/* Other */
const int FV_MMC_NO_FATH       = 0;   /* no father indicator */
const int FV_MMC_SAME_ORIENT   = 1;   /* indicator for the same orientation */
const int FV_MMC_OPP_ORIENT   = -1;   /* indicator for the opposite orientation */
//extern const int MMC_NO_FATH;		/* no father indicator */
//extern const int MMC_SAME_ORIEN;		 /* indicator for the same orientation */
//extern const int MMC_OPP_ORIENT;   /* indicator for the opposite orientation */

const int FV_MMC_BOUNDARY = 0; /* boundary indicator */
//#define APC_MAXELSD APC_MAXEQ*APC_MAXELVD
//const int APC_MAXELSD   = APC_MAXEQ*APC_MAXELVD;



#endif 




#endif /* _FV_CONST_H_
*/
