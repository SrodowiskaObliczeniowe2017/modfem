/*
 * defs.h
 *
 *  Created on: 2011-04-07
 *      Author: Paweł Macioł
 */
#ifndef _FV_DEFS_h__
#define _FV_DEFS_h__


#define FV_STRING(str)		#str
#define FV_STR(str)			FV_STRING(str)

#define FV_SIZEOF_ARR(arr) 	(sizeof(arr)/sizeof((arr)[0]))

/* This stuff is for special type element idx = el_id + attribs */
#define COUNT_BITS(w)		(sizeof(w) << 3)
#define BOUNDARY_BIT_POS(w)	(COUNT_BITS(w) -1)
#define ACTIVE_BIT_POS(w)	(COUNT_BITS(w) -2)
#define STFACE_BIT_POS(w)	(COUNT_BITS(w) -7)
#define ELTYPE_BIT_POS(w)	(COUNT_BITS(w) -8)
#define EL_ID_BIT_COUNT(w)	(COUNT_BITS(w) -9)

/* This is for element faces*/
#define FACE_TYPE_POS(w)	(COUNT_BITS(w) -1)

/* check status of bit */
#define IS_BIT_SET(w, b)     ((w) & (1U << (b)))

/* set bit (b) in word (w) */
#define SET_BIT(w, b)        ((w) |= (1U << (b)))

/* unset bit (b) in word (w) */
#define RESET_BIT(w, b)      ((w) &= ~(1U << (b)))

/* set/unset bit from position (b) with offset (o) in word */
#define SET_OFFSET_BIT(w, b, o) SET_BIT(w,(b + o))
#define USET_OFFSET_BIT(w, b, o) RSET_BIT(w,(b + o))

/* sign elemnet as active/inactive */
#define SET_ACTIVE(w)		SET_BIT(w, ACTIVE_BIT_POS(w))
#define UNSET_ACTIVE(w)		RESET_BIT(w, ACTIVE_BIT_POS(w))
#define IS_ACTIVE(w)		IS_BIT_SET(w, ACTIVE_BIT_POS(w))

/* sign elemnet as boundary */
#define SET_BOUND(w)		SET_BIT(w, BOUNDARY_BIT_POS(w))
#define UNSET_BOUND(w)		RESET_BIT(w, BOUNDARY_BIT_POS(w))
#define IS_BOUND(w)			IS_BIT_SET(w, BOUNDARY_BIT_POS(w))

/* activate state of face nr */
#define SET_STFACE0(w)		SET_OFFSET_BIT(w, STFACE_BIT_POS(w), 0)
#define SET_STFACE1(w)		SET_OFFSET_BIT(w, STFACE_BIT_POS(w), 1)
#define SET_STFACE2(w)		SET_OFFSET_BIT(w, STFACE_BIT_POS(w), 2)
#define SET_STFACE3(w)		SET_OFFSET_BIT(w, STFACE_BIT_POS(w), 3)
#define SET_STFACE4(w)		SET_OFFSET_BIT(w, STFACE_BIT_POS(w), 4)
#define SET_STFACEN(w, n)	SET_OFFSET_BIT(w, STFACE_BIT_POS(w), n)

/* inactivate state of face nr */
#define UNSET_STFACE0(w)	USET_OFFSET_BIT(w, STFACE_BIT_POS(w), 0)
#define UNSET_STFACE1(w)	USET_OFFSET_BIT(w, STFACE_BIT_POS(w), 1)
#define UNSET_STFACE2(w)	USET_OFFSET_BIT(w, STFACE_BIT_POS(w), 2)
#define UNSET_STFACE3(w)	USET_OFFSET_BIT(w, STFACE_BIT_POS(w), 3)
#define UNSET_STFACE4(w)	USET_OFFSET_BIT(w, STFACE_BIT_POS(w), 4)

#define IS_SET_FACEN(w, n)		IS_BIT_SET(w, STFACE_BIT_POS(w) + (n))

/* set element type */
#define SET_TETRA(w)		RESET_BIT(w, ELTYPE_BIT_POS(w))
#define SET_PRISM(w)		SET_BIT(w, ELTYPE_BIT_POS(w))
#define IS_PRISM(w)			IS_BIT_SET(w, ELTYPE_BIT_POS(w))

/* 1 byte of attribute info */
#define ATTRIB_MASK			0xFF

/* Extract element id */
#define ELEM_ID(w)			(w & ~(ATTRIB_MASK << EL_ID_BIT_COUNT(w)))

/* Extract attribute info */
#define ATTRIBS(w)			(w >> EL_ID_BIT_COUNT(w))

#define MIN(a,b) ((a) < (b)) ? (a) : (b)
#define MAX(a,b) ((a) > (b)) ? (a) : (b)

#define INVALID_LOCATION 0xFFFFFFFF
#define FV_SMALL 1e-10
#define FV_LARGE 10e10
#define FV_LARGE_FLOAT FV_LARGE##f
/* Define namespace*/

//#ifndef fv_min
//# define fv_min(x,y)	(x) > (y) ? (y) : (x)
//# define fv_max(x,y)	fv_min(y,x)
//#endif
template<typename T>
inline const T& fv_min(const T& lh,const T& rh) {
	return (lh < rh ? lh : rh );
}

template<typename T>
inline const T& fv_max(const T& lh,const T& rh) {
	return (lh > rh ? lh : rh );
}

#ifndef fv_abs
# define fv_abs(x)		((x) < 0) ? -(x) : (x)
#endif


#ifndef NULL
#define NULL (void *)0
#endif

template<typename T>
inline void FV_FREE_ARR(T* ptr) {
	if(ptr) { delete [] ptr; ptr = 0; }
}

template<typename T>
inline void FV_FREE_PTR(T*& ptr){
	if(ptr) { delete  ptr; ptr = 0; }
}

template<typename T>
inline void FV_FREE_MALOC(T*& ptr) {
	free(ptr); ptr = NULL;
}

typedef struct {
	int curr_mesh_type;
	int curr_approx_type;
} nodule_params;

#define LARGE_F 		10e10
#define MAX_ISO_VALUES	32




#endif /* _FV_DEFS_h__ */

