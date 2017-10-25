#ifndef _lin_alg_intf_
#define _lin_alg_intf_

#ifdef __cplusplus
extern "C" {
#endif 

/* BLAS and LAPACK names with or without underscore */
#ifdef WITHOUT_
#define daxpy_ daxpy
#define dcopy_ dcopy
#define ddot_ ddot
#define dnrm2_ dnrm2
#define drot_ drot
#define drotg_ drotg
#define dscal_ dscal
#define dgemv_ dgemv
#define dtrsv_ dtrsv
#define dgetrf_ dgetrf
#define dgetrs_ dgetrs
#endif

/* BLAS 1 */
void   daxpy_(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
void   dcopy_(int *n,double *x,int *incx,double *y,int *incy);
double ddot_(int *n,double *x,int *incx,double *y,int *incy);
double dnrm2_(int *n,double *x,int *incx);
void   drot_(int *n,double *x,int *incx,double *y,int *incy,double *c,double *s);
void   drotg_(double *a,double *b,double *c,double *s);
void   dscal_(int *n,double *a,double *x,int *incx);
void   DAXPY_(int *n,double *alpha,double *x,int *incx,double *y,int *incy);
void   DCOPY_(int *n,double *x,int *incx,double *y,int *incy);
double DDOT_(int *n,double *x,int *incx,double *y,int *incy);
double DNRM2_(int *n,double *x,int *incx);
void   DROT_(int *n,double *x,int *incx,double *y,int *incy,double *c,double *s);
void   DROTG_(double *a,double *b,double *c,double *s);
void   DSCAL_(int *n,double *a,double *x,int *incx);

/* BLAS 2 */
void dgemv_(char *trans,int *m,int *n,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
void dtrsv_(char *uplo,char *trans,char *diag,int *n,double *a,int *lda,double *x,int *incx);
void DGEMV_(char *trans,int *m,int *n,double *alpha,double *a,int *lda,double *x,int *incx,double *beta,double *y,int *incy);
void DTRSV_(char *uplo,char *trans,char *diag,int *n,double *a,int *lda,double *x,int *incx);

/* LAPACK */
void    DGETRF_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
void    dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info);
void    DGETRS_(char *trans,int *n,int *nrhs,double *a,int *lda,int *ipiv,double *b,int *ldb,int *info);
void    dgetrs_(char *trans,int *n,int *nrhs,double *a,int *lda,int *ipiv,double *b,int *ldb,int *info);

#ifdef __cplusplus
}
#endif 

#endif
