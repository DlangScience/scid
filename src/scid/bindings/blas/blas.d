/** BLAS bindings for D.

    Authors:    William V. Baxter III (with slight modifications by
                Lars Tandle Kyllingstad).
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.bindings.blas.blas;


public import scid.bindings.blas.types;




version (FORTRAN_FLOAT_FUNCTIONS_RETURN_DOUBLE) {
    alias f_double float_ret_t;
} else {
    alias f_float float_ret_t;
}


// Prototypes for the raw Fortran interface to BLAS
extern(C):

/* BLAS routines */

/** Level 1 BLAS */

/** Generate plane (Givens) rotation
    Given a and b, compute the elements of a rotation matrix such that
          _      _     _   _    _   _
          | c  s |     | a |    | r |
          |-s  c | *   | b | =  | 0 |
          -      -     -   -    -   -
     where
     r = +/- sqrt (a^2  + b^2 ) and c^2 + s^2  = 1   (real case)
     or
     r = (a/sqrt(conj(a)*a  + conj(b)*b)) * sqrt(conj(a)*a + conj(b)*b)
*/
void srotg_(f_float *a, f_float *b, f_float *c, f_float *s);
void drotg_(f_double *a, f_double *b, f_double *c, f_double *s);
void crotg_(f_cfloat *a, f_cfloat *b, f_float *c, f_cfloat *s);
void zrotg_(f_cdouble *a, f_cdouble *b, f_double *c, f_cdouble *s);

/// Generate modified plane (Givens) rotation
void drotmg_(f_double *d1, f_double *d2, f_double *b1, f_double *b2, f_double *param);
void srotmg_(f_float *d1, f_float *d2, f_float *b1, f_float *b2, f_float *param);

/// Apply plane (Givens) rotation
///             _      _
///     x_i  := | c  s | * x_i
///     y_i     |-s  c |   y_i
///             -      -
void srot_(f_int *n, f_float *x, f_int *incx, f_float *y, f_int *incy, f_float *c, f_float *s);
void drot_(f_int *n, f_double *x, f_int *incx, f_double *y, f_int *incy, f_double *c, f_double *s);
void csrot_(f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy, f_float *c, f_float *s);
void zdrot_(f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy, f_double *c, f_double *s);

/// Apply modified plane (Givens) rotation
void srotm_(f_int *n, f_float *x, f_int *incx, f_float *y, f_int *incy, f_float *param);
void drotm_(f_int *n, f_double *x, f_int *incx, f_double *y, f_int *incy, f_double *param);

/// Swap the values contained in x and y 
///     x <-> y
void sswap_(f_int *n, f_float *x, f_int *incx, f_float *y, f_int *incy);
void dswap_(f_int *n, f_double *x, f_int *incx, f_double *y, f_int *incy);
void cswap_(f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
void zswap_(f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);

/// x := alpha * x
void sscal_(f_int *n, f_float *alpha, f_float *x, f_int *incx);
void dscal_(f_int *n, f_double *alpha, f_double *x, f_int *incx);
void cscal_(f_int *n, f_cfloat *alpha, f_cfloat *x, f_int *incx);
void csscal_(f_int *n, f_float *alpha, f_cfloat *x, f_int *incx);
void zscal_(f_int *n, f_cdouble *alpha, f_cdouble *x, f_int *incx);
void zdscal_(f_int *n, f_double *alpha, f_cdouble *x, f_int *incx);

/// y := x
void scopy_(f_int *n, f_float *x, f_int *incx, f_float *y, f_int *incy);
void dcopy_(f_int *n, f_double *x, f_int *incx, f_double *y, f_int *incy);
void ccopy_(f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
void zcopy_(f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);

/// y := alpha * x + y
void saxpy_(f_int *n, f_float *alpha, f_float *x, f_int *incx, f_float *y, f_int *incy);
void daxpy_(f_int *n, f_double *alpha, f_double *x, f_int *incx, f_double *y, f_int *incy);
void caxpy_(f_int *n, f_cfloat *alpha, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
void zaxpy_(f_int *n, f_cdouble *alpha, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);


/// ret := x.T * y
float_ret_t sdot_(f_int *n, f_float *x, f_int *incx, f_float *y, f_int *incy);
f_double ddot_(f_int *n, f_double *x, f_int *incx, f_double *y, f_int *incy);
f_double dsdot_(f_int *n, f_float *sx, f_int *incx, f_float *sy, f_int *incy);
void cdotu_(f_cfloat *ret_val, f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
void zdotu_(f_cdouble *ret_val, f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);
//f_cfloat cdotu_(f_cfloat *ret_val, f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
//f_cdouble zdotu_(f_cdouble *ret_val, f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);

/// ret := x.H * y
void cdotc_(f_cfloat *ret_val, f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
void zdotc_(f_cdouble *ret_val, f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);
//f_cfloat cdotc_(f_cfloat *ret_val, f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
//f_cdouble zdotc_(f_cdouble *ret_val, f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);

/// ret := b + x.T * y
float_ret_t sdsdot_(f_int *n, f_float *b, f_float *x, f_int *incx, f_float *y, f_int *incy);

/// ret := sqrt( x.T * x )
float_ret_t scnrm2_(f_int *n, f_cfloat *x, f_int *incx);
float_ret_t snrm2_(f_int *n, f_float *x, f_int *incx);
f_double dnrm2_(f_int *n, f_double *x, f_int *incx);
f_double dznrm2_(f_int *n, f_cdouble *x, f_int *incx);

/// ret := |x|_1
float_ret_t sasum_(f_int *n, f_float *x, f_int *incx);
f_double dasum_(f_int *n, f_double *x, f_int *incx);


/// ret := |re(x)|_1 + |im(x)|_1
float_ret_t scasum_(f_int *n, f_cfloat *x, f_int *incx);
f_double dzasum_(f_int *n, f_cdouble *x, f_int *incx);

/// ret := argmax(abs(x_i))
f_int isamax_(f_int *n, f_float *x, f_int *incx);
f_int idamax_(f_int *n, f_double *x, f_int *incx);

/// ret := argmax( abs(re(x_i))+abs(im(x_i)) )
f_int icamax_(f_int *n, f_cfloat *x, f_int *incx);
f_int izamax_(f_int *n, f_cdouble *x, f_int *incx);


/// Level 2 BLAS

/** matrix vector multiply
        y = alpha*A*x + beta*y
   OR   y = alpha*A.T*x + beta*y
   OR   y = alpha*A.H*x + beta*y,  with A an mxn matrix
*/
void sgemv_(char *trans, f_int *m, f_int *n, f_float *alpha, f_float *A, f_int *lda, f_float *x, f_int *incx, f_float *beta, f_float *y, f_int *incy, f_int trans_len);
void dgemv_(char *trans, f_int *m, f_int *n, f_double *alpha, f_double *A, f_int *lda, f_double *x, f_int *incx, f_double *beta, f_double *y, f_int *incy, f_int trans_len);
void cgemv_(char *trans, f_int *m, f_int *n, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *x, f_int *incx, f_cfloat *beta, f_cfloat *y, f_int *incy, f_int trans_len);
void zgemv_(char *trans, f_int *m, f_int *n, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *x, f_int *incx, f_cdouble *beta, f_cdouble *y, f_int *incy, f_int trans_len);

/** banded matrix vector multiply
        y = alpha*A*x   + beta*y 
    OR  y = alpha*A.T*x + beta*y
    OR  y = alpha*A.H*x + beta*y,  with A a banded mxn matrix
*/
void sgbmv_(char *trans, f_int *m, f_int *n, f_int *kl, f_int *ku, f_float *alpha, f_float *A, f_int *lda, f_float *x, f_int *incx, f_float *beta, f_float *y, f_int *incy, f_int trans_len);
void dgbmv_(char *trans, f_int *m, f_int *n, f_int *kl, f_int *ku, f_double *alpha, f_double *A, f_int *lda, f_double *x, f_int *incx, f_double *beta, f_double *y, f_int *incy, f_int trans_len);
void cgbmv_(char *trans, f_int *m, f_int *n, f_int *kl, f_int *ku, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *x, f_int *incx, f_cfloat *beta, f_cfloat *y, f_int *incy, f_int trans_len);
void zgbmv_(char *trans, f_int *m, f_int *n, f_int *kl, f_int *ku, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *x, f_int *incx, f_cdouble *beta, f_cdouble *y, f_int *incy, f_int trans_len);

/** hermitian matrix vector multiply
 */
void chemv_(char *uplo, f_int *n, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *x, f_int *incx, f_cfloat *beta, f_cfloat *y, f_int *incy, f_int uplo_len);
void zhemv_(char *uplo, f_int *n, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *x, f_int *incx, f_cdouble *beta, f_cdouble *y, f_int *incy, f_int uplo_len);

/// hermitian banded matrix vector multiply
void chbmv_(char *uplo, f_int *n, f_int *k, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *x, f_int *incx, f_cfloat *beta, f_cfloat *y, f_int *incy, f_int uplo_len);
void zhbmv_(char *uplo, f_int *n, f_int *k, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *x, f_int *incx, f_cdouble *beta, f_cdouble *y, f_int *incy, f_int uplo_len);

/// hermitian packed matrix vector multiply
void chpmv_(char *uplo, f_int *n, f_cfloat *alpha, f_cfloat *A, f_cfloat *x, f_int *incx, f_cfloat *beta, f_cfloat *y, f_int *incy, f_int uplo_len);
void zhpmv_(char *uplo, f_int *n, f_cdouble *alpha, f_cdouble *A, f_cdouble *x, f_int *incx, f_cdouble *beta, f_cdouble *y, f_int *incy, f_int uplo_len);

/** symmetric matrix vector multiply
    y := alpha * A * x + beta * y
 */
void ssymv_(char *uplo, f_int *n, f_float *alpha, f_float *A, f_int *lda, f_float *x, f_int *incx, f_float *beta, f_float *y, f_int *incy, f_int uplo_len);
void dsymv_(char *uplo, f_int *n, f_double *alpha, f_double *A, f_int *lda, f_double *x, f_int *incx, f_double *beta, f_double *y, f_int *incy, f_int uplo_len);

/** symmetric banded matrix vector multiply
    y := alpha * A * x + beta * y
 */
void ssbmv_(char *uplo, f_int *n, f_int *k, f_float *alpha, f_float *A, f_int *lda, f_float *x, f_int *incx, f_float *beta, f_float *y, f_int *incy, f_int uplo_len);
void dsbmv_(char *uplo, f_int *n, f_int *k, f_double *alpha, f_double *A, f_int *lda, f_double *x, f_int *incx, f_double *beta, f_double *y, f_int *incy, f_int uplo_len);

/** symmetric packed matrix vector multiply
    y := alpha * A * x + beta * y
 */
void sspmv_(char *uplo, f_int *n, f_float *alpha, f_float *ap, f_float *x, f_int *incx, f_float *beta, f_float *y, f_int *incy, f_int uplo_len);
void dspmv_(char *uplo, f_int *n, f_double *alpha, f_double *ap, f_double *x, f_int *incx, f_double *beta, f_double *y, f_int *incy, f_int uplo_len);

/** triangular matrix vector multiply
        x := A * x
    OR  x := A.T * x
    OR  x := A.H * x
 */
void strmv_(char *uplo, char *trans, char *diag, f_int *n, f_float *A, f_int *lda, f_float *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtrmv_(char *uplo, char *trans, char *diag, f_int *n, f_double *A, f_int *lda, f_double *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctrmv_(char *uplo, char *trans, char *diag, f_int *n, f_cfloat *A, f_int *lda, f_cfloat *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztrmv_(char *uplo, char *trans, char *diag, f_int *n, f_cdouble *A, f_int *lda, f_cdouble *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);

/** triangular banded matrix vector multiply
        x := A * x
    OR  x := A.T * x
    OR  x := A.H * x
 */
void stbmv_(char *uplo, char *trans, char *diag, f_int *n, f_int *k, f_float *A, f_int *lda, f_float *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtbmv_(char *uplo, char *trans, char *diag, f_int *n, f_int *k, f_double *A, f_int *lda, f_double *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctbmv_(char *uplo, char *trans, char *diag, f_int *n, f_int *k, f_cfloat *A, f_int *lda, f_cfloat *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztbmv_(char *uplo, char *trans, char *diag, f_int *n, f_int *k, f_cdouble *A, f_int *lda, f_cdouble *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);

/** triangular packed matrix vector multiply
        x := A * x
    OR  x := A.T * x
    OR  x := A.H * x
 */
void stpmv_(char *uplo, char *trans, char *diag, f_int *n, f_float *ap, f_float *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtpmv_(char *uplo, char *trans, char *diag, f_int *n, f_double *ap, f_double *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctpmv_(char *uplo, char *trans, char *diag, f_int *n, f_cfloat *ap, f_cfloat *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztpmv_(char *uplo, char *trans, char *diag, f_int *n, f_cdouble *ap, f_cdouble *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);

/** solving triangular matrix problems
        x := A.inv * x
    OR  x := A.inv.T * x
    OR  x := A.inv.H * x
 */
void strsv_(char *uplo, char *trans, char *diag, f_int *n, f_float *A, f_int *lda, f_float *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtrsv_(char *uplo, char *trans, char *diag, f_int *n, f_double *A, f_int *lda, f_double *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctrsv_(char *uplo, char *trans, char *diag, f_int *n, f_cfloat *A, f_int *lda, f_cfloat *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztrsv_(char *uplo, char *trans, char *diag, f_int *n, f_cdouble *A, f_int *lda, f_cdouble *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);

/** solving triangular banded matrix problems
        x := A.inv * x
    OR  x := A.inv.T * x
    OR  x := A.inv.H * x
 */
void stbsv_(char *uplo, char *trans, char *diag, f_int *n, f_int *k, f_float *A, f_int *lda, f_float *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtbsv_(char *uplo, char *trans, char *diag, f_int *n, f_int *k, f_double *A, f_int *lda, f_double *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctbsv_(char *uplo, char *trans, char *diag, f_int *n, f_int *k, f_cfloat *A, f_int *lda, f_cfloat *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztbsv_(char *uplo, char *trans, char *diag, f_int *n, f_int *k, f_cdouble *A, f_int *lda, f_cdouble *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);

/** solving triangular packed matrix problems
        x := A.inv * x
    OR  x := A.inv.T * x
    OR  x := A.inv.H * x
 */
void stpsv_(char *uplo, char *trans, char *diag, f_int *n, f_float *ap, f_float *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtpsv_(char *uplo, char *trans, char *diag, f_int *n, f_double *ap, f_double *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctpsv_(char *uplo, char *trans, char *diag, f_int *n, f_cfloat *ap, f_cfloat *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztpsv_(char *uplo, char *trans, char *diag, f_int *n, f_cdouble *ap, f_cdouble *x, f_int *incx, f_int uplo_len, f_int trans_len, f_int diag_len);

/// performs the rank 1 operation 
///    A := A + alpha*x*y.T
void sger_(f_int *m, f_int *n, f_float *alpha, f_float *x, f_int *incx, f_float *y, f_int *incy, f_float *A, f_int *lda);
void dger_(f_int *m, f_int *n, f_double *alpha, f_double *x, f_int *incx, f_double *y, f_int *incy, f_double *A, f_int *lda);

/// performs the rank 1 operation 
///    A := A + alpha*x*y.T
void cgeru_(f_int *m, f_int *n, f_cfloat *alpha, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy, f_cfloat *A, f_int *lda);
void zgeru_(f_int *m, f_int *n, f_cdouble *alpha, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy, f_cdouble *A, f_int *lda);

/// performs the rank 1 operation 
///    A := A + alpha*x*y.H
void cgerc_(f_int *m, f_int *n, f_cfloat *alpha, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy, f_cfloat *A, f_int *lda);
void zgerc_(f_int *m, f_int *n, f_cdouble *alpha, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy, f_cdouble *A, f_int *lda);

/// hermitian rank 1 operation 
///    A := A + alpha*x*x.H
void cher_(char *uplo, f_int *n, f_float *alpha, f_cfloat *x, f_int *incx, f_cfloat *A, f_int *lda, f_int uplo_len);
void zher_(char *uplo, f_int *n, f_double *alpha, f_cdouble *x, f_int *incx, f_cdouble *A, f_int *lda, f_int uplo_len);

/// hermitian packed rank 1 operation
///    A := A + alpha*x*x.H
void chpr_(char *uplo, f_int *n, f_float *alpha, f_cfloat *x, f_int *incx, f_cfloat *A, f_int uplo_len);
void zhpr_(char *uplo, f_int *n, f_double *alpha, f_cdouble *x, f_int *incx, f_cdouble *A, f_int uplo_len);

/// hermitian rank 2 operation
///    A := A + alpha*x*y.H + alpha.conj * y * x.H
void cher2_(char *uplo, f_int *n, f_cfloat *alpha, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy, f_cfloat *A, f_int *lda, f_int uplo_len);
void zher2_(char *uplo, f_int *n, f_cdouble *alpha, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy, f_cdouble *A, f_int *lda, f_int uplo_len);

/// hermitian packed rank 2 operation
///    A := A + alpha*x*y.H + alpha.conj * y * x.H
void chpr2_(char *uplo, f_int *n, f_cfloat *alpha, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy, f_cfloat *A, f_int uplo_len);
void zhpr2_(char *uplo, f_int *n, f_cdouble *alpha, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy, f_cdouble *A, f_int uplo_len);

/// performs the symmetric rank 1 operation 
///    A := A + alpha*x*x.T
void ssyr_(char *uplo, f_int *n, f_float *alpha, f_float *x, f_int *incx, f_float *A, f_int *lda, f_int uplo_len);
void dsyr_(char *uplo, f_int *n, f_double *alpha, f_double *x, f_int *incx, f_double *A, f_int *lda, f_int uplo_len);

/// symmetric packed rank 1 operation  
///    A := A + alpha*x*x.T
void sspr_(char *uplo, f_int *n, f_float *alpha, f_float *x, f_int *incx, f_float *ap, f_int uplo_len);
void dspr_(char *uplo, f_int *n, f_double *alpha, f_double *x, f_int *incx, f_double *ap, f_int uplo_len);

/// performs the symmetric rank 2 operation
///    A := A + alpha * x * y.T  +  alpha * y * x.T
void ssyr2_(char *uplo, f_int *n, f_float *alpha, f_float *x, f_int *incx, f_float *y, f_int *incy, f_float *A, f_int *lda, f_int uplo_len);
void dsyr2_(char *uplo, f_int *n, f_double *alpha, f_double *x, f_int *incx, f_double *y, f_int *incy, f_double *A, f_int *lda, f_int uplo_len);

/// performs the symmetric packed rank 2 operation
///    A := A + alpha*x*y.T + alpha*y*x.T
void sspr2_(char *uplo, f_int *n, f_float *alpha, f_float *x, f_int *incx, f_float *y, f_int *incy, f_float *ap, f_int uplo_len);
void dspr2_(char *uplo, f_int *n, f_double *alpha, f_double *x, f_int *incx, f_double *y, f_int *incy, f_double *ap, f_int uplo_len);


/// Level 3 BLAS

/// matrix matrix multiply
///     C := alpha * transa(A) * transb(B) + beta * C
void sgemm_(char *transa, char *transb, f_int *m, f_int *n, f_int *k, f_float *alpha, f_float *A, f_int *lda, f_float *B, f_int *ldb, f_float *beta, f_float *C, f_int *ldc, f_int transa_len, f_int transb_len);
void dgemm_(char *transa, char *transb, f_int *m, f_int *n, f_int *k, f_double *alpha, f_double *A, f_int *lda, f_double *B, f_int *ldb, f_double *beta, f_double *C, f_int *ldc, f_int transa_len, f_int transb_len);
void cgemm_(char *transa, char *transb, f_int *m, f_int *n, f_int *k, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *B, f_int *ldb, f_cfloat *beta, f_cfloat *C, f_int *ldc, f_int transa_len, f_int transb_len);
void zgemm_(char *transa, char *transb, f_int *m, f_int *n, f_int *k, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *B, f_int *ldb, f_cdouble *beta, f_cdouble *C, f_int *ldc, f_int transa_len, f_int transb_len);

/// symmetric matrix matrix multiply
///     C := alpha * A * B + beta * C
/// OR  C := alpha * B * A + beta * C,    where A == A.T
void ssymm_(char *side, char *uplo, f_int *m, f_int *n, f_float *alpha, f_float *A, f_int *lda, f_float *B, f_int *ldb, f_float *beta, f_float *C, f_int *ldc, f_int side_len, f_int uplo_len);
void dsymm_(char *side, char *uplo, f_int *m, f_int *n, f_double *alpha, f_double *A, f_int *lda, f_double *B, f_int *ldb, f_double *beta, f_double *C, f_int *ldc, f_int side_len, f_int uplo_len);
void csymm_(char *side, char *uplo, f_int *m, f_int *n, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *B, f_int *ldb, f_cfloat *beta, f_cfloat *C, f_int *ldc, f_int side_len, f_int uplo_len);
void zsymm_(char *side, char *uplo, f_int *m, f_int *n, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *B, f_int *ldb, f_cdouble *beta, f_cdouble *C, f_int *ldc, f_int side_len, f_int uplo_len);

/// hermitian matrix matrix multiply
///     C := alpha * A * B + beta * C
/// OR  C := alpha * B * A + beta * C,    where A == A.H
void chemm_(char *side, char *uplo, f_int *m, f_int *n, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *B, f_int *ldb, f_cfloat *beta, f_cfloat *C, f_int *ldc, f_int side_len, f_int uplo_len);
void zhemm_(char *side, char *uplo, f_int *m, f_int *n, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *B, f_int *ldb, f_cdouble *beta, f_cdouble *C, f_int *ldc, f_int side_len, f_int uplo_len);

/// symmetric rank-k update to a matrix
///     C := alpha * A * A.T + beta * C
/// OR  C := alpha * A.T * A + beta * C
void ssyrk_(char *uplo, char *trans, f_int *n, f_int *k, f_float *alpha, f_float *A, f_int *lda, f_float *beta, f_float *C, f_int *ldc, f_int uplo_len, f_int trans_len);
void dsyrk_(char *uplo, char *trans, f_int *n, f_int *k, f_double *alpha, f_double *A, f_int *lda, f_double *beta, f_double *C, f_int *ldc, f_int uplo_len, f_int trans_len);
void csyrk_(char *uplo, char *trans, f_int *n, f_int *k, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *beta, f_cfloat *C, f_int *ldc, f_int uplo_len, f_int trans_len);
void zsyrk_(char *uplo, char *trans, f_int *n, f_int *k, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *beta, f_cdouble *C, f_int *ldc, f_int uplo_len, f_int trans_len);

/// hermitian rank-k update to a matrix
///     C := alpha * A * A.H + beta * C
/// OR  C := alpha * A.H * A + beta * C
void cherk_(char *uplo, char *trans, f_int *n, f_int *k, f_float *alpha, f_cfloat *A, f_int *lda, f_float *beta, f_cfloat *C, f_int *ldc, f_int uplo_len, f_int trans_len);
void zherk_(char *uplo, char *trans, f_int *n, f_int *k, f_double *alpha, f_cdouble *A, f_int *lda, f_double *beta, f_cdouble *C, f_int *ldc, f_int uplo_len, f_int trans_len);

/// symmetric rank-2k update to a matrix
///     C := alpha * A * B.T + alpha.conj * B * A.T + beta * C
/// OR  C := alpha * A.T * B + alpha.conj * B.T * A + beta * C
void ssyr2k_(char *uplo, char *trans, f_int *n, f_int *k, f_float *alpha, f_float *A, f_int *lda, f_float *B, f_int *ldb, f_float *beta, f_float *C, f_int *ldc, f_int uplo_len, f_int trans_len);
void dsyr2k_(char *uplo, char *trans, f_int *n, f_int *k, f_double *alpha, f_double *A, f_int *lda, f_double *B, f_int *ldb, f_double *beta, f_double *C, f_int *ldc, f_int uplo_len, f_int trans_len);
void csyr2k_(char *uplo, char *trans, f_int *n, f_int *k, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *B, f_int *ldb, f_cfloat *beta, f_cfloat *C, f_int *ldc, f_int uplo_len, f_int trans_len);
void zsyr2k_(char *uplo, char *trans, f_int *n, f_int *k, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *B, f_int *ldb, f_cdouble *beta, f_cdouble *C, f_int *ldc, f_int uplo_len, f_int trans_len);

/// hermitian rank-2k update to a matrix
///     C := alpha * A * B.H + alpha.conj * B * A.H + beta * C
/// OR  C := alpha * A.H * B + alpha.conj * B.H * A + beta * C
void cher2k_(char *uplo, char *trans, f_int *n, f_int *k, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *B, f_int *ldb, f_float *beta, f_cfloat *C, f_int *ldc, f_int uplo_len, f_int trans_len);
void zher2k_(char *uplo, char *trans, f_int *n, f_int *k, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *B, f_int *ldb, f_double *beta, f_cdouble *C, f_int *ldc, f_int uplo_len, f_int trans_len);

/// triangular matrix matrix multiply
///     B := alpha * transa(A) * B
/// OR  B := alpha * B * transa(A)
void strmm_(char *side, char *uplo, char *transa, char *diag, f_int *m, f_int *n, f_float *alpha, f_float *A, f_int *lda, f_float *B, f_int *ldb, f_int side_len, f_int uplo_len, f_int transa_len, f_int diag_len);
void dtrmm_(char *side, char *uplo, char *transa, char *diag, f_int *m, f_int *n, f_double *alpha, f_double *A, f_int *lda, f_double *B, f_int *ldb, f_int side_len, f_int uplo_len, f_int transa_len, f_int diag_len);
void ctrmm_(char *side, char *uplo, char *transa, char *diag, f_int *m, f_int *n, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *B, f_int *ldb, f_int side_len, f_int uplo_len, f_int transa_len, f_int diag_len);
void ztrmm_(char *side, char *uplo, char *transa, char *diag, f_int *m, f_int *n, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *B, f_int *ldb, f_int side_len, f_int uplo_len, f_int transa_len, f_int diag_len);

/// solving triangular matrix with multiple right hand sides
///     B := alpha * transa(A.inv) * B
/// OR  B := alpha * B * transa(A.inv)
void strsm_(char *side, char *uplo, char *transa, char *diag, f_int *m, f_int *n, f_float *alpha, f_float *A, f_int *lda, f_float *B, f_int *ldb, f_int side_len, f_int uplo_len, f_int transa_len, f_int diag_len);
void dtrsm_(char *side, char *uplo, char *transa, char *diag, f_int *m, f_int *n, f_double *alpha, f_double *A, f_int *lda, f_double *B, f_int *ldb, f_int side_len, f_int uplo_len, f_int transa_len, f_int diag_len);
void ctrsm_(char *side, char *uplo, char *transa, char *diag, f_int *m, f_int *n, f_cfloat *alpha, f_cfloat *A, f_int *lda, f_cfloat *B, f_int *ldb, f_int side_len, f_int uplo_len, f_int transa_len, f_int diag_len);
void ztrsm_(char *side, char *uplo, char *transa, char *diag, f_int *m, f_int *n, f_cdouble *alpha, f_cdouble *A, f_int *lda, f_cdouble *B, f_int *ldb, f_int side_len, f_int uplo_len, f_int transa_len, f_int diag_len);

/// Test if the characters are equal. (Auxiliary routine in Level 2 and 3 BLAS routines)
f_int lsame_(char *ca, char *cb, f_int ca_len, f_int cb_len);

/// Computes absolute values of a f_cdouble number. (Auxiliary routine for a few Level 1 BLAS routines)
f_double dcabs1_(f_cdouble *z);

/// Error handler for level 2 and 3 BLAS routines
void xerbla_(char *srname, f_int *info, f_int srname_len);

