/** BLAS bindings for D.

    Authors:    William V. Baxter III (with slight modifications by
                Lars Tandle Kyllingstad).
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.bindings.blas.dblas;


import scid.bindings.blas.blas;
public import scid.bindings.blas.types;
import scid.common.fortran;

// debug = blasCalls;

debug( blasCalls ) {
	import std.stdio, std.range, std.conv;
	private string stridedToString_( S )( const(S)* ptr, size_t len, size_t stride ) {
		if( len == 0 )
			return "[]";
	
		auto app = appender!string("[");
		app.put( to!string(*ptr) );
		auto e = ptr + len * stride;
		ptr += stride;
		for( ; ptr < e ; ptr += stride ) {
			app.put( ", " );
			app.put( to!string( *ptr ) );
		} 
		app.put(']');
		
		return app.data();
	}
}


/*  This unittest just checks whether the system is correctly
    set up to call BLAS routines.
*/
version(unittest)  import std.math: approxEqual;
unittest
{
    f_float[] fx = [1.0,2,3];
    f_float[] fy = [1.0,2,3];
    f_double[] dx = [1.0,2,3];
    f_double[] dy = [1.0,2,3];
    f_cfloat[] cx = [1.0+0i, 2+0i, 3+0i];
    f_cfloat[] cy = [1.0+0i, 2+0i, 3+0i];
    f_cdouble[] zx = [1.0+0i, 2+0i, 3+0i];
    f_cdouble[] zy = [1.0+0i, 2+0i, 3+0i];

    f_float fret = dot(toInt(fx.length), fx.ptr, 1,  fy.ptr, 1);
    assert (approxEqual(fret, 14.0, 1e-7));

    f_double dret = dot(toInt(dx.length), dx.ptr, 1,  dy.ptr, 1);
    assert (approxEqual(dret, 14.0, 1e-7));

    // This test causes segmentation fault, for reasons unknown.
    /*
    f_cfloat cret = dotu(toInt(cx.length), cx.ptr, 1,  cy.ptr, 1);
    assert (approxEqual(cret.re, 14.0, 1e-7)
            && approxEqual(cret.im, 0.0, 1e-7));
    */


    // This test fails on 64-bit, also for unknown reason.
    version(X86)
    {
    f_cdouble zret = dotu(toInt(zx.length), zx.ptr, 1,  zy.ptr, 1);
    assert (approxEqual(zret.re, 14.0, 1e-7)
            && approxEqual(zret.im, 0.0, 1e-7));
    }
}


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
void rotg(ref f_float a, ref f_float b, out f_float c, out f_float s) {
    srotg_(&a, &b, &c, &s);
}
void rotg(ref f_double a, ref f_double b, out f_double c, out f_double s) {
    drotg_(&a, &b, &c, &s);
}
void rotg(ref f_cfloat a, ref f_cfloat b, out f_float c, out f_cfloat s) {
    crotg_(&a, &b, &c, &s);
}
void rotg(ref f_cdouble a, ref f_cdouble b, out f_double c, out f_cdouble s) {
    zrotg_(&a, &b, &c, &s);
}

/// Generate modified plane (Givens) rotation
void rotmg(ref f_double d1, ref f_double d2, ref f_double b1, ref f_double b2, f_double *param) {
    drotmg_(&d1, &d2, &b1, &b2, param);
}
void rotmg(ref f_float d1, ref f_float d2, ref f_float b1, ref f_float b2, f_float *param) {
    srotmg_(&d1, &d2, &b1, &b2, param);
}

/// Apply plane (Givens) rotation
///             _      _
///     x_i  := | c  s | * x_i
///     y_i     |-s  c |   y_i
///             -      -
void rot(f_int n, f_float *x, f_int incx, f_float *y, f_int incy, f_float c, f_float s) {
    srot_(&n, x, &incx, y, &incy, &c, &s);
}
void rot(f_int n, f_double *x, f_int incx, f_double *y, f_int incy, f_double c, f_double s) {
    drot_(&n, x, &incx, y, &incy, &c, &s);
}
void rot(f_int n, f_cfloat *x, f_int incx, f_cfloat *y, f_int incy, f_float c, f_float s) {
    csrot_(&n, x, &incx, y, &incy, &c, &s);
}
void rot(f_int n, f_cdouble *x, f_int incx, f_cdouble *y, f_int incy, f_double c, f_double s) {
    zdrot_(&n, x, &incx, y, &incy, &c, &s);
}

/// Apply modified plane (Givens) rotation
void rotm(f_int n, f_float *x, f_int incx, f_float *y, f_int incy, f_float *param) {
    srotm_(&n, x, &incx, y, &incy, param);
}
void rotm(f_int n, f_double *x, f_int incx, f_double *y, f_int incy, f_double *param) {
    drotm_(&n, x, &incx, y, &incy, param);
}

/// Swap the values contained in x and y 
///     x <-> y
void swap(f_int n, f_float *x, f_int incx, f_float *y, f_int incy) {
    sswap_(&n, x, &incx, y, &incy);
}
void swap(f_int n, f_double *x, f_int incx, f_double *y, f_int incy) {
    dswap_(&n, x, &incx, y, &incy);
}
void swap(f_int n, f_cfloat *x, f_int incx, f_cfloat *y, f_int incy) {
    cswap_(&n, x, &incx, y, &incy);
}
void swap(f_int n, f_cdouble *x, f_int incx, f_cdouble *y, f_int incy) {
    zswap_(&n, x, &incx, y, &incy);
}

/// x := alpha * x
void scal(f_int n, f_float alpha, f_float *x, f_int incx) {
	debug( blasCalls ) write( "scal( ", alpha, ", ", stridedToString_(x, n, incx), " ) => " );
    sscal_(&n, &alpha, x, &incx);
	debug( blasCalls ) writeln( stridedToString_(x, n, incx) );
}
void scal(f_int n, f_double alpha, f_double *x, f_int incx) {
	debug( blasCalls ) write( "scal( ", alpha, ", ", stridedToString_(x, n, incx), " ) => " );
    dscal_(&n, &alpha, x, &incx);
	debug( blasCalls ) writeln( stridedToString_(x, n, incx) );
}
void scal(f_int n, f_cfloat alpha, f_cfloat *x, f_int incx) {
	debug( blasCalls ) write( "scal( ", alpha, ", ", stridedToString_(x, n, incx), " ) => " );
    cscal_(&n, &alpha, x, &incx);
	debug( blasCalls ) writeln( stridedToString_(x, n, incx) );
}
void scal(f_int n, f_float alpha, f_cfloat *x, f_int incx) {
	debug( blasCalls ) write( "scal( ", alpha, ", ", stridedToString_(x, n, incx), " ) => " );
    csscal_(&n, &alpha, x, &incx);
	debug( blasCalls ) writeln( stridedToString_(x, n, incx) );
}
void scal(f_int n, f_cdouble alpha, f_cdouble *x, f_int incx) {
	debug( blasCalls ) write( "scal( ", alpha, ", ", stridedToString_(x, n, incx), " ) => " );
    zscal_(&n, &alpha, x, &incx);
	debug( blasCalls ) writeln( stridedToString_(x, n, incx) );
}
void scal(f_int n, f_double alpha, f_cdouble *x, f_int incx) {
	debug( blasCalls ) write( "scal( ", alpha, ", ", stridedToString_(x, n, incx), " ) => " );
    zdscal_(&n, &alpha, x, &incx);
	debug( blasCalls ) writeln( stridedToString_(x, n, incx) );
}

/// y := x
void copy(f_int n, const f_float *x, f_int incx, f_float *y, f_int incy) {
	debug( blasCalls ) write( "copy( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    scopy_(&n, cast(f_float*) x, &incx, y, &incy);
	debug( blasCalls ) writeln( stridedToString_(y, n, incy) );
}
void copy(f_int n, const f_double *x, f_int incx, f_double *y, f_int incy) {
	debug( blasCalls ) write( "copy( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    dcopy_(&n, cast(f_double*) x, &incx, y, &incy);
	debug( blasCalls ) writeln( stridedToString_(y, n, incy) );
}
void copy(f_int n, const f_cfloat *x, f_int incx, f_cfloat *y, f_int incy) {
	debug( blasCalls ) write( "copy( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    ccopy_(&n, cast(f_cfloat*) x, &incx, y, &incy);
	debug( blasCalls ) writeln( stridedToString_(y, n, incy) );
}
void copy(f_int n, const f_cdouble *x, f_int incx, f_cdouble *y, f_int incy) {
	debug( blasCalls ) write( "copy( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    zcopy_(&n, cast(f_cdouble*) x, &incx, y, &incy);
	debug( blasCalls ) writeln( stridedToString_(y, n, incy) );
}

/// y := alpha * x + y
void axpy(f_int n, f_float alpha, const f_float *x, f_int incx, f_float *y, f_int incy) {
	debug( blasCalls ) write( "axpy( ", alpha, ", ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    saxpy_(&n, &alpha, cast(f_float*)x, &incx, y, &incy);
	debug( blasCalls ) writeln( stridedToString_(y, n, incy) );
}
void axpy(f_int n, f_double alpha, const f_double *x, f_int incx, f_double *y, f_int incy) {
	debug( blasCalls ) write( "axpy( ", alpha, ", ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    daxpy_(&n, &alpha, cast(f_double*)x, &incx, y, &incy);
	debug( blasCalls ) writeln( stridedToString_(y, n, incy) );
}
void axpy(f_int n, f_cfloat alpha, const f_cfloat *x, f_int incx, f_cfloat *y, f_int incy) {
	debug( blasCalls ) write( "axpy( ", alpha, ", ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    caxpy_(&n, &alpha, cast(f_cfloat*)x, &incx, y, &incy);
	debug( blasCalls ) writeln( stridedToString_(y, n, incy) );
}
void axpy(f_int n, f_cdouble alpha, const f_cdouble *x, f_int incx, f_cdouble *y, f_int incy) {
	debug( blasCalls ) write( "axpy( ", alpha, ", ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    zaxpy_(&n, &alpha, cast(f_cdouble*)x, &incx, y, &incy);
	debug( blasCalls ) writeln( stridedToString_(y, n, incy) );
}


/// ret := x.T * y
float_ret_t dot(f_int n, const f_float *x, f_int incx, const f_float *y, f_int incy) {
	debug( blasCalls ) write( "dot( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    auto ret_val = sdot_(&n, cast(f_float*)x, &incx, cast(f_float*)y, &incy);
	debug( blasCalls ) writeln( ret_val );
	return ret_val;
}
f_double dot(f_int n, const f_double *x, f_int incx, const f_double *y, f_int incy) {
	debug( blasCalls ) write( "ddot( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    auto ret_val =  ddot_(&n, cast(f_double*)x, &incx, cast(f_double*)y, &incy);
	debug( blasCalls ) writeln( ret_val );
	return ret_val;
}
f_double ddot(f_int n, const f_float *sx, f_int incx, const f_float *sy, f_int incy) {
    debug( blasCalls ) write( "dsdot( ", stridedToString_(sx, n, incx), ", ", stridedToString_(sy, n, incy), " ) => " );
	auto ret_val =  dsdot_(&n, cast(f_float*)sx, &incx, cast(f_float*)sy, &incy);
	debug( blasCalls ) writeln( stridedToString_(sy, n, incy) );
	return ret_val;
}
f_cfloat dotu(f_int n, const f_cfloat *x, f_int incx, const f_cfloat *y, f_int incy) {
	debug( blasCalls ) write( "dotu( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    f_cfloat ret_val;
    cdotu_(&ret_val, &n, cast(f_cfloat*)x, &incx,cast(f_cfloat*) y, &incy);
	debug( blasCalls ) writeln( ret_val );
    return ret_val;
	
}
f_cdouble dotu(f_int n, const f_cdouble *x, f_int incx, const f_cdouble *y, f_int incy) {
    debug( blasCalls ) write( "dotu( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
	f_cdouble ret_val;
    zdotu_(&ret_val, &n, cast(f_cdouble*)x, &incx, cast(f_cdouble*)y, &incy);
	debug( blasCalls ) writeln( ret_val );
    return ret_val;
}
//f_cfloat cdotu_(f_cfloat *ret_val, f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
//f_cdouble zdotu_(f_cdouble *ret_val, f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);

/// ret := x.H * y
f_cfloat dotc(f_int n, const f_cfloat *x, f_int incx, const f_cfloat *y, f_int incy) {
	debug( blasCalls ) write( "dotc( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    f_cfloat ret_val;
    cdotc_(&ret_val, &n, cast(f_cfloat*)x,  &incx, cast(f_cfloat*)y, &incy);
	debug( blasCalls ) writeln( ret_val );
    return ret_val;
}
f_cdouble dotc(f_int n, const f_cdouble *x, f_int incx, const f_cdouble *y, f_int incy) {
	debug( blasCalls ) write( "dotc( ", stridedToString_(x, n, incx), ", ", stridedToString_(y, n, incy), " ) => " );
    f_cdouble ret_val;
    zdotc_(&ret_val, &n, cast(f_cdouble*)x, &incx, cast(f_cdouble*)y, &incy);
	debug( blasCalls ) writeln( ret_val );
    return ret_val;
}
//f_cfloat cdotc_(f_cfloat *ret_val, f_int *n, f_cfloat *x, f_int *incx, f_cfloat *y, f_int *incy);
//f_cdouble zdotc_(f_cdouble *ret_val, f_int *n, f_cdouble *x, f_int *incx, f_cdouble *y, f_int *incy);

/// ret := b + x.T * y
float_ret_t dsdot(f_int n, f_float *b, f_float *x, f_int incx, f_float *y, f_int incy) {
    return sdsdot_(&n, b, x, &incx, y, &incy);
}

/// ret := sqrt( x.T * x )
float_ret_t nrm2(f_int n, f_cfloat *x, f_int incx) {
    return scnrm2_(&n, x, &incx);
}
float_ret_t nrm2(f_int n, f_float *x, f_int incx) {
    return snrm2_(&n, x, &incx);
}
f_double nrm2(f_int n, f_double *x, f_int incx) {
    return dnrm2_(&n, x, &incx);
}
f_double nrm2(f_int n, f_cdouble *x, f_int incx) {
    return dznrm2_(&n, x, &incx);
}

/// ret := |x|_1
float_ret_t asum(f_int n, f_float *x, f_int incx) {
    return sasum_(&n, x, &incx);
}
f_double asum(f_int n, f_double *x, f_int incx) {
    return dasum_(&n, x, &incx);
}


/// ret := |re(x)|_1 + |im(x)|_1
float_ret_t asum(f_int n, f_cfloat *x, f_int incx) {
    return scasum_(&n, x, &incx);
}
f_double asum(f_int n, f_cdouble *x, f_int incx) {
    return dzasum_(&n, x, &incx);
}

/// ret := argmax(abs(x_i))
f_int isamax(f_int n, f_float *x, f_int incx) {
    return isamax_(&n, x, &incx);
}
f_int idamax(f_int n, f_double *x, f_int incx) {
    return idamax_(&n, x, &incx);
}

/// ret := argmax( abs(re(x_i))+abs(im(x_i)) )
f_int icamax(f_int n, f_cfloat *x, f_int incx) {
    return icamax_(&n, x, &incx);
}
f_int izamax(f_int n, f_cdouble *x, f_int incx) {
    return izamax_(&n, x, &incx);
}


/// Level 2 BLAS

/** matrix vector multiply
        y = alpha*A*x + beta*y
   OR   y = alpha*A.T*x + beta*y
   OR   y = alpha*A.H*x + beta*y,  with A an mxn matrix
*/
void gemv(char trans, f_int m, f_int n, f_float alpha, f_float *A, f_int lda, f_float *x, f_int incx, f_float beta, f_float *y, f_int incy) {
    sgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void gemv(char trans, f_int m, f_int n, f_double alpha, f_double *A, f_int lda, f_double *x, f_int incx, f_double beta, f_double *y, f_int incy) {
    dgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void gemv(char trans, f_int m, f_int n, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *x, f_int incx, f_cfloat beta, f_cfloat *y, f_int incy) {
    cgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void gemv(char trans, f_int m, f_int n, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *x, f_int incx, f_cdouble beta, f_cdouble *y, f_int incy) {
    zgemv_(&trans, &m, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}

/** banded matrix vector multiply
        y = alpha*A*x   + beta*y 
    OR  y = alpha*A.T*x + beta*y
    OR  y = alpha*A.H*x + beta*y,  with A a banded mxn matrix
*/
void gbmv(char trans, f_int m, f_int n, f_int kl, f_int ku, f_float alpha, f_float *A, f_int lda, f_float *x, f_int incx, f_float beta, f_float *y, f_int incy) {
    sgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void gbmv(char trans, f_int m, f_int n, f_int kl, f_int ku, f_double alpha, f_double *A, f_int lda, f_double *x, f_int incx, f_double beta, f_double *y, f_int incy) {
    dgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void gbmv(char trans, f_int m, f_int n, f_int kl, f_int ku, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *x, f_int incx, f_cfloat beta, f_cfloat *y, f_int incy) {
    cgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void gbmv(char trans, f_int m, f_int n, f_int kl, f_int ku, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *x, f_int incx, f_cdouble beta, f_cdouble *y, f_int incy) {
    zgbmv_(&trans, &m, &n, &kl, &ku, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}

/** hermitian matrix vector multiply
 */
void hemv(char uplo, f_int n, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *x, f_int incx, f_cfloat beta, f_cfloat *y, f_int incy) {
    chemv_(&uplo, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void hemv(char uplo, f_int n, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *x, f_int incx, f_cdouble beta, f_cdouble *y, f_int incy) {
    zhemv_(&uplo, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}

/// hermitian banded matrix vector multiply
void hbmv(char uplo, f_int n, f_int k, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *x, f_int incx, f_cfloat beta, f_cfloat *y, f_int incy) {
    chbmv_(&uplo, &n, &k, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void hbmv(char uplo, f_int n, f_int k, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *x, f_int incx, f_cdouble beta, f_cdouble *y, f_int incy) {
    zhbmv_(&uplo, &n, &k, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}

/// hermitian packed matrix vector multiply
void hpmv(char uplo, f_int n, f_cfloat alpha, f_cfloat *A, f_cfloat *x, f_int incx, f_cfloat beta, f_cfloat *y, f_int incy) {
    chpmv_(&uplo, &n, &alpha, A, x, &incx, &beta, y, &incy, 1);
}
void hpmv(char uplo, f_int n, f_cdouble alpha, f_cdouble *A, f_cdouble *x, f_int incx, f_cdouble beta, f_cdouble *y, f_int incy) {
    zhpmv_(&uplo, &n, &alpha, A, x, &incx, &beta, y, &incy, 1);
}

/** symmetric matrix vector multiply
    y := alpha * A * x + beta * y
 */
void symv(char uplo, f_int n, f_float alpha, f_float *A, f_int lda, f_float *x, f_int incx, f_float beta, f_float *y, f_int incy) {
    ssymv_(&uplo, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void symv(char uplo, f_int n, f_double alpha, f_double *A, f_int lda, f_double *x, f_int incx, f_double beta, f_double *y, f_int incy) {
    dsymv_(&uplo, &n, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}

/** symmetric banded matrix vector multiply
    y := alpha * A * x + beta * y
 */
void sbmv(char uplo, f_int n, f_int k, f_float alpha, f_float *A, f_int lda, f_float *x, f_int incx, f_float beta, f_float *y, f_int incy) {
    ssbmv_(&uplo, &n, &k, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}
void sbmv(char uplo, f_int n, f_int k, f_double alpha, f_double *A, f_int lda, f_double *x, f_int incx, f_double beta, f_double *y, f_int incy) {
    dsbmv_(&uplo, &n, &k, &alpha, A, &lda, x, &incx, &beta, y, &incy, 1);
}

/** symmetric packed matrix vector multiply
    y := alpha * A * x + beta * y
 */
void spmv(char uplo, f_int n, f_float alpha, f_float *ap, f_float *x, f_int incx, f_float beta, f_float *y, f_int incy) {
    sspmv_(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy, 1);
}
void spmv(char uplo, f_int n, f_double alpha, f_double *ap, f_double *x, f_int incx, f_double beta, f_double *y, f_int incy) {
    dspmv_(&uplo, &n, &alpha, ap, x, &incx, &beta, y, &incy, 1);
}

/** triangular matrix vector multiply
        x := A * x
    OR  x := A.T * x
    OR  x := A.H * x
 */
void trmv(char uplo, char trans, char diag, f_int n, f_float *A, f_int lda, f_float *x, f_int incx) {
    strmv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx, 1, 1, 1);
}
void trmv(char uplo, char trans, char diag, f_int n, f_double *A, f_int lda, f_double *x, f_int incx) {
    dtrmv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx, 1, 1, 1);
}
void trmv(char uplo, char trans, char diag, f_int n, f_cfloat *A, f_int lda, f_cfloat *x, f_int incx) {
    ctrmv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx, 1, 1, 1);
}
void trmv(char uplo, char trans, char diag, f_int n, f_cdouble *A, f_int lda, f_cdouble *x, f_int incx) {
    ztrmv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx, 1, 1, 1);
}

/** triangular banded matrix vector multiply
        x := A * x
    OR  x := A.T * x
    OR  x := A.H * x
 */
void tbmv(char uplo, char trans, char diag, f_int n, f_int k, f_float *A, f_int lda, f_float *x, f_int incx) {
    stbmv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx, 1, 1, 1);
}
void tbmv(char uplo, char trans, char diag, f_int n, f_int k, f_double *A, f_int lda, f_double *x, f_int incx) {
    dtbmv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx, 1, 1, 1);
}
void tbmv(char uplo, char trans, char diag, f_int n, f_int k, f_cfloat *A, f_int lda, f_cfloat *x, f_int incx) {
    ctbmv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx, 1, 1, 1);
}
void tbmv(char uplo, char trans, char diag, f_int n, f_int k, f_cdouble *A, f_int lda, f_cdouble *x, f_int incx) {
    ztbmv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx, 1, 1, 1);
}

/** triangular packed matrix vector multiply
        x := A * x
    OR  x := A.T * x
    OR  x := A.H * x
 */
void tpmv(char uplo, char trans, char diag, f_int n, f_float *ap, f_float *x, f_int incx) {
    stpmv_(&uplo, &trans, &diag, &n, ap, x, &incx, 1, 1, 1);
}
void tpmv(char uplo, char trans, char diag, f_int n, f_double *ap, f_double *x, f_int incx) {
    dtpmv_(&uplo, &trans, &diag, &n, ap, x, &incx, 1, 1, 1);
}
void tpmv(char uplo, char trans, char diag, f_int n, f_cfloat *ap, f_cfloat *x, f_int incx) {
    ctpmv_(&uplo, &trans, &diag, &n, ap, x, &incx, 1, 1, 1);
}
void tpmv(char uplo, char trans, char diag, f_int n, f_cdouble *ap, f_cdouble *x, f_int incx) {
    ztpmv_(&uplo, &trans, &diag, &n, ap, x, &incx, 1, 1, 1);
}

/** solving triangular matrix problems
        x := A.inv * x
    OR  x := A.inv.T * x
    OR  x := A.inv.H * x
 */
void trsv(char uplo, char trans, char diag, f_int n, f_float *A, f_int lda, f_float *x, f_int incx) {
    strsv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx, 1, 1, 1);
}
void trsv(char uplo, char trans, char diag, f_int n, f_double *A, f_int lda, f_double *x, f_int incx) {
    dtrsv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx, 1, 1, 1);
}
void trsv(char uplo, char trans, char diag, f_int n, f_cfloat *A, f_int lda, f_cfloat *x, f_int incx) {
    ctrsv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx, 1, 1, 1);
}
void trsv(char uplo, char trans, char diag, f_int n, f_cdouble *A, f_int lda, f_cdouble *x, f_int incx) {
    ztrsv_(&uplo, &trans, &diag, &n, A, &lda, x, &incx, 1, 1, 1);
}

/** solving triangular banded matrix problems
        x := A.inv * x
    OR  x := A.inv.T * x
    OR  x := A.inv.H * x
 */
void tbsv(char uplo, char trans, char diag, f_int n, f_int k, f_float *A, f_int lda, f_float *x, f_int incx) {
    stbsv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx, 1, 1, 1);
}
void tbsv(char uplo, char trans, char diag, f_int n, f_int k, f_double *A, f_int lda, f_double *x, f_int incx) {
    dtbsv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx, 1, 1, 1);
}
void tbsv(char uplo, char trans, char diag, f_int n, f_int k, f_cfloat *A, f_int lda, f_cfloat *x, f_int incx) {
    ctbsv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx, 1, 1, 1);
}
void tbsv(char uplo, char trans, char diag, f_int n, f_int k, f_cdouble *A, f_int lda, f_cdouble *x, f_int incx) {
    ztbsv_(&uplo, &trans, &diag, &n, &k, A, &lda, x, &incx, 1, 1, 1);
}

/** solving triangular packed matrix problems
        x := A.inv * x
    OR  x := A.inv.T * x
    OR  x := A.inv.H * x
 */
void tpsv(char uplo, char trans, char diag, f_int n, f_float *ap, f_float *x, f_int incx) {
    stpsv_(&uplo, &trans, &diag, &n, ap, x, &incx, 1, 1, 1);
}
void tpsv(char uplo, char trans, char diag, f_int n, f_double *ap, f_double *x, f_int incx) {
    dtpsv_(&uplo, &trans, &diag, &n, ap, x, &incx, 1, 1, 1);
}
void tpsv(char uplo, char trans, char diag, f_int n, f_cfloat *ap, f_cfloat *x, f_int incx) {
    ctpsv_(&uplo, &trans, &diag, &n, ap, x, &incx, 1, 1, 1);
}
void tpsv(char uplo, char trans, char diag, f_int n, f_cdouble *ap, f_cdouble *x, f_int incx) {
    ztpsv_(&uplo, &trans, &diag, &n, ap, x, &incx, 1, 1, 1);
}

/// performs the rank 1 operation 
///    A := A + alpha*x*y.T
void ger(f_int m, f_int n, f_float alpha, f_float *x, f_int incx, f_float *y, f_int incy, f_float *A, f_int lda) {
    sger_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}
void ger(f_int m, f_int n, f_double alpha, f_double *x, f_int incx, f_double *y, f_int incy, f_double *A, f_int lda) {
    dger_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}

/// performs the rank 1 operation 
///    A := A + alpha*x*y.T
void geru(f_int m, f_int n, f_cfloat alpha, f_cfloat *x, f_int incx, f_cfloat *y, f_int incy, f_cfloat *A, f_int lda) {
    cgeru_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}
void geru(f_int m, f_int n, f_cdouble alpha, f_cdouble *x, f_int incx, f_cdouble *y, f_int incy, f_cdouble *A, f_int lda) {
    zgeru_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}

/// performs the rank 1 operation 
///    A := A + alpha*x*y.H
void gerc(f_int m, f_int n, f_cfloat alpha, f_cfloat *x, f_int incx, f_cfloat *y, f_int incy, f_cfloat *A, f_int lda) {
    cgerc_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}
void gerc(f_int m, f_int n, f_cdouble alpha, f_cdouble *x, f_int incx, f_cdouble *y, f_int incy, f_cdouble *A, f_int lda) {
    zgerc_(&m, &n, &alpha, x, &incx, y, &incy, A, &lda);
}

/// hermitian rank 1 operation 
///    A := A + alpha*x*x.H
void her(char uplo, f_int n, f_float alpha, f_cfloat *x, f_int incx, f_cfloat *A, f_int lda) {
    cher_(&uplo, &n, &alpha, x, &incx, A, &lda, 1);
}
void her(char uplo, f_int n, f_double alpha, f_cdouble *x, f_int incx, f_cdouble *A, f_int lda) {
    zher_(&uplo, &n, &alpha, x, &incx, A, &lda, 1);
}

/// hermitian packed rank 1 operation
///    A := A + alpha*x*x.H
void hpr(char uplo, f_int n, f_float alpha, f_cfloat *x, f_int incx, f_cfloat *A) {
    chpr_(&uplo, &n, &alpha, x, &incx, A, 1);
}
void hpr(char uplo, f_int n, f_double alpha, f_cdouble *x, f_int incx, f_cdouble *A) {
    zhpr_(&uplo, &n, &alpha, x, &incx, A, 1);
}

/// hermitian rank 2 operation
///    A := A + alpha*x*y.H + alpha.conj * y * x.H
void her2(char uplo, f_int n, f_cfloat alpha, f_cfloat *x, f_int incx, f_cfloat *y, f_int incy, f_cfloat *A, f_int lda) {
    cher2_(&uplo, &n, &alpha, x, &incx, y, &incy, A, &lda, 1);
}
void her2(char uplo, f_int n, f_cdouble alpha, f_cdouble *x, f_int incx, f_cdouble *y, f_int incy, f_cdouble *A, f_int lda) {
    zher2_(&uplo, &n, &alpha, x, &incx, y, &incy, A, &lda, 1);
}

/// hermitian packed rank 2 operation
///    A := A + alpha*x*y.H + alpha.conj * y * x.H
void hpr2(char uplo, f_int n, f_cfloat alpha, f_cfloat *x, f_int incx, f_cfloat *y, f_int incy, f_cfloat *A) {
    chpr2_(&uplo, &n, &alpha, x, &incx, y, &incy, A, 1);
}
void hpr2(char uplo, f_int n, f_cdouble alpha, f_cdouble *x, f_int incx, f_cdouble *y, f_int incy, f_cdouble *A) {
    zhpr2_(&uplo, &n, &alpha, x, &incx, y, &incy, A, 1);
}

/// performs the symmetric rank 1 operation 
///    A := A + alpha*x*x.T
void syr(char uplo, f_int n, f_float alpha, f_float *x, f_int incx, f_float *A, f_int lda) {
    ssyr_(&uplo, &n, &alpha, x, &incx, A, &lda, 1);
}
void syr(char uplo, f_int n, f_double alpha, f_double *x, f_int incx, f_double *A, f_int lda) {
    dsyr_(&uplo, &n, &alpha, x, &incx, A, &lda, 1);
}

/// symmetric packed rank 1 operation  
///    A := A + alpha*x*x.T
void spr(char uplo, f_int n, f_float alpha, f_float *x, f_int incx, f_float *ap) {
    sspr_(&uplo, &n, &alpha, x, &incx, ap, 1);
}
void spr(char uplo, f_int n, f_double alpha, f_double *x, f_int incx, f_double *ap) {
    dspr_(&uplo, &n, &alpha, x, &incx, ap, 1);
}

/// performs the symmetric rank 2 operation
///    A := A + alpha * x * y.T  +  alpha * y * x.T
void syr2(char uplo, f_int n, f_float alpha, f_float *x, f_int incx, f_float *y, f_int incy, f_float *A, f_int lda) {
    ssyr2_(&uplo, &n, &alpha, x, &incx, y, &incy, A, &lda, 1);
}
void syr2(char uplo, f_int n, f_double alpha, f_double *x, f_int incx, f_double *y, f_int incy, f_double *A, f_int lda) {
    dsyr2_(&uplo, &n, &alpha, x, &incx, y, &incy, A, &lda, 1);
}

/// performs the symmetric packed rank 2 operation
///    A := A + alpha*x*y.T + alpha*y*x.T
void spr2(char uplo, f_int n, f_float alpha, f_float *x, f_int incx, f_float *y, f_int incy, f_float *ap) {
    sspr2_(&uplo, &n, &alpha, x, &incx, y, &incy, ap, 1);
}
void spr2(char uplo, f_int n, f_double alpha, f_double *x, f_int incx, f_double *y, f_int incy, f_double *ap) {
    dspr2_(&uplo, &n, &alpha, x, &incx, y, &incy, ap, 1);
}


/// Level 3 BLAS

/// matrix matrix multiply
///     C := alpha * transa(A) * transb(B) + beta * C
void gemm(char transa, char transb, f_int m, f_int n, f_int k, f_float alpha, f_float *A, f_int lda, f_float *B, f_int ldb, f_float beta, f_float *C, f_int ldc) {
    sgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void gemm(char transa, char transb, f_int m, f_int n, f_int k, f_double alpha, f_double *A, f_int lda, f_double *B, f_int ldb, f_double beta, f_double *C, f_int ldc) {
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void gemm(char transa, char transb, f_int m, f_int n, f_int k, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *B, f_int ldb, f_cfloat beta, f_cfloat *C, f_int ldc) {
    cgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void gemm(char transa, char transb, f_int m, f_int n, f_int k, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *B, f_int ldb, f_cdouble beta, f_cdouble *C, f_int ldc) {
    zgemm_(&transa, &transb, &m, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}

/// symmetric matrix matrix multiply
///     C := alpha * A * B + beta * C
/// OR  C := alpha * B * A + beta * C,    where A == A.T
void symm(char side, char uplo, f_int m, f_int n, f_float alpha, f_float *A, f_int lda, f_float *B, f_int ldb, f_float beta, f_float *C, f_int ldc) {
    ssymm_(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void symm(char side, char uplo, f_int m, f_int n, f_double alpha, f_double *A, f_int lda, f_double *B, f_int ldb, f_double beta, f_double *C, f_int ldc) {
    dsymm_(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void symm(char side, char uplo, f_int m, f_int n, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *B, f_int ldb, f_cfloat beta, f_cfloat *C, f_int ldc) {
    csymm_(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void symm(char side, char uplo, f_int m, f_int n, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *B, f_int ldb, f_cdouble beta, f_cdouble *C, f_int ldc) {
    zsymm_(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}

/// hermitian matrix matrix multiply
///     C := alpha * A * B + beta * C
/// OR  C := alpha * B * A + beta * C,    where A == A.H
void hemm(char side, char uplo, f_int m, f_int n, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *B, f_int ldb, f_cfloat beta, f_cfloat *C, f_int ldc) {
    chemm_(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void hemm(char side, char uplo, f_int m, f_int n, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *B, f_int ldb, f_cdouble beta, f_cdouble *C, f_int ldc) {
    zhemm_(&side, &uplo, &m, &n, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}

/// symmetric rank-k update to a matrix
///     C := alpha * A * A.T + beta * C
/// OR  C := alpha * A.T * A + beta * C
void syrk(char uplo, char trans, f_int n, f_int k, f_float alpha, f_float *A, f_int lda, f_float beta, f_float *C, f_int ldc) {
    ssyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc, 1, 1);
}
void syrk(char uplo, char trans, f_int n, f_int k, f_double alpha, f_double *A, f_int lda, f_double beta, f_double *C, f_int ldc) {
    dsyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc, 1, 1);
}
void syrk(char uplo, char trans, f_int n, f_int k, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat beta, f_cfloat *C, f_int ldc) {
    csyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc, 1, 1);
}
void syrk(char uplo, char trans, f_int n, f_int k, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble beta, f_cdouble *C, f_int ldc) {
    zsyrk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc, 1, 1);
}

/// hermitian rank-k update to a matrix
///     C := alpha * A * A.H + beta * C
/// OR  C := alpha * A.H * A + beta * C
void herk(char uplo, char trans, f_int n, f_int k, f_float alpha, f_cfloat *A, f_int lda, f_float beta, f_cfloat *C, f_int ldc) {
    cherk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc, 1, 1);
}
void herk(char uplo, char trans, f_int n, f_int k, f_double alpha, f_cdouble *A, f_int lda, f_double beta, f_cdouble *C, f_int ldc) {
    zherk_(&uplo, &trans, &n, &k, &alpha, A, &lda, &beta, C, &ldc, 1, 1);
}

/// symmetric rank-2k update to a matrix
///     C := alpha * A * B.T + alpha.conj * B * A.T + beta * C
/// OR  C := alpha * A.T * B + alpha.conj * B.T * A + beta * C
void syr2k(char uplo, char trans, f_int n, f_int k, f_float alpha, f_float *A, f_int lda, f_float *B, f_int ldb, f_float beta, f_float *C, f_int ldc) {
    ssyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void syr2k(char uplo, char trans, f_int n, f_int k, f_double alpha, f_double *A, f_int lda, f_double *B, f_int ldb, f_double beta, f_double *C, f_int ldc) {
    dsyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void syr2k(char uplo, char trans, f_int n, f_int k, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *B, f_int ldb, f_cfloat beta, f_cfloat *C, f_int ldc) {
    csyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void syr2k(char uplo, char trans, f_int n, f_int k, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *B, f_int ldb, f_cdouble beta, f_cdouble *C, f_int ldc) {
    zsyr2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}

/// hermitian rank-2k update to a matrix
///     C := alpha * A * B.H + alpha.conj * B * A.H + beta * C
/// OR  C := alpha * A.H * B + alpha.conj * B.H * A + beta * C
void her2k(char uplo, char trans, f_int n, f_int k, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *B, f_int ldb, f_float beta, f_cfloat *C, f_int ldc) {
    cher2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}
void her2k(char uplo, char trans, f_int n, f_int k, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *B, f_int ldb, f_double beta, f_cdouble *C, f_int ldc) {
    zher2k_(&uplo, &trans, &n, &k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);
}

/// triangular matrix matrix multiply
///     B := alpha * transa(A) * B
/// OR  B := alpha * B * transa(A)
void trmm(char side, char uplo, char transa, char diag, f_int m, f_int n, f_float alpha, f_float *A, f_int lda, f_float *B, f_int ldb) {
    strmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
}
void trmm(char side, char uplo, char transa, char diag, f_int m, f_int n, f_double alpha, f_double *A, f_int lda, f_double *B, f_int ldb) {
    dtrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
}
void trmm(char side, char uplo, char transa, char diag, f_int m, f_int n, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *B, f_int ldb) {
    ctrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
}
void trmm(char side, char uplo, char transa, char diag, f_int m, f_int n, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *B, f_int ldb) {
    ztrmm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
}

/// solving triangular matrix with multiple right hand sides
///     B := alpha * transa(A.inv) * B
/// OR  B := alpha * B * transa(A.inv)
void trsm(char side, char uplo, char transa, char diag, f_int m, f_int n, f_float alpha, f_float *A, f_int lda, f_float *B, f_int ldb) {
    strsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
}
void trsm(char side, char uplo, char transa, char diag, f_int m, f_int n, f_double alpha, f_double *A, f_int lda, f_double *B, f_int ldb) {
    dtrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
}
void trsm(char side, char uplo, char transa, char diag, f_int m, f_int n, f_cfloat alpha, f_cfloat *A, f_int lda, f_cfloat *B, f_int ldb) {
    ctrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
}
void trsm(char side, char uplo, char transa, char diag, f_int m, f_int n, f_cdouble alpha, f_cdouble *A, f_int lda, f_cdouble *B, f_int ldb) {
    ztrsm_(&side, &uplo, &transa, &diag, &m, &n, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);
}

/// Test if the characters are equal. (Auxiliary routine in Level 2 and 3 BLAS routines)
// void lsame_() [no D interface]

/// Computes absolute values of a f_cdouble number. (Auxiliary routine for a few Level 1 BLAS routines)
// void dcabs1_() [no D interface]

/// Error handler for level 2 and 3 BLAS routines
// void xerbla_() [no D interface]

