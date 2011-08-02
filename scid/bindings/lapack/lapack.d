/** LAPACK bindings for D.

    Authors:    William V. Baxter III (with slight modifications by
                Lars Tandle Kyllingstad).
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.bindings.lapack.lapack;


public import scid.bindings.blas.types;


/*
  Copyright (C) 2006--2008 William V. Baxter III, OLM Digital, Inc.

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any
  damages arising from the use of this software.

  Permission is granted to anyone to use this software for any
  purpose, including commercial applications, and to alter it and
  redistribute it freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must
     not claim that you wrote the original software. If you use this
     software in a product, an acknowledgment in the product
     documentation would be appreciated but is not required.

  2. Altered source versions must be plainly marked as such, and must
     not be misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  William Baxter wbaxter@gmail.com
*/


// Prototypes for the raw Fortran interface to BLAS
extern(C):

alias f_int function(f_cfloat *) FCB_CGEES_SELECT;
alias f_int function(f_cfloat *) FCB_CGEESX_SELECT;
alias f_int function(f_cfloat *, f_cfloat *) FCB_CGGES_SELCTG;
alias f_int function(f_cfloat *, f_cfloat *) FCB_CGGESX_SELCTG;
alias f_int function(f_double *, f_double *) FCB_DGEES_SELECT;
alias f_int function(f_double *, f_double *) FCB_DGEESX_SELECT;
alias f_int function(f_double *, f_double *, f_double *) FCB_DGGES_DELCTG;
alias f_int function(f_double *, f_double *, f_double *) FCB_DGGESX_DELCTG;
alias f_int function(f_float *, f_float *) FCB_SGEES_SELECT;
alias f_int function(f_float *, f_float *) FCB_SGEESX_SELECT;
alias f_int function(f_float *, f_float *, f_float *) FCB_SGGES_SELCTG;
alias f_int function(f_float *, f_float *, f_float *) FCB_SGGESX_SELCTG;
alias f_int function(f_cdouble *) FCB_ZGEES_SELECT;
alias f_int function(f_cdouble *) FCB_ZGEESX_SELECT;
alias f_int function(f_cdouble *, f_cdouble *) FCB_ZGGES_DELCTG;
alias f_int function(f_cdouble *, f_cdouble *) FCB_ZGGESX_DELCTG;

version (FORTRAN_FLOAT_FUNCTIONS_RETURN_DOUBLE) {
    alias f_double lapack_float_ret_t;
} else {
    alias f_float lapack_float_ret_t;
}

/* LAPACK routines */

void dlaswp_( f_int *n, f_double  *a, f_int *lda, f_int *k1, f_int *k2, f_int *ipiv, f_int *incx );


//--------------------------------------------------------
// ---- SIMPLE and DIVIDE AND CONQUER DRIVER routines ----
//---------------------------------------------------------

/// Solves a general system of linear equations AX=B.
void sgesv_(f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_int *ipiv, f_float *b, f_int *ldb, f_int *info);
void dgesv_(f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_int *ipiv, f_double *b, f_int *ldb, f_int *info);
void cgesv_(f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info);
void zgesv_(f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info);

/// Solves a general banded system of linear equations AX=B.
void sgbsv_(f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_float *ab, f_int *ldab, f_int *ipiv, f_float *b, f_int *ldb, f_int *info);
void dgbsv_(f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_double *ab, f_int *ldab, f_int *ipiv, f_double *b, f_int *ldb, f_int *info);
void cgbsv_(f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info);
void zgbsv_(f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info);

/// Solves a general tridiagonal system of linear equations AX=B.
void sgtsv_(f_int *n, f_int *nrhs, f_float *dl, f_float *d, f_float *du, f_float *b, f_int *ldb, f_int *info);
void dgtsv_(f_int *n, f_int *nrhs, f_double *dl, f_double *d, f_double *du, f_double *b, f_int *ldb, f_int *info);
void cgtsv_(f_int *n, f_int *nrhs, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *b, f_int *ldb, f_int *info);
void zgtsv_(f_int *n, f_int *nrhs, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *b, f_int *ldb, f_int *info);

/// Solves a symmetric positive definite system of linear
/// equations AX=B.
void sposv_(char *uplo, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dposv_(char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void cposv_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zposv_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Solves a symmetric positive definite system of linear
/// equations AX=B, where A is held in packed storage.
void sppsv_(char *uplo, f_int *n, f_int *nrhs, f_float *ap, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dppsv_(char *uplo, f_int *n, f_int *nrhs, f_double *ap, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void cppsv_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zppsv_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Solves a symmetric positive definite banded system
/// of linear equations AX=B.
void spbsv_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_float *ab, f_int *ldab, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dpbsv_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_double *ab, f_int *ldab, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void cpbsv_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zpbsv_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Solves a symmetric positive definite tridiagonal system
/// of linear equations AX=B.
void sptsv_(f_int *n, f_int *nrhs, f_float *d, f_float *e, f_float *b, f_int *ldb, f_int *info);
void dptsv_(f_int *n, f_int *nrhs, f_double *d, f_double *e, f_double *b, f_int *ldb, f_int *info);
void cptsv_(f_int *n, f_int *nrhs, f_float *d, f_cfloat *e, f_cfloat *b, f_int *ldb, f_int *info);
void zptsv_(f_int *n, f_int *nrhs, f_double *d, f_cdouble *e, f_cdouble *b, f_int *ldb, f_int *info);


/// Solves a real symmetric indefinite system of linear equations AX=B.
void ssysv_(char *uplo, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_int *ipiv, f_float *b, f_int *ldb, f_float *work, f_int *lwork, f_int *info, f_int uplo_len);
void dsysv_(char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_int *ipiv, f_double *b, f_int *ldb, f_double *work, f_int *lwork, f_int *info, f_int uplo_len);
void csysv_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *work, f_int *lwork, f_int *info, f_int uplo_len);
void zsysv_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *work, f_int *lwork, f_int *info, f_int uplo_len);

/// Solves a complex Hermitian indefinite system of linear equations AX=B.
void chesv_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *work, f_int *lwork, f_int *info, f_int uplo_len);
void zhesv_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *work, f_int *lwork, f_int *info, f_int uplo_len);

/// Solves a real symmetric indefinite system of linear equations AX=B,
/// where A is held in packed storage.
void sspsv_(char *uplo, f_int *n, f_int *nrhs, f_float *ap, f_int *ipiv, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dspsv_(char *uplo, f_int *n, f_int *nrhs, f_double *ap, f_int *ipiv, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void cspsv_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zspsv_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Solves a complex Hermitian indefinite system of linear equations AX=B,
/// where A is held in packed storage.
void chpsv_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zhpsv_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Computes the least squares solution to an over-determined system
/// of linear equations, A X=B or A**H X=B,  or the minimum norm
/// solution of an under-determined system, where A is a general
/// rectangular matrix of full rank,  using a QR or LQ factorization
/// of A.
void sgels_(char *trans, f_int *m, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *work, f_int *lwork, f_int *info, f_int trans_len);
void dgels_(char *trans, f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *work, f_int *lwork, f_int *info, f_int trans_len);
void cgels_(char *trans, f_int *m, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *work, f_int *lwork, f_int *info, f_int trans_len);
void zgels_(char *trans, f_int *m, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *work, f_int *lwork, f_int *info, f_int trans_len);

/// Computes the least squares solution to an over-determined system
/// of linear equations, A X=B or A**H X=B,  or the minimum norm
/// solution of an under-determined system, using a divide and conquer
/// method, where A is a general rectangular matrix of full rank,
/// using a QR or LQ factorization of A.
void sgelsd_(f_int *m, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *s, f_float *rcond, f_int *rank, f_float *work, f_int *lwork, f_int *iwork, f_int *info);
void dgelsd_(f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *s, f_double *rcond, f_int *rank, f_double *work, f_int *lwork, f_int *iwork, f_int *info);
void cgelsd_(f_int *m, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_float *s, f_float *rcond, f_int *rank, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *iwork, f_int *info);
void zgelsd_(f_int *m, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_double *s, f_double *rcond, f_int *rank, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *iwork, f_int *info);

/// Solves the LSE (Constrained Linear Least Squares Problem) using
/// the GRQ (Generalized RQ) factorization
void sgglse_(f_int *m, f_int *n, f_int *p, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *c, f_float *d, f_float *x, f_float *work, f_int *lwork, f_int *info);
void dgglse_(f_int *m, f_int *n, f_int *p, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *c, f_double *d, f_double *x, f_double *work, f_int *lwork, f_int *info);
void cgglse_(f_int *m, f_int *n, f_int *p, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *c, f_cfloat *d, f_cfloat *x, f_cfloat *work, f_int *lwork, f_int *info);
void zgglse_(f_int *m, f_int *n, f_int *p, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *c, f_cdouble *d, f_cdouble *x, f_cdouble *work, f_int *lwork, f_int *info);

/// Solves the GLM (Generalized Linear Regression Model) using
/// the GQR (Generalized QR) factorization
void sggglm_(f_int *n, f_int *m, f_int *p, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *d, f_float *x, f_float *y, f_float *work, f_int *lwork, f_int *info);
void dggglm_(f_int *n, f_int *m, f_int *p, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *d, f_double *x, f_double *y, f_double *work, f_int *lwork, f_int *info);
void cggglm_(f_int *n, f_int *m, f_int *p, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *d, f_cfloat *x, f_cfloat *y, f_cfloat *work, f_int *lwork, f_int *info);
void zggglm_(f_int *n, f_int *m, f_int *p, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *d, f_cdouble *x, f_cdouble *y, f_cdouble *work, f_int *lwork, f_int *info);

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix.
void ssyev_(char *jobz, char *uplo, f_int *n, f_float *a, f_int *lda, f_float *w, f_float *work, f_int *lwork, f_int *info, f_int jobz_len, f_int uplo_len);
void dsyev_(char *jobz, char *uplo, f_int *n, f_double *a, f_int *lda, f_double *w, f_double *work, f_int *lwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian matrix.
void cheev_(char *jobz, char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_float *w, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zheev_(char *jobz, char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_double *w, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int jobz_len, f_int uplo_len);


/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix.  If eigenvectors are desired, it uses a divide
/// and conquer algorithm.
void ssyevd_(char *jobz, char *uplo, f_int *n, f_float *a, f_int *lda, f_float *w, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void dsyevd_(char *jobz, char *uplo, f_int *n, f_double *a, f_int *lda, f_double *w, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian matrix.  If eigenvectors are desired, it uses a divide
/// and conquer algorithm.
void cheevd_(char *jobz, char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_float *w, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zheevd_(char *jobz, char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_double *w, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix in packed storage.
void sspev_(char *jobz, char *uplo, f_int *n, f_float *ap, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *info, f_int jobz_len, f_int uplo_len);
void dspev_(char *jobz, char *uplo, f_int *n, f_double *ap, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes selected eigenvalues, and optionally, eigenvectors of a complex
/// Hermitian matrix.  Eigenvalues are computed by the dqds
/// algorithm, and eigenvectors are computed from various "good" LDL^T
/// representations (also known as Relatively Robust Representations).
/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian matrix in packed storage.
void chpev_(char *jobz, char *uplo, f_int *n, f_cfloat *ap, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_float *rwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhpev_(char *jobz, char *uplo, f_int *n, f_cdouble *ap, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_double *rwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix in packed storage.  If eigenvectors are desired,
/// it uses a divide and conquer algorithm.
void sspevd_(char *jobz, char *uplo, f_int *n, f_float *ap, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void dspevd_(char *jobz, char *uplo, f_int *n, f_double *ap, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian matrix in packed storage.  If eigenvectors are desired, it
/// uses a divide and conquer algorithm.
void chpevd_(char *jobz, char *uplo, f_int *n, f_cfloat *ap, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhpevd_(char *jobz, char *uplo, f_int *n, f_cdouble *ap, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric band matrix.
void ssbev_(char *jobz, char *uplo, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *info, f_int jobz_len, f_int uplo_len);
void dsbev_(char *jobz, char *uplo, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian band matrix.
void chbev_(char *jobz, char *uplo, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_float *rwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhbev_(char *jobz, char *uplo, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_double *rwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric band matrix.  If eigenvectors are desired, it uses a
/// divide and conquer algorithm.
void ssbevd_(char *jobz, char *uplo, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void dsbevd_(char *jobz, char *uplo, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian band matrix.  If eigenvectors are desired, it uses a divide
/// and conquer algorithm.
void chbevd_(char *jobz, char *uplo, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhbevd_(char *jobz, char *uplo, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric tridiagonal matrix.
void sstev_(char *jobz, f_int *n, f_float *d, f_float *e, f_float *z, f_int *ldz, f_float *work, f_int *info, f_int jobz_len);
void dstev_(char *jobz, f_int *n, f_double *d, f_double *e, f_double *z, f_int *ldz, f_double *work, f_int *info, f_int jobz_len);

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric tridiagonal matrix.  If eigenvectors are desired, it uses
/// a divide and conquer algorithm.
void sstevd_(char *jobz, f_int *n, f_float *d, f_float *e, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len);
void dstevd_(char *jobz, f_int *n, f_double *d, f_double *e, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len);

/// Computes the eigenvalues and Schur factorization of a general
/// matrix, and orders the factorization so that selected eigenvalues
/// are at the top left of the Schur form.
void sgees_(char *jobvs, char *sort, FCB_SGEES_SELECT select, f_int *n, f_float *a, f_int *lda, f_int *sdim, f_float *wr, f_float *wi, f_float *vs, f_int *ldvs, f_float *work, f_int *lwork, f_int *bwork, f_int *info, f_int jobvs_len, f_int sort_len);
void dgees_(char *jobvs, char *sort, FCB_DGEES_SELECT select, f_int *n, f_double *a, f_int *lda, f_int *sdim, f_double *wr, f_double *wi, f_double *vs, f_int *ldvs, f_double *work, f_int *lwork, f_int *bwork, f_int *info, f_int jobvs_len, f_int sort_len);
void cgees_(char *jobvs, char *sort, FCB_CGEES_SELECT select, f_int *n, f_cfloat *a, f_int *lda, f_int *sdim, f_cfloat *w, f_cfloat *vs, f_int *ldvs, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *bwork, f_int *info, f_int jobvs_len, f_int sort_len);
void zgees_(char *jobvs, char *sort, FCB_ZGEES_SELECT select, f_int *n, f_cdouble *a, f_int *lda, f_int *sdim, f_cdouble *w, f_cdouble *vs, f_int *ldvs, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *bwork, f_int *info, f_int jobvs_len, f_int sort_len);

/// Computes the eigenvalues and left and right eigenvectors of
/// a general matrix.
void sgeev_(char *jobvl, char *jobvr, f_int *n, f_float *a, f_int *lda, f_float *wr, f_float *wi, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_float *work, f_int *lwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void dgeev_(char *jobvl, char *jobvr, f_int *n, f_double *a, f_int *lda, f_double *wr, f_double *wi, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_double *work, f_int *lwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void cgeev_(char *jobvl, char *jobvr, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *w, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void zgeev_(char *jobvl, char *jobvr, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *w, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int jobvl_len, f_int jobvr_len);

/// Computes the singular value decomposition (SVD) of a general
/// rectangular matrix.
void sgesvd_(char *jobu, char *jobvt, f_int *m, f_int *n, f_float *a, f_int *lda, f_float *s, f_float *u, f_int *ldu, f_float *vt, f_int *ldvt, f_float *work, f_int *lwork, f_int *info, f_int jobu_len, f_int jobvt_len);
void dgesvd_(char *jobu, char *jobvt, f_int *m, f_int *n, f_double *a, f_int *lda, f_double *s, f_double *u, f_int *ldu, f_double *vt, f_int *ldvt, f_double *work, f_int *lwork, f_int *info, f_int jobu_len, f_int jobvt_len);
void cgesvd_(char *jobu, char *jobvt, f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_float *s, f_cfloat *u, f_int *ldu, f_cfloat *vt, f_int *ldvt, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int jobu_len, f_int jobvt_len);
void zgesvd_(char *jobu, char *jobvt, f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_double *s, f_cdouble *u, f_int *ldu, f_cdouble *vt, f_int *ldvt, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int jobu_len, f_int jobvt_len);

/// Computes the singular value decomposition (SVD) of a general
/// rectangular matrix using divide-and-conquer.
void sgesdd_(char *jobz, f_int *m, f_int *n, f_float *a, f_int *lda, f_float *s, f_float *u, f_int *ldu, f_float *vt, f_int *ldvt, f_float *work, f_int *lwork, f_int *iwork, f_int *info, f_int jobz_len);
void dgesdd_(char *jobz, f_int *m, f_int *n, f_double *a, f_int *lda, f_double *s, f_double *u, f_int *ldu, f_double *vt, f_int *ldvt, f_double *work, f_int *lwork, f_int *iwork, f_int *info, f_int jobz_len);
void cgesdd_(char *jobz, f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_float *s, f_cfloat *u, f_int *ldu, f_cfloat *vt, f_int *ldvt, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *iwork, f_int *info, f_int jobz_len);
void zgesdd_(char *jobz, f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_double *s, f_cdouble *u, f_int *ldu, f_cdouble *vt, f_int *ldvt, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *iwork, f_int *info, f_int jobz_len);

/// Computes all eigenvalues and the eigenvectors of  a generalized
/// symmetric-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
void ssygv_(f_int *itype, char *jobz, char *uplo, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *w, f_float *work, f_int *lwork, f_int *info, f_int jobz_len, f_int uplo_len);
void dsygv_(f_int *itype, char *jobz, char *uplo, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *w, f_double *work, f_int *lwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and the eigenvectors of  a generalized
/// Hermitian-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
void chegv_(f_int *itype, char *jobz, char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_float *w, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhegv_(f_int *itype, char *jobz, char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_double *w, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and the eigenvectors of  a generalized
/// symmetric-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void ssygvd_(f_int *itype, char *jobz, char *uplo, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *w, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void dsygvd_(f_int *itype, char *jobz, char *uplo, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *w, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
/// Computes all eigenvalues and the eigenvectors of  a generalized
/// Hermitian-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void chegvd_(f_int *itype, char *jobz, char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_float *w, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhegvd_(f_int *itype, char *jobz, char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_double *w, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and eigenvectors of  a generalized
/// symmetric-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
void sspgv_(f_int *itype, char *jobz, char *uplo, f_int *n, f_float *ap, f_float *bp, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *info, f_int jobz_len, f_int uplo_len);
void dspgv_(f_int *itype, char *jobz, char *uplo, f_int *n, f_double *ap, f_double *bp, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and eigenvectors of  a generalized
/// Hermitian-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
void chpgv_(f_int *itype, char *jobz, char *uplo, f_int *n, f_cfloat *ap, f_cfloat *bp, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_float *rwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhpgv_(f_int *itype, char *jobz, char *uplo, f_int *n, f_cdouble *ap, f_cdouble *bp, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_double *rwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and eigenvectors of  a generalized
/// symmetric-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void sspgvd_(f_int *itype, char *jobz, char *uplo, f_int *n, f_float *ap, f_float *bp, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void dspgvd_(f_int *itype, char *jobz, char *uplo, f_int *n, f_double *ap, f_double *bp, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all eigenvalues and eigenvectors of  a generalized
/// Hermitian-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void chpgvd_(f_int *itype, char *jobz, char *uplo, f_int *n, f_cfloat *ap, f_cfloat *bp, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhpgvd_(f_int *itype, char *jobz, char *uplo, f_int *n, f_cdouble *ap, f_cdouble *bp, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all the eigenvalues, and optionally, the eigenvectors
/// of a real generalized symmetric-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be symmetric
/// and banded, and B is also positive definite.
void ssbgv_(char *jobz, char *uplo, f_int *n, f_int *ka, f_int *kb, f_float *ab, f_int *ldab, f_float *bb, f_int *ldbb, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *info, f_int jobz_len, f_int uplo_len);
void dsbgv_(char *jobz, char *uplo, f_int *n, f_int *ka, f_int *kb, f_double *ab, f_int *ldab, f_double *bb, f_int *ldbb, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all the eigenvalues, and optionally, the eigenvectors
/// of a complex generalized Hermitian-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be Hermitian
/// and banded, and B is also positive definite.
void chbgv_(char *jobz, char *uplo, f_int *n, f_int *ka, f_int *kb, f_cfloat *ab, f_int *ldab, f_cfloat *bb, f_int *ldbb, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_float *rwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhbgv_(char *jobz, char *uplo, f_int *n, f_int *ka, f_int *kb, f_cdouble *ab, f_int *ldab, f_cdouble *bb, f_int *ldbb, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_double *rwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all the eigenvalues, and optionally, the eigenvectors
/// of a real generalized symmetric-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be symmetric
/// and banded, and B is also positive definite.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void ssbgvd_(char *jobz, char *uplo, f_int *n, f_int *ka, f_int *kb, f_float *ab, f_int *ldab, f_float *bb, f_int *ldbb, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void dsbgvd_(char *jobz, char *uplo, f_int *n, f_int *ka, f_int *kb, f_double *ab, f_int *ldab, f_double *bb, f_int *ldbb, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes all the eigenvalues, and optionally, the eigenvectors
/// of a complex generalized Hermitian-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be Hermitian
/// and banded, and B is also positive definite.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void chbgvd_(char *jobz, char *uplo, f_int *n, f_int *ka, f_int *kb, f_cfloat *ab, f_int *ldab, f_cfloat *bb, f_int *ldbb, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);
void zhbgvd_(char *jobz, char *uplo, f_int *n, f_int *ka, f_int *kb, f_cdouble *ab, f_int *ldab, f_cdouble *bb, f_int *ldbb, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int uplo_len);

/// Computes the generalized eigenvalues, Schur form, and left and/or
/// right Schur vectors for a pair of nonsymmetric matrices
void sgegs_(char *jobvsl, char *jobvsr, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *vsl, f_int *ldvsl, f_float *vsr, f_int *ldvsr, f_float *work, f_int *lwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len);
void dgegs_(char *jobvsl, char *jobvsr, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *vsl, f_int *ldvsl, f_double *vsr, f_int *ldvsr, f_double *work, f_int *lwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len);
void cgegs_(char *jobvsl, char *jobvsr, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vsl, f_int *ldvsl, f_cfloat *vsr, f_int *ldvsr, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len);
void zgegs_(char *jobvsl, char *jobvsr, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vsl, f_int *ldvsl, f_cdouble *vsr, f_int *ldvsr, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len);

/// Computes the generalized eigenvalues, Schur form, and left and/or
/// right Schur vectors for a pair of nonsymmetric matrices
void sgges_(char *jobvsl, char *jobvsr, char *sort, FCB_SGGES_SELCTG selctg, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *sdim, f_float *alphar, f_float *alphai, f_float *betav, f_float *vsl, f_int *ldvsl, f_float *vsr, f_int *ldvsr, f_float *work, f_int *lwork, f_int *bwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len, f_int sort_len);
void dgges_(char *jobvsl, char *jobvsr, char *sort, FCB_DGGES_DELCTG delctg, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *sdim, f_double *alphar, f_double *alphai, f_double *betav, f_double *vsl, f_int *ldvsl, f_double *vsr, f_int *ldvsr, f_double *work, f_int *lwork, f_int *bwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len, f_int sort_len);
void cgges_(char *jobvsl, char *jobvsr, char *sort, FCB_CGGES_SELCTG selctg, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *sdim, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vsl, f_int *ldvsl, f_cfloat *vsr, f_int *ldvsr, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *bwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len, f_int sort_len);
void zgges_(char *jobvsl, char *jobvsr, char *sort, FCB_ZGGES_DELCTG delctg, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *sdim, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vsl, f_int *ldvsl, f_cdouble *vsr, f_int *ldvsr, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *bwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len, f_int sort_len);

/// Computes the generalized eigenvalues, and left and/or right
/// generalized eigenvectors for a pair of nonsymmetric matrices
void sgegv_(char *jobvl, char *jobvr, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_float *work, f_int *lwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void dgegv_(char *jobvl, char *jobvr, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_double *work, f_int *lwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void cgegv_(char *jobvl, char *jobvr, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *alphar, f_cfloat *betav, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void zgegv_(char *jobvl, char *jobvr, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *alphar, f_cdouble *betav, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int jobvl_len, f_int jobvr_len);

/// Computes the generalized eigenvalues, and left and/or right
/// generalized eigenvectors for a pair of nonsymmetric matrices
void sggev_(char *jobvl, char *jobvr, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_float *work, f_int *lwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void dggev_(char *jobvl, char *jobvr, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_double *work, f_int *lwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void cggev_(char *jobvl, char *jobvr, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int jobvl_len, f_int jobvr_len);
void zggev_(char *jobvl, char *jobvr, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int jobvl_len, f_int jobvr_len);

/// Computes the Generalized Singular Value Decomposition
void sggsvd_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *n, f_int *p, f_int *k, f_int *l, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *alphav, f_float *betav, f_float *u, f_int *ldu, f_float *v, f_int *ldv, f_float *q, f_int *ldq, f_float *work, f_int *iwork, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void dggsvd_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *n, f_int *p, f_int *k, f_int *l, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *alphav, f_double *betav, f_double *u, f_int *ldu, f_double *v, f_int *ldv, f_double *q, f_int *ldq, f_double *work, f_int *iwork, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void cggsvd_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *n, f_int *p, f_int *k, f_int *l, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_float *alphav, f_float *betav, f_cfloat *u, f_int *ldu, f_cfloat *v, f_int *ldv, f_cfloat *q, f_int *ldq, f_cfloat *work, f_float *rwork, f_int *iwork, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void zggsvd_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *n, f_int *p, f_int *k, f_int *l, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_double *alphav, f_double *betav, f_cdouble *u, f_int *ldu, f_cdouble *v, f_int *ldv, f_cdouble *q, f_int *ldq, f_cdouble *work, f_double *rwork, f_int *iwork, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);

//-----------------------------------------------------
//       ---- EXPERT and RRR DRIVER routines ----
//-----------------------------------------------------

/// Solves a general system of linear equations AX=B, A**T X=B
/// or A**H X=B, and provides an estimate of the condition number
/// and error bounds on the solution.
void sgesvx_(char *fact, char *trans, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *af, f_int *ldaf, f_int *ipiv, char *equed, f_float *r, f_float *c, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int fact_len, f_int trans_len, f_int equed_len);
void dgesvx_(char *fact, char *trans, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *af, f_int *ldaf, f_int *ipiv, char *equed, f_double *r, f_double *c, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int fact_len, f_int trans_len, f_int equed_len);
void cgesvx_(char *fact, char *trans, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *af, f_int *ldaf, f_int *ipiv, char *equed, f_float *r, f_float *c, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len, f_int trans_len, f_int equed_len);
void zgesvx_(char *fact, char *trans, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *af, f_int *ldaf, f_int *ipiv, char *equed, f_double *r, f_double *c, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len, f_int trans_len, f_int equed_len);

/// Solves a general banded system of linear equations AX=B,
/// A**T X=B or A**H X=B, and provides an estimate of the condition
/// number and error bounds on the solution.
void sgbsvx_(char *fact, char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_float *ab, f_int *ldab, f_float *afb, f_int *ldafb, f_int *ipiv, char *equed, f_float *r, f_float *c, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int fact_len, f_int trans_len, f_int equed_len);
void dgbsvx_(char *fact, char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_double *ab, f_int *ldab, f_double *afb, f_int *ldafb, f_int *ipiv, char *equed, f_double *r, f_double *c, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int fact_len, f_int trans_len, f_int equed_len);
void cgbsvx_(char *fact, char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_cfloat *afb, f_int *ldafb, f_int *ipiv, char *equed, f_float *r, f_float *c, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len, f_int trans_len, f_int equed_len);
void zgbsvx_(char *fact, char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_cdouble *afb, f_int *ldafb, f_int *ipiv, char *equed, f_double *r, f_double *c, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len, f_int trans_len, f_int equed_len);

/// Solves a general tridiagonal system of linear equations AX=B,
/// A**T X=B or A**H X=B, and provides an estimate of the condition
/// number  and error bounds on the solution.
void sgtsvx_(char *fact, char *trans, f_int *n, f_int *nrhs, f_float *dl, f_float *d, f_float *du, f_float *dlf, f_float *df, f_float *duf, f_float *du2, f_int *ipiv, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int fact_len, f_int trans_len);
void dgtsvx_(char *fact, char *trans, f_int *n, f_int *nrhs, f_double *dl, f_double *d, f_double *du, f_double *dlf, f_double *df, f_double *duf, f_double *du2, f_int *ipiv, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int fact_len, f_int trans_len);
void cgtsvx_(char *fact, char *trans, f_int *n, f_int *nrhs, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *dlf, f_cfloat *df, f_cfloat *duf, f_cfloat *du2, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len, f_int trans_len);
void zgtsvx_(char *fact, char *trans, f_int *n, f_int *nrhs, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *dlf, f_cdouble *df, f_cdouble *duf, f_cdouble *du2, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len, f_int trans_len);

/// Solves a symmetric positive definite system of linear
/// equations AX=B, and provides an estimate of the condition number
/// and error bounds on the solution.
void sposvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *af, f_int *ldaf, char *equed, f_float *s, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void dposvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *af, f_int *ldaf, char *equed, f_double *s, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void cposvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *af, f_int *ldaf, char *equed, f_float *s, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void zposvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *af, f_int *ldaf, char *equed, f_double *s, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);

/// Solves a symmetric positive definite system of linear
/// equations AX=B, where A is held in packed storage, and provides
/// an estimate of the condition number and error bounds on the
/// solution.
void sppsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_float *ap, f_float *afp, char *equed, f_float *s, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void dppsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_double *ap, f_double *afp, char *equed, f_double *s, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void cppsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *afp, char *equed, f_float *s, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void zppsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *afp, char *equed, f_double *s, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);

/// Solves a symmetric positive definite banded system
/// of linear equations AX=B, and provides an estimate of the condition
/// number and error bounds on the solution.
void spbsvx_(char *fact, char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_float *ab, f_int *ldab, f_float *afb, f_int *ldafb, char *equed, f_float *s, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void dpbsvx_(char *fact, char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_double *ab, f_int *ldab, f_double *afb, f_int *ldafb, char *equed, f_double *s, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void cpbsvx_(char *fact, char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_cfloat *afb, f_int *ldafb, char *equed, f_float *s, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);
void zpbsvx_(char *fact, char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_cdouble *afb, f_int *ldafb, char *equed, f_double *s, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len, f_int uplo_len, f_int equed_len);

/// Solves a symmetric positive definite tridiagonal
/// system of linear equations AX=B, and provides an estimate of
/// the condition number and error bounds on the solution.
void sptsvx_(char *fact, f_int *n, f_int *nrhs, f_float *d, f_float *e, f_float *df, f_float *ef, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *info, f_int fact_len);
void dptsvx_(char *fact, f_int *n, f_int *nrhs, f_double *d, f_double *e, f_double *df, f_double *ef, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *info, f_int fact_len);
void cptsvx_(char *fact, f_int *n, f_int *nrhs, f_float *d, f_cfloat *e, f_float *df, f_cfloat *ef, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len);
void zptsvx_(char *fact, f_int *n, f_int *nrhs, f_double *d, f_cdouble *e, f_double *df, f_cdouble *ef, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len);

/// Solves a real symmetric
/// indefinite system  of linear equations AX=B, and provides an
/// estimate of the condition number and error bounds on the solution.
void ssysvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *af, f_int *ldaf, f_int *ipiv, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *lwork, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len);
void dsysvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *af, f_int *ldaf, f_int *ipiv, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *lwork, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len);
void csysvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *af, f_int *ldaf, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int fact_len, f_int uplo_len);
void zsysvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *af, f_int *ldaf, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int fact_len, f_int uplo_len);

/// Solves a complex Hermitian
/// indefinite system  of linear equations AX=B, and provides an
/// estimate of the condition number and error bounds on the solution.
void chesvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *af, f_int *ldaf, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int fact_len, f_int uplo_len);
void zhesvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *af, f_int *ldaf, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int fact_len, f_int uplo_len);

/// Solves a real symmetric
/// indefinite system of linear equations AX=B, where A is held
/// in packed storage, and provides an estimate of the condition
/// number and error bounds on the solution.
void sspsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_float *ap, f_float *afp, f_int *ipiv, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len);
void dspsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_double *ap, f_double *afp, f_int *ipiv, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int fact_len, f_int uplo_len);
void cspsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *afp, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len, f_int uplo_len);
void zspsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *afp, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len, f_int uplo_len);

/// Solves a complex Hermitian
/// indefinite system of linear equations AX=B, where A is held
/// in packed storage, and provides an estimate of the condition
/// number and error bounds on the solution.
void chpsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *afp, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int fact_len, f_int uplo_len);
void zhpsvx_(char *fact, char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *afp, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int fact_len, f_int uplo_len);

/// Computes the minimum norm least squares solution to an over-
/// or under-determined system of linear equations A X=B, using a
/// complete orthogonal factorization of A.
void sgelsx_(f_int *m, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *jpvt, f_float *rcond, f_int *rank, f_float *work, f_int *info);
void dgelsx_(f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *jpvt, f_double *rcond, f_int *rank, f_double *work, f_int *info);
void cgelsx_(f_int *m, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *jpvt, f_float *rcond, f_int *rank, f_cfloat *work, f_float *rwork, f_int *info);
void zgelsx_(f_int *m, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *jpvt, f_double *rcond, f_int *rank, f_cdouble *work, f_double *rwork, f_int *info);

/// Computes the minimum norm least squares solution to an over-
/// or under-determined system of linear equations A X=B, using a
/// complete orthogonal factorization of A.
void sgelsy_(f_int *m, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *jpvt, f_float *rcond, f_int *rank, f_float *work, f_int *lwork, f_int *info);
void dgelsy_(f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *jpvt, f_double *rcond, f_int *rank, f_double *work, f_int *lwork, f_int *info);
void cgelsy_(f_int *m, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *jpvt, f_float *rcond, f_int *rank, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info);
void zgelsy_(f_int *m, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *jpvt, f_double *rcond, f_int *rank, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info);

/// Computes the minimum norm least squares solution to an over-
/// or under-determined system of linear equations A X=B,  using
/// the singular value decomposition of A.
void sgelss_(f_int *m, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *s, f_float *rcond, f_int *rank, f_float *work, f_int *lwork, f_int *info);
void dgelss_(f_int *m, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *s, f_double *rcond, f_int *rank, f_double *work, f_int *lwork, f_int *info);
void cgelss_(f_int *m, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_float *s, f_float *rcond, f_int *rank, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info);
void zgelss_(f_int *m, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_double *s, f_double *rcond, f_int *rank, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info);

/// Computes selected eigenvalues and eigenvectors of a symmetric matrix.
void ssyevx_(char *jobz, char *range, char *uplo, f_int *n, f_float *a, f_int *lda, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void dsyevx_(char *jobz, char *range, char *uplo, f_int *n, f_double *a, f_int *lda, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues and eigenvectors of a Hermitian matrix.
void cheevx_(char *jobz, char *range, char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void zheevx_(char *jobz, char *range, char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix.  Eigenvalues are computed by the dqds
/// algorithm, and eigenvectors are computed from various "good" LDL^T
/// representations (also known as Relatively Robust Representations).
void ssyevr_(char *jobz, char *range, char *uplo, f_int *n, f_float *a, f_int *lda, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_int *isuppz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void dsyevr_(char *jobz, char *range, char *uplo, f_int *n, f_double *a, f_int *lda, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_int *isuppz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues, and optionally, eigenvectors of a complex
/// Hermitian matrix.  Eigenvalues are computed by the dqds
/// algorithm, and eigenvectors are computed from various "good" LDL^T
/// representations (also known as Relatively Robust Representations).
void cheevr_(char *jobz, char *range, char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_cfloat *z, f_int *ldz, f_int *isuppz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void zheevr_(char *jobz, char *range, char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_cdouble *z, f_int *ldz, f_int *isuppz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);


/// Computes selected eigenvalues, and optionally, the eigenvectors of
/// a generalized symmetric-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
void ssygvx_(f_int *itype, char *jobz, char *range, char *uplo, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void dsygvx_(f_int *itype, char *jobz, char *range, char *uplo, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues, and optionally, the eigenvectors of
/// a generalized Hermitian-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
void chegvx_(f_int *itype, char *jobz, char *range, char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void zhegvx_(f_int *itype, char *jobz, char *range, char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues and eigenvectors of a
/// symmetric matrix in packed storage.
void sspevx_(char *jobz, char *range, char *uplo, f_int *n, f_float *ap, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void dspevx_(char *jobz, char *range, char *uplo, f_int *n, f_double *ap, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues and eigenvectors of a
/// Hermitian matrix in packed storage.
void chpevx_(char *jobz, char *range, char *uplo, f_int *n, f_cfloat *ap, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_float *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void zhpevx_(char *jobz, char *range, char *uplo, f_int *n, f_cdouble *ap, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_double *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues, and optionally, eigenvectors of
/// a generalized symmetric-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
void sspgvx_(f_int *itype, char *jobz, char *range, char *uplo, f_int *n, f_float *ap, f_float *bp, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void dspgvx_(f_int *itype, char *jobz, char *range, char *uplo, f_int *n, f_double *ap, f_double *bp, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues, and optionally, the eigenvectors of
/// a generalized Hermitian-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
void chpgvx_(f_int *itype, char *jobz, char *range, char *uplo, f_int *n, f_cfloat *ap, f_cfloat *bp, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_float *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void zhpgvx_(f_int *itype, char *jobz, char *range, char *uplo, f_int *n, f_cdouble *ap, f_cdouble *bp, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_double *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues and eigenvectors of a
/// symmetric band matrix.
void ssbevx_(char *jobz, char *range, char *uplo, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_float *q, f_int *ldq, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void dsbevx_(char *jobz, char *range, char *uplo, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_double *q, f_int *ldq, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues and eigenvectors of a
/// Hermitian band matrix.
void chbevx_(char *jobz, char *range, char *uplo, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_cfloat *q, f_int *ldq, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_float *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void zhbevx_(char *jobz, char *range, char *uplo, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_cdouble *q, f_int *ldq, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_double *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues, and optionally, the eigenvectors
/// of a real generalized symmetric-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be symmetric
/// and banded, and B is also positive definite.
void ssbgvx_(char *jobz, char *range, char *uplo, f_int *n, f_int *ka, f_int *kb, f_float *ab, f_int *ldab, f_float *bb, f_int *ldbb, f_float *q, f_int *ldq, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void dsbgvx_(char *jobz, char *range, char *uplo, f_int *n, f_int *ka, f_int *kb, f_double *ab, f_int *ldab, f_double *bb, f_int *ldbb, f_double *q, f_int *ldq, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues, and optionally, the eigenvectors
/// of a complex generalized Hermitian-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be Hermitian
/// and banded, and B is also positive definite.
void chbgvx_(char *jobz, char *range, char *uplo, f_int *n, f_int *ka, f_int *kb, f_cfloat *ab, f_int *ldab, f_cfloat *bb, f_int *ldbb, f_cfloat *q, f_int *ldq, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_float *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);
void zhbgvx_(char *jobz, char *range, char *uplo, f_int *n, f_int *ka, f_int *kb, f_cdouble *ab, f_int *ldab, f_cdouble *bb, f_int *ldbb, f_cdouble *q, f_int *ldq, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_double *rwork, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len, f_int uplo_len);

/// Computes selected eigenvalues and eigenvectors of a real
/// symmetric tridiagonal matrix.
void sstevx_(char *jobz, char *range, f_int *n, f_float *d, f_float *e, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_float *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len);
void dstevx_(char *jobz, char *range, f_int *n, f_double *d, f_double *e, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_double *work, f_int *iwork, f_int *ifail, f_int *info, f_int jobz_len, f_int range_len);

/// Computes selected eigenvalues, and optionally, eigenvectors of a real
/// symmetric tridiagonal matrix.  Eigenvalues are computed by the dqds
/// algorithm, and eigenvectors are computed from various "good" LDL^T
/// representations (also known as Relatively Robust Representations).
void sstevr_(char *jobz, char *range, f_int *n, f_float *d, f_float *e, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_int *isuppz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len);
void dstevr_(char *jobz, char *range, f_int *n, f_double *d, f_double *e, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_int *isuppz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len);

/// Computes the eigenvalues and Schur factorization of a general
/// matrix, orders the factorization so that selected eigenvalues
/// are at the top left of the Schur form, and computes reciprocal
/// condition numbers for the average of the selected eigenvalues,
/// and for the associated right invariant subspace.
void sgeesx_(char *jobvs, char *sort, FCB_SGEESX_SELECT select, char *sense, f_int *n, f_float *a, f_int *lda, f_int *sdim, f_float *wr, f_float *wi, f_float *vs, f_int *ldvs, f_float *rconde, f_float *rcondv, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *bwork, f_int *info, f_int jobvs_len, f_int sort_len, f_int sense_len);
void dgeesx_(char *jobvs, char *sort, FCB_DGEESX_SELECT select, char *sense, f_int *n, f_double *a, f_int *lda, f_int *sdim, f_double *wr, f_double *wi, f_double *vs, f_int *ldvs, f_double *rconde, f_double *rcondv, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *bwork, f_int *info, f_int jobvs_len, f_int sort_len, f_int sense_len);
void cgeesx_(char *jobvs, char *sort, FCB_CGEESX_SELECT select, char *sense, f_int *n, f_cfloat *a, f_int *lda, f_int *sdim, f_cfloat *w, f_cfloat *vs, f_int *ldvs, f_float *rconde, f_float *rcondv, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *bwork, f_int *info, f_int jobvs_len, f_int sort_len, f_int sense_len);
void zgeesx_(char *jobvs, char *sort, FCB_ZGEESX_SELECT select, char *sense, f_int *n, f_cdouble *a, f_int *lda, f_int *sdim, f_cdouble *w, f_cdouble *vs, f_int *ldvs, f_double *rconde, f_double *rcondv, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *bwork, f_int *info, f_int jobvs_len, f_int sort_len, f_int sense_len);

/// Computes the generalized eigenvalues, the real Schur form, and,
/// optionally, the left and/or right matrices of Schur vectors.
void sggesx_(char *jobvsl, char *jobvsr, char *sort, FCB_SGGESX_SELCTG selctg, char *sense, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *sdim, f_float *alphar, f_float *alphai, f_float *betav, f_float *vsl, f_int *ldvsl, f_float *vsr, f_int *ldvsr, f_float *rconde, f_float *rcondv, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *bwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len, f_int sort_len, f_int sense_len);
void dggesx_(char *jobvsl, char *jobvsr, char *sort, FCB_DGGESX_DELCTG delctg, char *sense, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *sdim, f_double *alphar, f_double *alphai, f_double *betav, f_double *vsl, f_int *ldvsl, f_double *vsr, f_int *ldvsr, f_double *rconde, f_double *rcondv, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *bwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len, f_int sort_len, f_int sense_len);
void cggesx_(char *jobvsl, char *jobvsr, char *sort, FCB_CGGESX_SELCTG selctg, char *sense, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *sdim, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vsl, f_int *ldvsl, f_cfloat *vsr, f_int *ldvsr, f_float *rconde, f_float *rcondv, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *iwork, f_int *liwork, f_int *bwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len, f_int sort_len, f_int sense_len);
void zggesx_(char *jobvsl, char *jobvsr, char *sort, FCB_ZGGESX_DELCTG delctg, char *sense, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *sdim, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vsl, f_int *ldvsl, f_cdouble *vsr, f_int *ldvsr, f_double *rconde, f_double *rcondv, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *iwork, f_int *liwork, f_int *bwork, f_int *info, f_int jobvsl_len, f_int jobvsr_len, f_int sort_len, f_int sense_len);

/// Computes the eigenvalues and left and right eigenvectors of
/// a general matrix,  with preliminary balancing of the matrix,
/// and computes reciprocal condition numbers for the eigenvalues
/// and right eigenvectors.
void sgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, f_int *n, f_float *a, f_int *lda, f_float *wr, f_float *wi, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_int *ilo, f_int *ihi, f_float *scale, f_float *abnrm, f_float *rconde, f_float *rcondv, f_float *work, f_int *lwork, f_int *iwork, f_int *info, f_int balanc_len, f_int jobvl_len, f_int jobvr_len, f_int sense_len);
void dgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, f_int *n, f_double *a, f_int *lda, f_double *wr, f_double *wi, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_int *ilo, f_int *ihi, f_double *scale, f_double *abnrm, f_double *rconde, f_double *rcondv, f_double *work, f_int *lwork, f_int *iwork, f_int *info, f_int balanc_len, f_int jobvl_len, f_int jobvr_len, f_int sense_len);
void cgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *w, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_int *ilo, f_int *ihi, f_float *scale, f_float *abnrm, f_float *rconde, f_float *rcondv, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int balanc_len, f_int jobvl_len, f_int jobvr_len, f_int sense_len);
void zgeevx_(char *balanc, char *jobvl, char *jobvr, char *sense, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *w, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_int *ilo, f_int *ihi, f_double *scale, f_double *abnrm, f_double *rconde, f_double *rcondv, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int balanc_len, f_int jobvl_len, f_int jobvr_len, f_int sense_len);

/// Computes the generalized eigenvalues, and optionally, the left
/// and/or right generalized eigenvectors.
void sggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_int *ilo, f_int *ihi, f_float *lscale, f_float *rscale, f_float *abnrm, f_float *bbnrm, f_float *rconde, f_float *rcondv, f_float *work, f_int *lwork, f_int *iwork, f_int *bwork, f_int *info, f_int balanc_len, f_int jobvl_len, f_int jobvr_len, f_int sense_len);
void dggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_int *ilo, f_int *ihi, f_double *lscale, f_double *rscale, f_double *abnrm, f_double *bbnrm, f_double *rconde, f_double *rcondv, f_double *work, f_int *lwork, f_int *iwork, f_int *bwork, f_int *info, f_int balanc_len, f_int jobvl_len, f_int jobvr_len, f_int sense_len);
void cggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_int *ilo, f_int *ihi, f_float *lscale, f_float *rscale, f_float *abnrm, f_float *bbnrm, f_float *rconde, f_float *rcondv, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *iwork, f_int *bwork, f_int *info, f_int balanc_len, f_int jobvl_len, f_int jobvr_len, f_int sense_len);
void zggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_int *ilo, f_int *ihi, f_double *lscale, f_double *rscale, f_double *abnrm, f_double *bbnrm, f_double *rconde, f_double *rcondv, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *iwork, f_int *bwork, f_int *info, f_int balanc_len, f_int jobvl_len, f_int jobvr_len, f_int sense_len);



//----------------------------------------
//    ---- COMPUTATIONAL routines ----
//----------------------------------------


/// Computes the singular value decomposition (SVD) of a real bidiagonal
/// matrix, using a divide and conquer method.
void sbdsdc_(char *uplo, char *compq, f_int *n, f_float *d, f_float *e, f_float *u, f_int *ldu, f_float *vt, f_int *ldvt, f_float *q, f_int *iq, f_float *work, f_int *iwork, f_int *info, f_int uplo_len, f_int compq_len);
void dbdsdc_(char *uplo, char *compq, f_int *n, f_double *d, f_double *e, f_double *u, f_int *ldu, f_double *vt, f_int *ldvt, f_double *q, f_int *iq, f_double *work, f_int *iwork, f_int *info, f_int uplo_len, f_int compq_len);

/// Computes the singular value decomposition (SVD) of a real bidiagonal
/// matrix, using the bidiagonal QR algorithm.
void sbdsqr_(char *uplo, f_int *n, f_int *ncvt, f_int *nru, f_int *ncc, f_float *d, f_float *e, f_float *vt, f_int *ldvt, f_float *u, f_int *ldu, f_float *c, f_int *ldc, f_float *work, f_int *info, f_int uplo_len);
void dbdsqr_(char *uplo, f_int *n, f_int *ncvt, f_int *nru, f_int *ncc, f_double *d, f_double *e, f_double *vt, f_int *ldvt, f_double *u, f_int *ldu, f_double *c, f_int *ldc, f_double *work, f_int *info, f_int uplo_len);
void cbdsqr_(char *uplo, f_int *n, f_int *ncvt, f_int *nru, f_int *ncc, f_float *d, f_float *e, f_cfloat *vt, f_int *ldvt, f_cfloat *u, f_int *ldu, f_cfloat *c, f_int *ldc, f_float *rwork, f_int *info, f_int uplo_len);
void zbdsqr_(char *uplo, f_int *n, f_int *ncvt, f_int *nru, f_int *ncc, f_double *d, f_double *e, f_cdouble *vt, f_int *ldvt, f_cdouble *u, f_int *ldu, f_cdouble *c, f_int *ldc, f_double *rwork, f_int *info, f_int uplo_len);

/// Computes the reciprocal condition numbers for the eigenvectors of a
/// real symmetric or complex Hermitian matrix or for the left or right
/// singular vectors of a general matrix.
void sdisna_(char *job, f_int *m, f_int *n, f_float *d, f_float *sep, f_int *info, f_int job_len);
void ddisna_(char *job, f_int *m, f_int *n, f_double *d, f_double *sep, f_int *info, f_int job_len);

/// Reduces a general band matrix to real upper bidiagonal form
/// by an orthogonal transformation.
void sgbbrd_(char *vect, f_int *m, f_int *n, f_int *ncc, f_int *kl, f_int *ku, f_float *ab, f_int *ldab, f_float *d, f_float *e, f_float *q, f_int *ldq, f_float *pt, f_int *ldpt, f_float *c, f_int *ldc, f_float *work, f_int *info, f_int vect_len);
void dgbbrd_(char *vect, f_int *m, f_int *n, f_int *ncc, f_int *kl, f_int *ku, f_double *ab, f_int *ldab, f_double *d, f_double *e, f_double *q, f_int *ldq, f_double *pt, f_int *ldpt, f_double *c, f_int *ldc, f_double *work, f_int *info, f_int vect_len);
void cgbbrd_(char *vect, f_int *m, f_int *n, f_int *ncc, f_int *kl, f_int *ku, f_cfloat *ab, f_int *ldab, f_float *d, f_float *e, f_cfloat *q, f_int *ldq, f_cfloat *pt, f_int *ldpt, f_cfloat *c, f_int *ldc, f_cfloat *work, f_float *rwork, f_int *info, f_int vect_len);
void zgbbrd_(char *vect, f_int *m, f_int *n, f_int *ncc, f_int *kl, f_int *ku, f_cdouble *ab, f_int *ldab, f_double *d, f_double *e, f_cdouble *q, f_int *ldq, f_cdouble *pt, f_int *ldpt, f_cdouble *c, f_int *ldc, f_cdouble *work, f_double *rwork, f_int *info, f_int vect_len);

/// Estimates the reciprocal of the condition number of a general
/// band matrix, in either the 1-norm or the infinity-norm, using
/// the LU factorization computed by SGBTRF.
void sgbcon_(char *norm, f_int *n, f_int *kl, f_int *ku, f_float *ab, f_int *ldab, f_int *ipiv, f_float *anorm, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int norm_len);
void dgbcon_(char *norm, f_int *n, f_int *kl, f_int *ku, f_double *ab, f_int *ldab, f_int *ipiv, f_double *anorm, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int norm_len);
void cgbcon_(char *norm, f_int *n, f_int *kl, f_int *ku, f_cfloat *ab, f_int *ldab, f_int *ipiv, f_float *anorm, f_float *rcond, f_cfloat *work, f_float *rwork, f_int *info, f_int norm_len);
void zgbcon_(char *norm, f_int *n, f_int *kl, f_int *ku, f_cdouble *ab, f_int *ldab, f_int *ipiv, f_double *anorm, f_double *rcond, f_cdouble *work, f_double *rwork, f_int *info, f_int norm_len);

/// Computes row and column scalings to equilibrate a general band
/// matrix and reduce its condition number.
void sgbequ_(f_int *m, f_int *n, f_int *kl, f_int *ku, f_float *ab, f_int *ldab, f_float *r, f_float *c, f_float *rowcnd, f_float *colcnd, f_float *amax, f_int *info);
void dgbequ_(f_int *m, f_int *n, f_int *kl, f_int *ku, f_double *ab, f_int *ldab, f_double *r, f_double *c, f_double *rowcnd, f_double *colcnd, f_double *amax, f_int *info);
void cgbequ_(f_int *m, f_int *n, f_int *kl, f_int *ku, f_cfloat *ab, f_int *ldab, f_float *r, f_float *c, f_float *rowcnd, f_float *colcnd, f_float *amax, f_int *info);
void zgbequ_(f_int *m, f_int *n, f_int *kl, f_int *ku, f_cdouble *ab, f_int *ldab, f_double *r, f_double *c, f_double *rowcnd, f_double *colcnd, f_double *amax, f_int *info);

/// Improves the computed solution to a general banded system of
/// linear equations AX=B, A**T X=B or A**H X=B, and provides forward
/// and backward error bounds for the solution.
void sgbrfs_(char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_float *ab, f_int *ldab, f_float *afb, f_int *ldafb, f_int *ipiv, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int trans_len);
void dgbrfs_(char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_double *ab, f_int *ldab, f_double *afb, f_int *ldafb, f_int *ipiv, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int trans_len);
void cgbrfs_(char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_cfloat *afb, f_int *ldafb, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int trans_len);
void zgbrfs_(char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_cdouble *afb, f_int *ldafb, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int trans_len);

/// Computes an LU factorization of a general band matrix, using
/// partial pivoting with row interchanges.
void sgbtrf_(f_int *m, f_int *n, f_int *kl, f_int *ku, f_float *ab, f_int *ldab, f_int *ipiv, f_int *info);
void dgbtrf_(f_int *m, f_int *n, f_int *kl, f_int *ku, f_double *ab, f_int *ldab, f_int *ipiv, f_int *info);
void cgbtrf_(f_int *m, f_int *n, f_int *kl, f_int *ku, f_cfloat *ab, f_int *ldab, f_int *ipiv, f_int *info);
void zgbtrf_(f_int *m, f_int *n, f_int *kl, f_int *ku, f_cdouble *ab, f_int *ldab, f_int *ipiv, f_int *info);

/// Solves a general banded system of linear equations AX=B,
/// A**T X=B or A**H X=B, using the LU factorization computed
/// by SGBTRF.
void sgbtrs_(char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_float *ab, f_int *ldab, f_int *ipiv, f_float *b, f_int *ldb, f_int *info, f_int trans_len);
void dgbtrs_(char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_double *ab, f_int *ldab, f_int *ipiv, f_double *b, f_int *ldb, f_int *info, f_int trans_len);
void cgbtrs_(char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int trans_len);
void zgbtrs_(char *trans, f_int *n, f_int *kl, f_int *ku, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int trans_len);

/// Transforms eigenvectors of a balanced matrix to those of the
/// original matrix supplied to SGEBAL.
void sgebak_(char *job, char *side, f_int *n, f_int *ilo, f_int *ihi, f_float *scale, f_int *m, f_float *v, f_int *ldv, f_int *info, f_int job_len, f_int side_len);
void dgebak_(char *job, char *side, f_int *n, f_int *ilo, f_int *ihi, f_double *scale, f_int *m, f_double *v, f_int *ldv, f_int *info, f_int job_len, f_int side_len);
void cgebak_(char *job, char *side, f_int *n, f_int *ilo, f_int *ihi, f_float *scale, f_int *m, f_cfloat *v, f_int *ldv, f_int *info, f_int job_len, f_int side_len);
void zgebak_(char *job, char *side, f_int *n, f_int *ilo, f_int *ihi, f_double *scale, f_int *m, f_cdouble *v, f_int *ldv, f_int *info, f_int job_len, f_int side_len);

/// Balances a general matrix in order to improve the accuracy
/// of computed eigenvalues.
void sgebal_(char *job, f_int *n, f_float *a, f_int *lda, f_int *ilo, f_int *ihi, f_float *scale, f_int *info, f_int job_len);
void dgebal_(char *job, f_int *n, f_double *a, f_int *lda, f_int *ilo, f_int *ihi, f_double *scale, f_int *info, f_int job_len);
void cgebal_(char *job, f_int *n, f_cfloat *a, f_int *lda, f_int *ilo, f_int *ihi, f_float *scale, f_int *info, f_int job_len);
void zgebal_(char *job, f_int *n, f_cdouble *a, f_int *lda, f_int *ilo, f_int *ihi, f_double *scale, f_int *info, f_int job_len);

/// Reduces a general rectangular matrix to real bidiagonal form
/// by an orthogonal transformation.
void sgebrd_(f_int *m, f_int *n, f_float *a, f_int *lda, f_float *d, f_float *e, f_float *tauq, f_float *taup, f_float *work, f_int *lwork, f_int *info);
void dgebrd_(f_int *m, f_int *n, f_double *a, f_int *lda, f_double *d, f_double *e, f_double *tauq, f_double *taup, f_double *work, f_int *lwork, f_int *info);
void cgebrd_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_float *d, f_float *e, f_cfloat *tauq, f_cfloat *taup, f_cfloat *work, f_int *lwork, f_int *info);
void zgebrd_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_double *d, f_double *e, f_cdouble *tauq, f_cdouble *taup, f_cdouble *work, f_int *lwork, f_int *info);

/// Estimates the reciprocal of the condition number of a general
/// matrix, in either the 1-norm or the infinity-norm, using the
/// LU factorization computed by SGETRF.
void sgecon_(char *norm, f_int *n, f_float *a, f_int *lda, f_float *anorm, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int norm_len);
void dgecon_(char *norm, f_int *n, f_double *a, f_int *lda, f_double *anorm, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int norm_len);
void cgecon_(char *norm, f_int *n, f_cfloat *a, f_int *lda, f_float *anorm, f_float *rcond, f_cfloat *work, f_float *rwork, f_int *info, f_int norm_len);
void zgecon_(char *norm, f_int *n, f_cdouble *a, f_int *lda, f_double *anorm, f_double *rcond, f_cdouble *work, f_double *rwork, f_int *info, f_int norm_len);

/// Computes row and column scalings to equilibrate a general
/// rectangular matrix and reduce its condition number.
void sgeequ_(f_int *m, f_int *n, f_float *a, f_int *lda, f_float *r, f_float *c, f_float *rowcnd, f_float *colcnd, f_float *amax, f_int *info);
void dgeequ_(f_int *m, f_int *n, f_double *a, f_int *lda, f_double *r, f_double *c, f_double *rowcnd, f_double *colcnd, f_double *amax, f_int *info);
void cgeequ_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_float *r, f_float *c, f_float *rowcnd, f_float *colcnd, f_float *amax, f_int *info);
void zgeequ_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_double *r, f_double *c, f_double *rowcnd, f_double *colcnd, f_double *amax, f_int *info);

/// Reduces a general matrix to upper Hessenberg form by an
/// orthogonal similarity transformation.
void sgehrd_(f_int *n, f_int *ilo, f_int *ihi, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dgehrd_(f_int *n, f_int *ilo, f_int *ihi, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);
void cgehrd_(f_int *n, f_int *ilo, f_int *ihi, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zgehrd_(f_int *n, f_int *ilo, f_int *ihi, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Computes an LQ factorization of a general rectangular matrix.
void sgelqf_(f_int *m, f_int *n, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dgelqf_(f_int *m, f_int *n, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);
void cgelqf_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zgelqf_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Computes a QL factorization of a general rectangular matrix.
void sgeqlf_(f_int *m, f_int *n, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dgeqlf_(f_int *m, f_int *n, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);
void cgeqlf_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zgeqlf_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Computes a QR factorization with column pivoting of a general
/// rectangular matrix using Level 3 BLAS.
void sgeqp3_(f_int *m, f_int *n, f_float *a, f_int *lda, f_int *jpvt, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dgeqp3_(f_int *m, f_int *n, f_double *a, f_int *lda, f_int *jpvt, f_double *tau, f_double *work, f_int *lwork, f_int *info);
void cgeqp3_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_int *jpvt, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info);
void zgeqp3_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_int *jpvt, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info);

/// Computes a QR factorization with column pivoting of a general
/// rectangular matrix.
void sgeqpf_(f_int *m, f_int *n, f_float *a, f_int *lda, f_int *jpvt, f_float *tau, f_float *work, f_int *info);
void dgeqpf_(f_int *m, f_int *n, f_double *a, f_int *lda, f_int *jpvt, f_double *tau, f_double *work, f_int *info);
void cgeqpf_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_int *jpvt, f_cfloat *tau, f_cfloat *work, f_float *rwork, f_int *info);
void zgeqpf_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_int *jpvt, f_cdouble *tau, f_cdouble *work, f_double *rwork, f_int *info);

/// Computes a QR factorization of a general rectangular matrix.
void sgeqrf_(f_int *m, f_int *n, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dgeqrf_(f_int *m, f_int *n, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);
void cgeqrf_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zgeqrf_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Improves the computed solution to a general system of linear
/// equations AX=B, A**T X=B or A**H X=B, and provides forward and
/// backward error bounds for the solution.
void sgerfs_(char *trans, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *af, f_int *ldaf, f_int *ipiv, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int trans_len);
void dgerfs_(char *trans, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *af, f_int *ldaf, f_int *ipiv, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int trans_len);
void cgerfs_(char *trans, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *af, f_int *ldaf, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int trans_len);
void zgerfs_(char *trans, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *af, f_int *ldaf, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int trans_len);

/// Computes an RQ factorization of a general rectangular matrix.
void sgerqf_(f_int *m, f_int *n, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dgerqf_(f_int *m, f_int *n, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);
void cgerqf_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zgerqf_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Computes an LU factorization of a general matrix, using partial
/// pivoting with row interchanges.
void sgetrf_(f_int *m, f_int *n, f_float *a, f_int *lda, f_int *ipiv, f_int *info);
void dgetrf_(f_int *m, f_int *n, f_double *a, f_int *lda, f_int *ipiv, f_int *info);
void cgetrf_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_int *ipiv, f_int *info);
void zgetrf_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_int *ipiv, f_int *info);

/// Computes the inverse of a general matrix, using the LU factorization
/// computed by SGETRF.
void sgetri_(f_int *n, f_float *a, f_int *lda, f_int *ipiv, f_float *work, f_int *lwork, f_int *info);
void dgetri_(f_int *n, f_double *a, f_int *lda, f_int *ipiv, f_double *work, f_int *lwork, f_int *info);
void cgetri_(f_int *n, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *work, f_int *lwork, f_int *info);
void zgetri_(f_int *n, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *work, f_int *lwork, f_int *info);

/// Solves a general system of linear equations AX=B, A**T X=B
/// or A**H X=B, using the LU factorization computed by SGETRF.
void sgetrs_(char *trans, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_int *ipiv, f_float *b, f_int *ldb, f_int *info, f_int trans_len);
void dgetrs_(char *trans, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_int *ipiv, f_double *b, f_int *ldb, f_int *info, f_int trans_len);
void cgetrs_(char *trans, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int trans_len);
void zgetrs_(char *trans, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int trans_len);

/// Forms the right or left eigenvectors of the generalized eigenvalue
/// problem by backward transformation on the computed eigenvectors of
/// the balanced pair of matrices output by SGGBAL.
void sggbak_(char *job, char *side, f_int *n, f_int *ilo, f_int *ihi, f_float *lscale, f_float *rscale, f_int *m, f_float *v, f_int *ldv, f_int *info, f_int job_len, f_int side_len);
void dggbak_(char *job, char *side, f_int *n, f_int *ilo, f_int *ihi, f_double *lscale, f_double *rscale, f_int *m, f_double *v, f_int *ldv, f_int *info, f_int job_len, f_int side_len);
void cggbak_(char *job, char *side, f_int *n, f_int *ilo, f_int *ihi, f_float *lscale, f_float *rscale, f_int *m, f_cfloat *v, f_int *ldv, f_int *info, f_int job_len, f_int side_len);
void zggbak_(char *job, char *side, f_int *n, f_int *ilo, f_int *ihi, f_double *lscale, f_double *rscale, f_int *m, f_cdouble *v, f_int *ldv, f_int *info, f_int job_len, f_int side_len);

/// Balances a pair of general real matrices for the generalized
/// eigenvalue problem A x = lambda B x.
void sggbal_(char *job, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *ilo, f_int *ihi, f_float *lscale, f_float *rscale, f_float *work, f_int *info, f_int job_len);
void dggbal_(char *job, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *ilo, f_int *ihi, f_double *lscale, f_double *rscale, f_double *work, f_int *info, f_int job_len);
void cggbal_(char *job, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *ilo, f_int *ihi, f_float *lscale, f_float *rscale, f_float *work, f_int *info, f_int job_len);
void zggbal_(char *job, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *ilo, f_int *ihi, f_double *lscale, f_double *rscale, f_double *work, f_int *info, f_int job_len);

/// Reduces a pair of real matrices to generalized upper
/// Hessenberg form using orthogonal transformations 
void sgghrd_(char *compq, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *q, f_int *ldq, f_float *z, f_int *ldz, f_int *info, f_int compq_len, f_int compz_len);
void dgghrd_(char *compq, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *q, f_int *ldq, f_double *z, f_int *ldz, f_int *info, f_int compq_len, f_int compz_len);
void cgghrd_(char *compq, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *q, f_int *ldq, f_cfloat *z, f_int *ldz, f_int *info, f_int compq_len, f_int compz_len);
void zgghrd_(char *compq, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *q, f_int *ldq, f_cdouble *z, f_int *ldz, f_int *info, f_int compq_len, f_int compz_len);

/// Computes a generalized QR factorization of a pair of matrices. 
void sggqrf_(f_int *n, f_int *m, f_int *p, f_float *a, f_int *lda, f_float *taua, f_float *b, f_int *ldb, f_float *taub, f_float *work, f_int *lwork, f_int *info);
void dggqrf_(f_int *n, f_int *m, f_int *p, f_double *a, f_int *lda, f_double *taua, f_double *b, f_int *ldb, f_double *taub, f_double *work, f_int *lwork, f_int *info);
void cggqrf_(f_int *n, f_int *m, f_int *p, f_cfloat *a, f_int *lda, f_cfloat *taua, f_cfloat *b, f_int *ldb, f_cfloat *taub, f_cfloat *work, f_int *lwork, f_int *info);
void zggqrf_(f_int *n, f_int *m, f_int *p, f_cdouble *a, f_int *lda, f_cdouble *taua, f_cdouble *b, f_int *ldb, f_cdouble *taub, f_cdouble *work, f_int *lwork, f_int *info);

/// Computes a generalized RQ factorization of a pair of matrices.
void sggrqf_(f_int *m, f_int *p, f_int *n, f_float *a, f_int *lda, f_float *taua, f_float *b, f_int *ldb, f_float *taub, f_float *work, f_int *lwork, f_int *info);
void dggrqf_(f_int *m, f_int *p, f_int *n, f_double *a, f_int *lda, f_double *taua, f_double *b, f_int *ldb, f_double *taub, f_double *work, f_int *lwork, f_int *info);
void cggrqf_(f_int *m, f_int *p, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *taua, f_cfloat *b, f_int *ldb, f_cfloat *taub, f_cfloat *work, f_int *lwork, f_int *info);
void zggrqf_(f_int *m, f_int *p, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *taua, f_cdouble *b, f_int *ldb, f_cdouble *taub, f_cdouble *work, f_int *lwork, f_int *info);

/// Computes orthogonal matrices as a preprocessing step
/// for computing the generalized singular value decomposition
void sggsvp_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *p, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *tola, f_float *tolb, f_int *k, f_int *l, f_float *u, f_int *ldu, f_float *v, f_int *ldv, f_float *q, f_int *ldq, f_int *iwork, f_float *tau, f_float *work, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void dggsvp_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *p, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *tola, f_double *tolb, f_int *k, f_int *l, f_double *u, f_int *ldu, f_double *v, f_int *ldv, f_double *q, f_int *ldq, f_int *iwork, f_double *tau, f_double *work, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void cggsvp_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *p, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_float *tola, f_float *tolb, f_int *k, f_int *l, f_cfloat *u, f_int *ldu, f_cfloat *v, f_int *ldv, f_cfloat *q, f_int *ldq, f_int *iwork, f_float *rwork, f_cfloat *tau, f_cfloat *work, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void zggsvp_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *p, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_double *tola, f_double *tolb, f_int *k, f_int *l, f_cdouble *u, f_int *ldu, f_cdouble *v, f_int *ldv, f_cdouble *q, f_int *ldq, f_int *iwork, f_double *rwork, f_cdouble *tau, f_cdouble *work, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);

/// Estimates the reciprocal of the condition number of a general
/// tridiagonal matrix, in either the 1-norm or the infinity-norm,
/// using the LU factorization computed by SGTTRF.
void sgtcon_(char *norm, f_int *n, f_float *dl, f_float *d, f_float *du, f_float *du2, f_int *ipiv, f_float *anorm, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int norm_len);
void dgtcon_(char *norm, f_int *n, f_double *dl, f_double *d, f_double *du, f_double *du2, f_int *ipiv, f_double *anorm, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int norm_len);
void cgtcon_(char *norm, f_int *n, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *du2, f_int *ipiv, f_float *anorm, f_float *rcond, f_cfloat *work, f_int *info, f_int norm_len);
void zgtcon_(char *norm, f_int *n, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *du2, f_int *ipiv, f_double *anorm, f_double *rcond, f_cdouble *work, f_int *info, f_int norm_len);

/// Improves the computed solution to a general tridiagonal system
/// of linear equations AX=B, A**T X=B or A**H X=B, and provides
/// forward and backward error bounds for the solution.
void sgtrfs_(char *trans, f_int *n, f_int *nrhs, f_float *dl, f_float *d, f_float *du, f_float *dlf, f_float *df, f_float *duf, f_float *du2, f_int *ipiv, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int trans_len);
void dgtrfs_(char *trans, f_int *n, f_int *nrhs, f_double *dl, f_double *d, f_double *du, f_double *dlf, f_double *df, f_double *duf, f_double *du2, f_int *ipiv, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int trans_len);
void cgtrfs_(char *trans, f_int *n, f_int *nrhs, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *dlf, f_cfloat *df, f_cfloat *duf, f_cfloat *du2, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int trans_len);
void zgtrfs_(char *trans, f_int *n, f_int *nrhs, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *dlf, f_cdouble *df, f_cdouble *duf, f_cdouble *du2, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int trans_len);

/// Computes an LU factorization of a general tridiagonal matrix,
/// using partial pivoting with row interchanges.
void sgttrf_(f_int *n, f_float *dl, f_float *d, f_float *du, f_float *du2, f_int *ipiv, f_int *info);
void dgttrf_(f_int *n, f_double *dl, f_double *d, f_double *du, f_double *du2, f_int *ipiv, f_int *info);
void cgttrf_(f_int *n, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *du2, f_int *ipiv, f_int *info);
void zgttrf_(f_int *n, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *du2, f_int *ipiv, f_int *info);

/// Solves a general tridiagonal system of linear equations AX=B,
/// A**T X=B or A**H X=B, using the LU factorization computed by
/// SGTTRF.
void sgttrs_(char *trans, f_int *n, f_int *nrhs, f_float *dl, f_float *d, f_float *du, f_float *du2, f_int *ipiv, f_float *b, f_int *ldb, f_int *info, f_int trans_len);
void dgttrs_(char *trans, f_int *n, f_int *nrhs, f_double *dl, f_double *d, f_double *du, f_double *du2, f_int *ipiv, f_double *b, f_int *ldb, f_int *info, f_int trans_len);
void cgttrs_(char *trans, f_int *n, f_int *nrhs, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *du2, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int trans_len);
void zgttrs_(char *trans, f_int *n, f_int *nrhs, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *du2, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int trans_len);

/// Implements a single-/f_double-shift version of the QZ method for
/// finding the generalized eigenvalues of the equation 
/// det(A - w(i) B) = 0
void shgeqz_(char *job, char *compq, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *q, f_int *ldq, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *info, f_int job_len, f_int compq_len, f_int compz_len);
void dhgeqz_(char *job, char *compq, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *q, f_int *ldq, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *info, f_int job_len, f_int compq_len, f_int compz_len);
void chgeqz_(char *job, char *compq, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *q, f_int *ldq, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *info, f_int job_len, f_int compq_len, f_int compz_len);
void zhgeqz_(char *job, char *compq, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *q, f_int *ldq, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *info, f_int job_len, f_int compq_len, f_int compz_len);

/// Computes specified right and/or left eigenvectors of an upper
/// Hessenberg matrix by inverse iteration.
void shsein_(char *side, char *eigsrc, char *initv, f_int *select, f_int *n, f_float *h, f_int *ldh, f_float *wr, f_float *wi, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_int *mm, f_int *m, f_float *work, f_int *ifaill, f_int *ifailr, f_int *info, f_int side_len, f_int eigsrc_len, f_int initv_len);
void dhsein_(char *side, char *eigsrc, char *initv, f_int *select, f_int *n, f_double *h, f_int *ldh, f_double *wr, f_double *wi, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_int *mm, f_int *m, f_double *work, f_int *ifaill, f_int *ifailr, f_int *info, f_int side_len, f_int eigsrc_len, f_int initv_len);
void chsein_(char *side, char *eigsrc, char *initv, f_int *select, f_int *n, f_cfloat *h, f_int *ldh, f_cfloat *w, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_int *mm, f_int *m, f_cfloat *work, f_float *rwork, f_int *ifaill, f_int *ifailr, f_int *info, f_int side_len, f_int eigsrc_len, f_int initv_len);
void zhsein_(char *side, char *eigsrc, char *initv, f_int *select, f_int *n, f_cdouble *h, f_int *ldh, f_cdouble *w, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_int *mm, f_int *m, f_cdouble *work, f_double *rwork, f_int *ifaill, f_int *ifailr, f_int *info, f_int side_len, f_int eigsrc_len, f_int initv_len);

/// Computes the eigenvalues and Schur factorization of an upper
/// Hessenberg matrix, using the multishift QR algorithm.
void shseqr_(char *job, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_float *h, f_int *ldh, f_float *wr, f_float *wi, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *info, f_int job_len, f_int compz_len);
void dhseqr_(char *job, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_double *h, f_int *ldh, f_double *wr, f_double *wi, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *info, f_int job_len, f_int compz_len);
void chseqr_(char *job, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_cfloat *h, f_int *ldh, f_cfloat *w, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_int *info, f_int job_len, f_int compz_len);
void zhseqr_(char *job, char *compz, f_int *n, f_int *ilo, f_int *ihi, f_cdouble *h, f_int *ldh, f_cdouble *w, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_int *info, f_int job_len, f_int compz_len);

/// Generates the orthogonal transformation matrix from
/// a reduction to tridiagonal form determined by SSPTRD.
void sopgtr_(char *uplo, f_int *n, f_float *ap, f_float *tau, f_float *q, f_int *ldq, f_float *work, f_int *info, f_int uplo_len);
void dopgtr_(char *uplo, f_int *n, f_double *ap, f_double *tau, f_double *q, f_int *ldq, f_double *work, f_int *info, f_int uplo_len);

/// Generates the unitary transformation matrix from
/// a reduction to tridiagonal form determined by CHPTRD.
void cupgtr_(char *uplo, f_int *n, f_cfloat *ap, f_cfloat *tau, f_cfloat *q, f_int *ldq, f_cfloat *work, f_int *info, f_int uplo_len);
void zupgtr_(char *uplo, f_int *n, f_cdouble *ap, f_cdouble *tau, f_cdouble *q, f_int *ldq, f_cdouble *work, f_int *info, f_int uplo_len);


/// Multiplies a general matrix by the orthogonal
/// transformation matrix from a reduction to tridiagonal form
/// determined by SSPTRD.
void sopmtr_(char *side, char *uplo, char *trans, f_int *m, f_int *n, f_float *ap, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *info, f_int side_len, f_int uplo_len, f_int trans_len);
void dopmtr_(char *side, char *uplo, char *trans, f_int *m, f_int *n, f_double *ap, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *info, f_int side_len, f_int uplo_len, f_int trans_len);

/// Generates the orthogonal transformation matrices from
/// a reduction to bidiagonal form determined by SGEBRD.
void sorgbr_(char *vect, f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info, f_int vect_len);
void dorgbr_(char *vect, f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info, f_int vect_len);

/// Generates the unitary transformation matrices from
/// a reduction to bidiagonal form determined by CGEBRD.
void cungbr_(char *vect, f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info, f_int vect_len);
void zungbr_(char *vect, f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info, f_int vect_len);

/// Generates the orthogonal transformation matrix from
/// a reduction to Hessenberg form determined by SGEHRD.
void sorghr_(f_int *n, f_int *ilo, f_int *ihi, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dorghr_(f_int *n, f_int *ilo, f_int *ihi, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);

/// Generates the unitary transformation matrix from
/// a reduction to Hessenberg form determined by CGEHRD.
void cunghr_(f_int *n, f_int *ilo, f_int *ihi, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zunghr_(f_int *n, f_int *ilo, f_int *ihi, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Generates all or part of the orthogonal matrix Q from
/// an LQ factorization determined by SGELQF.
void sorglq_(f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dorglq_(f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);

/// Generates all or part of the unitary matrix Q from
/// an LQ factorization determined by CGELQF.
void cunglq_(f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zunglq_(f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Generates all or part of the orthogonal matrix Q from
/// a QL factorization determined by SGEQLF.
void sorgql_(f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dorgql_(f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);

/// Generates all or part of the unitary matrix Q from
/// a QL factorization determined by CGEQLF.
void cungql_(f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zungql_(f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Generates all or part of the orthogonal matrix Q from
/// a QR factorization determined by SGEQRF.
void sorgqr_(f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dorgqr_(f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);

/// Generates all or part of the unitary matrix Q from
/// a QR factorization determined by CGEQRF.
void cungqr_(f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zungqr_(f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Generates all or part of the orthogonal matrix Q from
/// an RQ factorization determined by SGERQF.
void sorgrq_(f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dorgrq_(f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);

/// Generates all or part of the unitary matrix Q from
/// an RQ factorization determined by CGERQF.
void cungrq_(f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void zungrq_(f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);

/// Generates the orthogonal transformation matrix from
/// a reduction to tridiagonal form determined by SSYTRD.
void sorgtr_(char *uplo, f_int *n, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info, f_int uplo_len);
void dorgtr_(char *uplo, f_int *n, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info, f_int uplo_len);

/// Generates the unitary transformation matrix from
/// a reduction to tridiagonal form determined by CHETRD.
void cungtr_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info, f_int uplo_len);
void zungtr_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info, f_int uplo_len);

/// Multiplies a general matrix by one of the orthogonal
/// transformation  matrices from a reduction to bidiagonal form
/// determined by SGEBRD.
void sormbr_(char *vect, char *side, char *trans, f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *lwork, f_int *info, f_int vect_len, f_int side_len, f_int trans_len);
void dormbr_(char *vect, char *side, char *trans, f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *lwork, f_int *info, f_int vect_len, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by one of the unitary
/// transformation matrices from a reduction to bidiagonal form
/// determined by CGEBRD.
void cunmbr_(char *vect, char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *lwork, f_int *info, f_int vect_len, f_int side_len, f_int trans_len);
void zunmbr_(char *vect, char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *lwork, f_int *info, f_int vect_len, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the orthogonal transformation
/// matrix from a reduction to Hessenberg form determined by SGEHRD.
void sormhr_(char *side, char *trans, f_int *m, f_int *n, f_int *ilo, f_int *ihi, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void dormhr_(char *side, char *trans, f_int *m, f_int *n, f_int *ilo, f_int *ihi, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the unitary transformation
/// matrix from a reduction to Hessenberg form determined by CGEHRD.
void cunmhr_(char *side, char *trans, f_int *m, f_int *n, f_int *ilo, f_int *ihi, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void zunmhr_(char *side, char *trans, f_int *m, f_int *n, f_int *ilo, f_int *ihi, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the orthogonal matrix
/// from an LQ factorization determined by SGELQF.
void sormlq_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void dormlq_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the unitary matrix
/// from an LQ factorization determined by CGELQF.
void cunmlq_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void zunmlq_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the orthogonal matrix
/// from a QL factorization determined by SGEQLF.
void sormql_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void dormql_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the unitary matrix
/// from a QL factorization determined by CGEQLF.
void cunmql_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void zunmql_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the orthogonal matrix
/// from a QR factorization determined by SGEQRF.
void sormqr_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void dormqr_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the unitary matrix
/// from a QR factorization determined by CGEQRF.
void cunmqr_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void zunmqr_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiples a general matrix by the orthogonal matrix
/// from an RZ factorization determined by STZRZF.
void sormr3_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_int *l, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *info, f_int side_len, f_int trans_len);
void dormr3_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_int *l, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *info, f_int side_len, f_int trans_len);

/// Multiples a general matrix by the unitary matrix
/// from an RZ factorization determined by CTZRZF.
void cunmr3_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_int *l, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *info, f_int side_len, f_int trans_len);
void zunmr3_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_int *l, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the orthogonal matrix
/// from an RQ factorization determined by SGERQF.
void sormrq_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void dormrq_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the unitary matrix
/// from an RQ factorization determined by CGERQF.
void cunmrq_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void zunmrq_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiples a general matrix by the orthogonal matrix
/// from an RZ factorization determined by STZRZF.
void sormrz_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_int *l, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void dormrz_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_int *l, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiples a general matrix by the unitary matrix
/// from an RZ factorization determined by CTZRZF.
void cunmrz_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_int *l, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);
void zunmrz_(char *side, char *trans, f_int *m, f_int *n, f_int *k, f_int *l, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *lwork, f_int *info, f_int side_len, f_int trans_len);

/// Multiplies a general matrix by the orthogonal
/// transformation matrix from a reduction to tridiagonal form
/// determined by SSYTRD.
void sormtr_(char *side, char *uplo, char *trans, f_int *m, f_int *n, f_float *a, f_int *lda, f_float *tau, f_float *c, f_int *ldc, f_float *work, f_int *lwork, f_int *info, f_int side_len, f_int uplo_len, f_int trans_len);
void dormtr_(char *side, char *uplo, char *trans, f_int *m, f_int *n, f_double *a, f_int *lda, f_double *tau, f_double *c, f_int *ldc, f_double *work, f_int *lwork, f_int *info, f_int side_len, f_int uplo_len, f_int trans_len);

/// Multiplies a general matrix by the unitary
/// transformation matrix from a reduction to tridiagonal form
/// determined by CHETRD.
void cunmtr_(char *side, char *uplo, char *trans, f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *lwork, f_int *info, f_int side_len, f_int uplo_len, f_int trans_len);
void zunmtr_(char *side, char *uplo, char *trans, f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *lwork, f_int *info, f_int side_len, f_int uplo_len, f_int trans_len);

/// Estimates the reciprocal of the condition number of a
/// symmetric positive definite band matrix, using the
/// Cholesky factorization computed by SPBTRF.
void spbcon_(char *uplo, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_float *anorm, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dpbcon_(char *uplo, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_double *anorm, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void cpbcon_(char *uplo, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_float *anorm, f_float *rcond, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zpbcon_(char *uplo, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_double *anorm, f_double *rcond, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Computes row and column scalings to equilibrate a symmetric
/// positive definite band matrix and reduce its condition number.
void spbequ_(char *uplo, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_float *s, f_float *scond, f_float *amax, f_int *info, f_int uplo_len);
void dpbequ_(char *uplo, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_double *s, f_double *scond, f_double *amax, f_int *info, f_int uplo_len);
void cpbequ_(char *uplo, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_float *s, f_float *scond, f_float *amax, f_int *info, f_int uplo_len);
void zpbequ_(char *uplo, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_double *s, f_double *scond, f_double *amax, f_int *info, f_int uplo_len);

/// Improves the computed solution to a symmetric positive
/// definite banded system of linear equations AX=B, and provides
/// forward and backward error bounds for the solution.
void spbrfs_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_float *ab, f_int *ldab, f_float *afb, f_int *ldafb, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dpbrfs_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_double *ab, f_int *ldab, f_double *afb, f_int *ldafb, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void cpbrfs_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_cfloat *afb, f_int *ldafb, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zpbrfs_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_cdouble *afb, f_int *ldafb, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Computes a split Cholesky factorization of a real symmetric positive
/// definite band matrix.
void spbstf_(char *uplo, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_int *info, f_int uplo_len);
void dpbstf_(char *uplo, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_int *info, f_int uplo_len);
void cpbstf_(char *uplo, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_int *info, f_int uplo_len);
void zpbstf_(char *uplo, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_int *info, f_int uplo_len);

/// Computes the Cholesky factorization of a symmetric
/// positive definite band matrix.
void spbtrf_(char *uplo, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_int *info, f_int uplo_len);
void dpbtrf_(char *uplo, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_int *info, f_int uplo_len);
void cpbtrf_(char *uplo, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_int *info, f_int uplo_len);
void zpbtrf_(char *uplo, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_int *info, f_int uplo_len);

/// Solves a symmetric positive definite banded system
/// of linear equations AX=B, using the Cholesky factorization
/// computed by SPBTRF.
void spbtrs_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_float *ab, f_int *ldab, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dpbtrs_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_double *ab, f_int *ldab, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void cpbtrs_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zpbtrs_(char *uplo, f_int *n, f_int *kd, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Estimates the reciprocal of the condition number of a
/// symmetric positive definite matrix, using the
/// Cholesky factorization computed by SPOTRF.
void spocon_(char *uplo, f_int *n, f_float *a, f_int *lda, f_float *anorm, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dpocon_(char *uplo, f_int *n, f_double *a, f_int *lda, f_double *anorm, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void cpocon_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_float *anorm, f_float *rcond, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zpocon_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_double *anorm, f_double *rcond, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Computes row and column scalings to equilibrate a symmetric
/// positive definite matrix and reduce its condition number.
void spoequ_(f_int *n, f_float *a, f_int *lda, f_float *s, f_float *scond, f_float *amax, f_int *info);
void dpoequ_(f_int *n, f_double *a, f_int *lda, f_double *s, f_double *scond, f_double *amax, f_int *info);
void cpoequ_(f_int *n, f_cfloat *a, f_int *lda, f_float *s, f_float *scond, f_float *amax, f_int *info);
void zpoequ_(f_int *n, f_cdouble *a, f_int *lda, f_double *s, f_double *scond, f_double *amax, f_int *info);

/// Improves the computed solution to a symmetric positive
/// definite system of linear equations AX=B, and provides forward
/// and backward error bounds for the solution.
void sporfs_(char *uplo, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *af, f_int *ldaf, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dporfs_(char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *af, f_int *ldaf, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void cporfs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *af, f_int *ldaf, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zporfs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *af, f_int *ldaf, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Computes the Cholesky factorization of a symmetric
/// positive definite matrix.
void spotrf_(char *uplo, f_int *n, f_float *a, f_int *lda, f_int *info, f_int uplo_len);
void dpotrf_(char *uplo, f_int *n, f_double *a, f_int *lda, f_int *info, f_int uplo_len);
void cpotrf_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_int *info, f_int uplo_len);
void zpotrf_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_int *info, f_int uplo_len);

/// Computes the inverse of a symmetric positive definite
/// matrix, using the Cholesky factorization computed by SPOTRF.
void spotri_(char *uplo, f_int *n, f_float *a, f_int *lda, f_int *info, f_int uplo_len);
void dpotri_(char *uplo, f_int *n, f_double *a, f_int *lda, f_int *info, f_int uplo_len);
void cpotri_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_int *info, f_int uplo_len);
void zpotri_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_int *info, f_int uplo_len);

/// Solves a symmetric positive definite system of linear
/// equations AX=B, using the Cholesky factorization computed by
/// SPOTRF.
void spotrs_(char *uplo, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dpotrs_(char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void cpotrs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zpotrs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Estimates the reciprocal of the condition number of a
/// symmetric positive definite matrix in packed storage,
/// using the Cholesky factorization computed by SPPTRF.
void sppcon_(char *uplo, f_int *n, f_float *ap, f_float *anorm, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dppcon_(char *uplo, f_int *n, f_double *ap, f_double *anorm, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void cppcon_(char *uplo, f_int *n, f_cfloat *ap, f_float *anorm, f_float *rcond, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zppcon_(char *uplo, f_int *n, f_cdouble *ap, f_double *anorm, f_double *rcond, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Computes row and column scalings to equilibrate a symmetric
/// positive definite matrix in packed storage and reduce its condition
/// number.
void sppequ_(char *uplo, f_int *n, f_float *ap, f_float *s, f_float *scond, f_float *amax, f_int *info, f_int uplo_len);
void dppequ_(char *uplo, f_int *n, f_double *ap, f_double *s, f_double *scond, f_double *amax, f_int *info, f_int uplo_len);
void cppequ_(char *uplo, f_int *n, f_cfloat *ap, f_float *s, f_float *scond, f_float *amax, f_int *info, f_int uplo_len);
void zppequ_(char *uplo, f_int *n, f_cdouble *ap, f_double *s, f_double *scond, f_double *amax, f_int *info, f_int uplo_len);

/// Improves the computed solution to a symmetric positive
/// definite system of linear equations AX=B, where A is held in
/// packed storage, and provides forward and backward error bounds
/// for the solution.
void spprfs_(char *uplo, f_int *n, f_int *nrhs, f_float *ap, f_float *afp, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dpprfs_(char *uplo, f_int *n, f_int *nrhs, f_double *ap, f_double *afp, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void cpprfs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *afp, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zpprfs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *afp, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Computes the Cholesky factorization of a symmetric
/// positive definite matrix in packed storage.
void spptrf_(char *uplo, f_int *n, f_float *ap, f_int *info, f_int uplo_len);
void dpptrf_(char *uplo, f_int *n, f_double *ap, f_int *info, f_int uplo_len);
void cpptrf_(char *uplo, f_int *n, f_cfloat *ap, f_int *info, f_int uplo_len);
void zpptrf_(char *uplo, f_int *n, f_cdouble *ap, f_int *info, f_int uplo_len);

/// Computes the inverse of a symmetric positive definite
/// matrix in packed storage, using the Cholesky factorization computed
/// by SPPTRF.
void spptri_(char *uplo, f_int *n, f_float *ap, f_int *info, f_int uplo_len);
void dpptri_(char *uplo, f_int *n, f_double *ap, f_int *info, f_int uplo_len);
void cpptri_(char *uplo, f_int *n, f_cfloat *ap, f_int *info, f_int uplo_len);
void zpptri_(char *uplo, f_int *n, f_cdouble *ap, f_int *info, f_int uplo_len);

/// Solves a symmetric positive definite system of linear
/// equations AX=B, where A is held in packed storage, using the
/// Cholesky factorization computed by SPPTRF.
void spptrs_(char *uplo, f_int *n, f_int *nrhs, f_float *ap, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dpptrs_(char *uplo, f_int *n, f_int *nrhs, f_double *ap, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void cpptrs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zpptrs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Computes the reciprocal of the condition number of a
/// symmetric positive definite tridiagonal matrix,
/// using the LDL**H factorization computed by SPTTRF.
void sptcon_(f_int *n, f_float *d, f_float *e, f_float *anorm, f_float *rcond, f_float *work, f_int *info);
void dptcon_(f_int *n, f_double *d, f_double *e, f_double *anorm, f_double *rcond, f_double *work, f_int *info);
void cptcon_(f_int *n, f_float *d, f_cfloat *e, f_float *anorm, f_float *rcond, f_float *rwork, f_int *info);
void zptcon_(f_int *n, f_double *d, f_cdouble *e, f_double *anorm, f_double *rcond, f_double *rwork, f_int *info);

/// Computes all eigenvalues and eigenvectors of a real symmetric
/// positive definite tridiagonal matrix, by computing the SVD of
/// its bidiagonal Cholesky factor.
void spteqr_(char *compz, f_int *n, f_float *d, f_float *e, f_float *z, f_int *ldz, f_float *work, f_int *info, f_int compz_len);
void dpteqr_(char *compz, f_int *n, f_double *d, f_double *e, f_double *z, f_int *ldz, f_double *work, f_int *info, f_int compz_len);
void cpteqr_(char *compz, f_int *n, f_float *d, f_float *e, f_cfloat *z, f_int *ldz, f_float *work, f_int *info, f_int compz_len);
void zpteqr_(char *compz, f_int *n, f_double *d, f_double *e, f_cdouble *z, f_int *ldz, f_double *work, f_int *info, f_int compz_len);

/// Improves the computed solution to a symmetric positive
/// definite tridiagonal system of linear equations AX=B, and provides
/// forward and backward error bounds for the solution.
void sptrfs_(f_int *n, f_int *nrhs, f_float *d, f_float *e, f_float *df, f_float *ef, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *info);
void dptrfs_(f_int *n, f_int *nrhs, f_double *d, f_double *e, f_double *df, f_double *ef, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *info);
void cptrfs_(char *uplo, f_int *n, f_int *nrhs, f_float *d, f_cfloat *e, f_float *df, f_cfloat *ef, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zptrfs_(char *uplo, f_int *n, f_int *nrhs, f_double *d, f_cdouble *e, f_double *df, f_cdouble *ef, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Computes the LDL**H factorization of a symmetric
/// positive definite tridiagonal matrix.
void spttrf_(f_int *n, f_float *d, f_float *e, f_int *info);
void dpttrf_(f_int *n, f_double *d, f_double *e, f_int *info);
void cpttrf_(f_int *n, f_float *d, f_cfloat *e, f_int *info);
void zpttrf_(f_int *n, f_double *d, f_cdouble *e, f_int *info);

/// Solves a symmetric positive definite tridiagonal
/// system of linear equations, using the LDL**H factorization
/// computed by SPTTRF.
void spttrs_(f_int *n, f_int *nrhs, f_float *d, f_float *e, f_float *b, f_int *ldb, f_int *info);
void dpttrs_(f_int *n, f_int *nrhs, f_double *d, f_double *e, f_double *b, f_int *ldb, f_int *info);
void cpttrs_(char *uplo, f_int *n, f_int *nrhs, f_float *d, f_cfloat *e, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zpttrs_(char *uplo, f_int *n, f_int *nrhs, f_double *d, f_cdouble *e, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Reduces a real symmetric-definite banded generalized eigenproblem
/// A x = lambda B x to standard form, where B has been factorized by
/// SPBSTF (Crawford's algorithm).
void ssbgst_(char *vect, char *uplo, f_int *n, f_int *ka, f_int *kb, f_float *ab, f_int *ldab, f_float *bb, f_int *ldbb, f_float *x, f_int *ldx, f_float *work, f_int *info, f_int vect_len, f_int uplo_len);
void dsbgst_(char *vect, char *uplo, f_int *n, f_int *ka, f_int *kb, f_double *ab, f_int *ldab, f_double *bb, f_int *ldbb, f_double *x, f_int *ldx, f_double *work, f_int *info, f_int vect_len, f_int uplo_len);

/// Reduces a complex Hermitian-definite banded generalized eigenproblem
/// A x = lambda B x to standard form, where B has been factorized by
/// CPBSTF (Crawford's algorithm).
void chbgst_(char *vect, char *uplo, f_int *n, f_int *ka, f_int *kb, f_cfloat *ab, f_int *ldab, f_cfloat *bb, f_int *ldbb, f_cfloat *x, f_int *ldx, f_cfloat *work, f_float *rwork, f_int *info, f_int vect_len, f_int uplo_len);
void zhbgst_(char *vect, char *uplo, f_int *n, f_int *ka, f_int *kb, f_cdouble *ab, f_int *ldab, f_cdouble *bb, f_int *ldbb, f_cdouble *x, f_int *ldx, f_cdouble *work, f_double *rwork, f_int *info, f_int vect_len, f_int uplo_len);

/// Reduces a symmetric band matrix to real symmetric
/// tridiagonal form by an orthogonal similarity transformation.
void ssbtrd_(char *vect, char *uplo, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_float *d, f_float *e, f_float *q, f_int *ldq, f_float *work, f_int *info, f_int vect_len, f_int uplo_len);
void dsbtrd_(char *vect, char *uplo, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_double *d, f_double *e, f_double *q, f_int *ldq, f_double *work, f_int *info, f_int vect_len, f_int uplo_len);

/// Reduces a Hermitian band matrix to real symmetric
/// tridiagonal form by a unitary similarity transformation.
void chbtrd_(char *vect, char *uplo, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_float *d, f_float *e, f_cfloat *q, f_int *ldq, f_cfloat *work, f_int *info, f_int vect_len, f_int uplo_len);
void zhbtrd_(char *vect, char *uplo, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_double *d, f_double *e, f_cdouble *q, f_int *ldq, f_cdouble *work, f_int *info, f_int vect_len, f_int uplo_len);

/// Estimates the reciprocal of the condition number of a
/// real symmetric indefinite
/// matrix in packed storage, using the factorization computed
/// by SSPTRF.
void sspcon_(char *uplo, f_int *n, f_float *ap, f_int *ipiv, f_float *anorm, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dspcon_(char *uplo, f_int *n, f_double *ap, f_int *ipiv, f_double *anorm, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void cspcon_(char *uplo, f_int *n, f_cfloat *ap, f_int *ipiv, f_float *anorm, f_float *rcond, f_cfloat *work, f_int *info, f_int uplo_len);
void zspcon_(char *uplo, f_int *n, f_cdouble *ap, f_int *ipiv, f_double *anorm, f_double *rcond, f_cdouble *work, f_int *info, f_int uplo_len);

/// Estimates the reciprocal of the condition number of a
/// complex Hermitian indefinite
/// matrix in packed storage, using the factorization computed
/// by CHPTRF.
void chpcon_(char *uplo, f_int *n, f_cfloat *ap, f_int *ipiv, f_float *anorm, f_float *rcond, f_cfloat *work, f_int *info, f_int uplo_len);
void zhpcon_(char *uplo, f_int *n, f_cdouble *ap, f_int *ipiv, f_double *anorm, f_double *rcond, f_cdouble *work, f_int *info, f_int uplo_len);

/// Reduces a symmetric-definite generalized eigenproblem
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard
/// form,  where A and B are held in packed storage, and B has been
/// factorized by SPPTRF.
void sspgst_(f_int *itype, char *uplo, f_int *n, f_float *ap, f_float *bp, f_int *info, f_int uplo_len);
void dspgst_(f_int *itype, char *uplo, f_int *n, f_double *ap, f_double *bp, f_int *info, f_int uplo_len);

/// Reduces a Hermitian-definite generalized eigenproblem
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard
/// form,  where A and B are held in packed storage, and B has been
/// factorized by CPPTRF.
void chpgst_(f_int *itype, char *uplo, f_int *n, f_cfloat *ap, f_cfloat *bp, f_int *info, f_int uplo_len);
void zhpgst_(f_int *itype, char *uplo, f_int *n, f_cdouble *ap, f_cdouble *bp, f_int *info, f_int uplo_len);

/// Improves the computed solution to a real
/// symmetric indefinite system of linear equations
/// AX=B, where A is held in packed storage, and provides forward
/// and backward error bounds for the solution.
void ssprfs_(char *uplo, f_int *n, f_int *nrhs, f_float *ap, f_float *afp, f_int *ipiv, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dsprfs_(char *uplo, f_int *n, f_int *nrhs, f_double *ap, f_double *afp, f_int *ipiv, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void csprfs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *afp, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zsprfs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *afp, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Improves the computed solution to a complex
/// Hermitian indefinite system of linear equations
/// AX=B, where A is held in packed storage, and provides forward
/// and backward error bounds for the solution.
void chprfs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *afp, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zhprfs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *afp, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Reduces a symmetric matrix in packed storage to real
/// symmetric tridiagonal form by an orthogonal similarity
/// transformation.
void ssptrd_(char *uplo, f_int *n, f_float *ap, f_float *d, f_float *e, f_float *tau, f_int *info, f_int uplo_len);
void dsptrd_(char *uplo, f_int *n, f_double *ap, f_double *d, f_double *e, f_double *tau, f_int *info, f_int uplo_len);

/// Reduces a Hermitian matrix in packed storage to real
/// symmetric tridiagonal form by a unitary similarity
/// transformation.
void chptrd_(char *uplo, f_int *n, f_cfloat *ap, f_float *d, f_float *e, f_cfloat *tau, f_int *info, f_int uplo_len);
void zhptrd_(char *uplo, f_int *n, f_cdouble *ap, f_double *d, f_double *e, f_cdouble *tau, f_int *info, f_int uplo_len);

/// Computes the factorization of a real
/// symmetric-indefinite matrix in packed storage,
/// using the diagonal pivoting method.
void ssptrf_(char *uplo, f_int *n, f_float *ap, f_int *ipiv, f_int *info, f_int uplo_len);
void dsptrf_(char *uplo, f_int *n, f_double *ap, f_int *ipiv, f_int *info, f_int uplo_len);
void csptrf_(char *uplo, f_int *n, f_cfloat *ap, f_int *ipiv, f_int *info, f_int uplo_len);
void zsptrf_(char *uplo, f_int *n, f_cdouble *ap, f_int *ipiv, f_int *info, f_int uplo_len);

/// Computes the factorization of a complex
/// Hermitian-indefinite matrix in packed storage,
/// using the diagonal pivoting method.
void chptrf_(char *uplo, f_int *n, f_cfloat *ap, f_int *ipiv, f_int *info, f_int uplo_len);
void zhptrf_(char *uplo, f_int *n, f_cdouble *ap, f_int *ipiv, f_int *info, f_int uplo_len);

/// Computes the inverse of a real symmetric
/// indefinite matrix in packed storage, using the factorization
/// computed by SSPTRF.
void ssptri_(char *uplo, f_int *n, f_float *ap, f_int *ipiv, f_float *work, f_int *info, f_int uplo_len);
void dsptri_(char *uplo, f_int *n, f_double *ap, f_int *ipiv, f_double *work, f_int *info, f_int uplo_len);
void csptri_(char *uplo, f_int *n, f_cfloat *ap, f_int *ipiv, f_cfloat *work, f_int *info, f_int uplo_len);
void zsptri_(char *uplo, f_int *n, f_cdouble *ap, f_int *ipiv, f_cdouble *work, f_int *info, f_int uplo_len);

/// Computes the inverse of a complex
/// Hermitian indefinite matrix in packed storage, using the factorization
/// computed by CHPTRF.
void chptri_(char *uplo, f_int *n, f_cfloat *ap, f_int *ipiv, f_cfloat *work, f_int *info, f_int uplo_len);
void zhptri_(char *uplo, f_int *n, f_cdouble *ap, f_int *ipiv, f_cdouble *work, f_int *info, f_int uplo_len);

/// Solves a real symmetric
/// indefinite system of linear equations AX=B, where A is held
/// in packed storage, using the factorization computed
/// by SSPTRF.
void ssptrs_(char *uplo, f_int *n, f_int *nrhs, f_float *ap, f_int *ipiv, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dsptrs_(char *uplo, f_int *n, f_int *nrhs, f_double *ap, f_int *ipiv, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void csptrs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zsptrs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Solves a complex Hermitian
/// indefinite system of linear equations AX=B, where A is held
/// in packed storage, using the factorization computed
/// by CHPTRF.
void chptrs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *ap, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zhptrs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *ap, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Computes selected eigenvalues of a real symmetric tridiagonal
/// matrix by bisection.
void sstebz_(char *range, char *order, f_int *n, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_float *d, f_float *e, f_int *m, f_int *nsplit, f_float *w, f_int *iblock, f_int *isplit, f_float *work, f_int *iwork, f_int *info, f_int range_len, f_int order_len);
void dstebz_(char *range, char *order, f_int *n, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_double *d, f_double *e, f_int *m, f_int *nsplit, f_double *w, f_int *iblock, f_int *isplit, f_double *work, f_int *iwork, f_int *info, f_int range_len, f_int order_len);

/// Computes all eigenvalues and, optionally, eigenvectors of a
/// symmetric tridiagonal matrix using the divide and conquer algorithm.
void sstedc_(char *compz, f_int *n, f_float *d, f_float *e, f_float *z, f_int *ldz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int compz_len);
void dstedc_(char *compz, f_int *n, f_double *d, f_double *e, f_double *z, f_int *ldz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int compz_len);
void cstedc_(char *compz, f_int *n, f_float *d, f_float *e, f_cfloat *z, f_int *ldz, f_cfloat *work, f_int *lwork, f_float *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int compz_len);
void zstedc_(char *compz, f_int *n, f_double *d, f_double *e, f_cdouble *z, f_int *ldz, f_cdouble *work, f_int *lwork, f_double *rwork, f_int *lrwork, f_int *iwork, f_int *liwork, f_int *info, f_int compz_len);

/// Computes selected eigenvalues and, optionally, eigenvectors of a
/// symmetric tridiagonal matrix.  The eigenvalues are computed by the
/// dqds algorithm, while eigenvectors are computed from various "good"
/// LDL^T representations (also known as Relatively Robust Representations.)
void sstegr_(char *jobz, char *range, f_int *n, f_float *d, f_float *e, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_float *z, f_int *ldz, f_int *isuppz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len);
void dstegr_(char *jobz, char *range, f_int *n, f_double *d, f_double *e, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_double *z, f_int *ldz, f_int *isuppz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len);
void cstegr_(char *jobz, char *range, f_int *n, f_float *d, f_float *e, f_float *vl, f_float *vu, f_int *il, f_int *iu, f_float *abstol, f_int *m, f_float *w, f_cfloat *z, f_int *ldz, f_int *isuppz, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len);
void zstegr_(char *jobz, char *range, f_int *n, f_double *d, f_double *e, f_double *vl, f_double *vu, f_int *il, f_int *iu, f_double *abstol, f_int *m, f_double *w, f_cdouble *z, f_int *ldz, f_int *isuppz, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int jobz_len, f_int range_len);

/// Computes selected eigenvectors of a real symmetric tridiagonal
/// matrix by inverse iteration.
void sstein_(f_int *n, f_float *d, f_float *e, f_int *m, f_float *w, f_int *iblock, f_int *isplit, f_float *z, f_int *ldz, f_float *work, f_int *iwork, f_int *ifail, f_int *info);
void dstein_(f_int *n, f_double *d, f_double *e, f_int *m, f_double *w, f_int *iblock, f_int *isplit, f_double *z, f_int *ldz, f_double *work, f_int *iwork, f_int *ifail, f_int *info);
void cstein_(f_int *n, f_float *d, f_float *e, f_int *m, f_float *w, f_int *iblock, f_int *isplit, f_cfloat *z, f_int *ldz, f_float *work, f_int *iwork, f_int *ifail, f_int *info);
void zstein_(f_int *n, f_double *d, f_double *e, f_int *m, f_double *w, f_int *iblock, f_int *isplit, f_cdouble *z, f_int *ldz, f_double *work, f_int *iwork, f_int *ifail, f_int *info);

/// Computes all eigenvalues and eigenvectors of a real symmetric
/// tridiagonal matrix, using the implicit QL or QR algorithm.
void ssteqr_(char *compz, f_int *n, f_float *d, f_float *e, f_float *z, f_int *ldz, f_float *work, f_int *info, f_int compz_len);
void dsteqr_(char *compz, f_int *n, f_double *d, f_double *e, f_double *z, f_int *ldz, f_double *work, f_int *info, f_int compz_len);
void csteqr_(char *compz, f_int *n, f_float *d, f_float *e, f_cfloat *z, f_int *ldz, f_float *work, f_int *info, f_int compz_len);
void zsteqr_(char *compz, f_int *n, f_double *d, f_double *e, f_cdouble *z, f_int *ldz, f_double *work, f_int *info, f_int compz_len);

/// Computes all eigenvalues of a real symmetric tridiagonal matrix,
/// using a root-free variant of the QL or QR algorithm.
void ssterf_(f_int *n, f_float *d, f_float *e, f_int *info);
void dsterf_(f_int *n, f_double *d, f_double *e, f_int *info);

/// Estimates the reciprocal of the condition number of a
/// real symmetric indefinite matrix,
/// using the factorization computed by SSYTRF.
void ssycon_(char *uplo, f_int *n, f_float *a, f_int *lda, f_int *ipiv, f_float *anorm, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dsycon_(char *uplo, f_int *n, f_double *a, f_int *lda, f_int *ipiv, f_double *anorm, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void csycon_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_int *ipiv, f_float *anorm, f_float *rcond, f_cfloat *work, f_int *info, f_int uplo_len);
void zsycon_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_int *ipiv, f_double *anorm, f_double *rcond, f_cdouble *work, f_int *info, f_int uplo_len);

/// Estimates the reciprocal of the condition number of a
/// complex Hermitian indefinite matrix,
/// using the factorization computed by CHETRF.
void checon_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_int *ipiv, f_float *anorm, f_float *rcond, f_cfloat *work, f_int *info, f_int uplo_len);
void zhecon_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_int *ipiv, f_double *anorm, f_double *rcond, f_cdouble *work, f_int *info, f_int uplo_len);

/// Reduces a symmetric-definite generalized eigenproblem
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard
/// form, where B has been factorized by SPOTRF.
void ssygst_(f_int *itype, char *uplo, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dsygst_(f_int *itype, char *uplo, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Reduces a Hermitian-definite generalized eigenproblem
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard
/// form, where B has been factorized by CPOTRF.
void chegst_(f_int *itype, char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zhegst_(f_int *itype, char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Improves the computed solution to a real
/// symmetric indefinite system of linear equations
/// AX=B, and provides forward and backward error bounds for the
/// solution.
void ssyrfs_(char *uplo, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *af, f_int *ldaf, f_int *ipiv, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int uplo_len);
void dsyrfs_(char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *af, f_int *ldaf, f_int *ipiv, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int uplo_len);
void csyrfs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *af, f_int *ldaf, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zsyrfs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *af, f_int *ldaf, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Improves the computed solution to a complex
/// Hermitian indefinite system of linear equations
/// AX=B, and provides forward and backward error bounds for the
/// solution.
void cherfs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *af, f_int *ldaf, f_int *ipiv, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len);
void zherfs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *af, f_int *ldaf, f_int *ipiv, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len);

/// Reduces a symmetric matrix to real symmetric tridiagonal
/// form by an orthogonal similarity transformation.
void ssytrd_(char *uplo, f_int *n, f_float *a, f_int *lda, f_float *d, f_float *e, f_float *tau, f_float *work, f_int *lwork, f_int *info, f_int uplo_len);
void dsytrd_(char *uplo, f_int *n, f_double *a, f_int *lda, f_double *d, f_double *e, f_double *tau, f_double *work, f_int *lwork, f_int *info, f_int uplo_len);

/// Reduces a Hermitian matrix to real symmetric tridiagonal
/// form by an orthogonal/unitary similarity transformation.
void chetrd_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_float *d, f_float *e, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info, f_int uplo_len);
void zhetrd_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_double *d, f_double *e, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info, f_int uplo_len);

/// Computes the factorization of a real symmetric-indefinite matrix,
/// using the diagonal pivoting method.
void ssytrf_(char *uplo, f_int *n, f_float *a, f_int *lda, f_int *ipiv, f_float *work, f_int *lwork, f_int *info, f_int uplo_len);
void dsytrf_(char *uplo, f_int *n, f_double *a, f_int *lda, f_int *ipiv, f_double *work, f_int *lwork, f_int *info, f_int uplo_len);
void csytrf_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *work, f_int *lwork, f_int *info, f_int uplo_len);
void zsytrf_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *work, f_int *lwork, f_int *info, f_int uplo_len);

/// Computes the factorization of a complex Hermitian-indefinite matrix,
/// using the diagonal pivoting method.
void chetrf_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *work, f_int *lwork, f_int *info, f_int uplo_len);
void zhetrf_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *work, f_int *lwork, f_int *info, f_int uplo_len);

/// Computes the inverse of a real symmetric indefinite matrix,
/// using the factorization computed by SSYTRF.
void ssytri_(char *uplo, f_int *n, f_float *a, f_int *lda, f_int *ipiv, f_float *work, f_int *info, f_int uplo_len);
void dsytri_(char *uplo, f_int *n, f_double *a, f_int *lda, f_int *ipiv, f_double *work, f_int *info, f_int uplo_len);
void csytri_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *work, f_int *info, f_int uplo_len);
void zsytri_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *work, f_int *info, f_int uplo_len);

/// Computes the inverse of a complex Hermitian indefinite matrix,
/// using the factorization computed by CHETRF.
void chetri_(char *uplo, f_int *n, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *work, f_int *info, f_int uplo_len);
void zhetri_(char *uplo, f_int *n, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *work, f_int *info, f_int uplo_len);

/// Solves a real symmetric indefinite system of linear equations AX=B,
/// using the factorization computed by SSPTRF.
void ssytrs_(char *uplo, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_int *ipiv, f_float *b, f_int *ldb, f_int *info, f_int uplo_len);
void dsytrs_(char *uplo, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_int *ipiv, f_double *b, f_int *ldb, f_int *info, f_int uplo_len);
void csytrs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zsytrs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Solves a complex Hermitian indefinite system of linear equations AX=B,
/// using the factorization computed by CHPTRF.
void chetrs_(char *uplo, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_int *ipiv, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len);
void zhetrs_(char *uplo, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_int *ipiv, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len);

/// Estimates the reciprocal of the condition number of a triangular
/// band matrix, in either the 1-norm or the infinity-norm.
void stbcon_(char *norm, char *uplo, char *diag, f_int *n, f_int *kd, f_float *ab, f_int *ldab, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void dtbcon_(char *norm, char *uplo, char *diag, f_int *n, f_int *kd, f_double *ab, f_int *ldab, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void ctbcon_(char *norm, char *uplo, char *diag, f_int *n, f_int *kd, f_cfloat *ab, f_int *ldab, f_float *rcond, f_cfloat *work, f_float *rwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void ztbcon_(char *norm, char *uplo, char *diag, f_int *n, f_int *kd, f_cdouble *ab, f_int *ldab, f_double *rcond, f_cdouble *work, f_double *rwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);

/// Provides forward and backward error bounds for the solution
/// of a triangular banded system of linear equations AX=B,
/// A**T X=B or A**H X=B.
void stbrfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *kd, f_int *nrhs, f_float *ab, f_int *ldab, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtbrfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *kd, f_int *nrhs, f_double *ab, f_int *ldab, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctbrfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *kd, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztbrfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *kd, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);

/// Solves a triangular banded system of linear equations AX=B,
/// A**T X=B or A**H X=B.
void stbtrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *kd, f_int *nrhs, f_float *ab, f_int *ldab, f_float *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtbtrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *kd, f_int *nrhs, f_double *ab, f_int *ldab, f_double *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctbtrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *kd, f_int *nrhs, f_cfloat *ab, f_int *ldab, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztbtrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *kd, f_int *nrhs, f_cdouble *ab, f_int *ldab, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);

/// Computes some or all of the right and/or left generalized eigenvectors
/// of a pair of upper triangular matrices.
void stgevc_(char *side, char *howmny, f_int *select, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_int *mm, f_int *m, f_float *work, f_int *info, f_int side_len, f_int howmny_len);
void dtgevc_(char *side, char *howmny, f_int *select, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_int *mm, f_int *m, f_double *work, f_int *info, f_int side_len, f_int howmny_len);
void ctgevc_(char *side, char *howmny, f_int *select, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_int *mm, f_int *m, f_cfloat *work, f_float *rwork, f_int *info, f_int side_len, f_int howmny_len);
void ztgevc_(char *side, char *howmny, f_int *select, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_int *mm, f_int *m, f_cdouble *work, f_double *rwork, f_int *info, f_int side_len, f_int howmny_len);

/// Reorders the generalized real Schur decomposition of a real
/// matrix pair (A,B) using an orthogonal equivalence transformation
/// so that the diagonal block of (A,B) with row index IFST is moved
/// to row ILST.
void stgexc_(f_int *wantq, f_int *wantz, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *q, f_int *ldq, f_float *z, f_int *ldz, f_int *ifst, f_int *ilst, f_float *work, f_int *lwork, f_int *info);
void dtgexc_(f_int *wantq, f_int *wantz, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *q, f_int *ldq, f_double *z, f_int *ldz, f_int *ifst, f_int *ilst, f_double *work, f_int *lwork, f_int *info);
void ctgexc_(f_int *wantq, f_int *wantz, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *q, f_int *ldq, f_cfloat *z, f_int *ldz, f_int *ifst, f_int *ilst, f_int *info);
void ztgexc_(f_int *wantq, f_int *wantz, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *q, f_int *ldq, f_cdouble *z, f_int *ldz, f_int *ifst, f_int *ilst, f_int *info);

/// Reorders the generalized real Schur decomposition of a real
/// matrix pair (A, B) so that a selected cluster of eigenvalues
/// appears in the leading diagonal blocks of the upper quasi-triangular
/// matrix A and the upper triangular B.
void stgsen_(f_int *ijob, f_int *wantq, f_int *wantz, f_int *select, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *q, f_int *ldq, f_float *z, f_int *ldz, f_int *m, f_float *pl, f_float *pr, f_float *dif, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info);
void dtgsen_(f_int *ijob, f_int *wantq, f_int *wantz, f_int *select, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *q, f_int *ldq, f_double *z, f_int *ldz, f_int *m, f_double *pl, f_double *pr, f_double *dif, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info);
void ctgsen_(f_int *ijob, f_int *wantq, f_int *wantz, f_int *select, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *q, f_int *ldq, f_cfloat *z, f_int *ldz, f_int *m, f_float *pl, f_float *pr, f_float *dif, f_cfloat *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info);
void ztgsen_(f_int *ijob, f_int *wantq, f_int *wantz, f_int *select, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *q, f_int *ldq, f_cdouble *z, f_int *ldz, f_int *m, f_double *pl, f_double *pr, f_double *dif, f_cdouble *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info);

/// Computes the generalized singular value decomposition of two real
/// upper triangular (or trapezoidal) matrices as output by SGGSVP.
void stgsja_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *p, f_int *n, f_int *k, f_int *l, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *tola, f_float *tolb, f_float *alphav, f_float *betav, f_float *u, f_int *ldu, f_float *v, f_int *ldv, f_float *q, f_int *ldq, f_float *work, f_int *ncycle, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void dtgsja_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *p, f_int *n, f_int *k, f_int *l, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *tola, f_double *tolb, f_double *alphav, f_double *betav, f_double *u, f_int *ldu, f_double *v, f_int *ldv, f_double *q, f_int *ldq, f_double *work, f_int *ncycle, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void ctgsja_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *p, f_int *n, f_int *k, f_int *l, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_float *tola, f_float *tolb, f_float *alphav, f_float *betav, f_cfloat *u, f_int *ldu, f_cfloat *v, f_int *ldv, f_cfloat *q, f_int *ldq, f_cfloat *work, f_int *ncycle, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);
void ztgsja_(char *jobu, char *jobv, char *jobq, f_int *m, f_int *p, f_int *n, f_int *k, f_int *l, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_double *tola, f_double *tolb, f_double *alphav, f_double *betav, f_cdouble *u, f_int *ldu, f_cdouble *v, f_int *ldv, f_cdouble *q, f_int *ldq, f_cdouble *work, f_int *ncycle, f_int *info, f_int jobu_len, f_int jobv_len, f_int jobq_len);

/// Estimates reciprocal condition numbers for specified
/// eigenvalues and/or eigenvectors of a matrix pair (A, B) in
/// generalized real Schur canonical form, as returned by SGGES.
void stgsna_(char *job, char *howmny, f_int *select, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_float *s, f_float *dif, f_int *mm, f_int *m, f_float *work, f_int *lwork, f_int *iwork, f_int *info, f_int job_len, f_int howmny_len);
void dtgsna_(char *job, char *howmny, f_int *select, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_double *s, f_double *dif, f_int *mm, f_int *m, f_double *work, f_int *lwork, f_int *iwork, f_int *info, f_int job_len, f_int howmny_len);
void ctgsna_(char *job, char *howmny, f_int *select, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_float *s, f_float *dif, f_int *mm, f_int *m, f_cfloat *work, f_int *lwork, f_int *iwork, f_int *info, f_int job_len, f_int howmny_len);
void ztgsna_(char *job, char *howmny, f_int *select, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_double *s, f_double *dif, f_int *mm, f_int *m, f_cdouble *work, f_int *lwork, f_int *iwork, f_int *info, f_int job_len, f_int howmny_len);

/// Solves the generalized Sylvester equation.
void stgsyl_(char *trans, f_int *ijob, f_int *m, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *c, f_int *ldc, f_float *d, f_int *ldd, f_float *e, f_int *lde, f_float *f, f_int *ldf, f_float *scale, f_float *dif, f_float *work, f_int *lwork, f_int *iwork, f_int *info, f_int trans_len);
void dtgsyl_(char *trans, f_int *ijob, f_int *m, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *c, f_int *ldc, f_double *d, f_int *ldd, f_double *e, f_int *lde, f_double *f, f_int *ldf, f_double *scale, f_double *dif, f_double *work, f_int *lwork, f_int *iwork, f_int *info, f_int trans_len);
void ctgsyl_(char *trans, f_int *ijob, f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *c, f_int *ldc, f_cfloat *d, f_int *ldd, f_cfloat *e, f_int *lde, f_cfloat *f, f_int *ldf, f_float *scale, f_float *dif, f_cfloat *work, f_int *lwork, f_int *iwork, f_int *info, f_int trans_len);
void ztgsyl_(char *trans, f_int *ijob, f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *c, f_int *ldc, f_cdouble *d, f_int *ldd, f_cdouble *e, f_int *lde, f_cdouble *f, f_int *ldf, f_double *scale, f_double *dif, f_cdouble *work, f_int *lwork, f_int *iwork, f_int *info, f_int trans_len);

/// Estimates the reciprocal of the condition number of a triangular
/// matrix in packed storage, in either the 1-norm or the infinity-norm.
void stpcon_(char *norm, char *uplo, char *diag, f_int *n, f_float *ap, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void dtpcon_(char *norm, char *uplo, char *diag, f_int *n, f_double *ap, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void ctpcon_(char *norm, char *uplo, char *diag, f_int *n, f_cfloat *ap, f_float *rcond, f_cfloat *work, f_float *rwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void ztpcon_(char *norm, char *uplo, char *diag, f_int *n, f_cdouble *ap, f_double *rcond, f_cdouble *work, f_double *rwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);

/// Provides forward and backward error bounds for the solution
/// of a triangular system of linear equations AX=B, A**T X=B or
/// A**H X=B, where A is held in packed storage.
void stprfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_float *ap, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtprfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_double *ap, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctprfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztprfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);

///  Computes the inverse of a triangular matrix in packed storage.
void stptri_(char *uplo, char *diag, f_int *n, f_float *ap, f_int *info, f_int uplo_len, f_int diag_len);
void dtptri_(char *uplo, char *diag, f_int *n, f_double *ap, f_int *info, f_int uplo_len, f_int diag_len);
void ctptri_(char *uplo, char *diag, f_int *n, f_cfloat *ap, f_int *info, f_int uplo_len, f_int diag_len);
void ztptri_(char *uplo, char *diag, f_int *n, f_cdouble *ap, f_int *info, f_int uplo_len, f_int diag_len);

/// Solves a triangular system of linear equations AX=B,
/// A**T X=B or A**H X=B, where A is held in packed storage.
void stptrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_float *ap, f_float *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtptrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_double *ap, f_double *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctptrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_cfloat *ap, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztptrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_cdouble *ap, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);

/// Estimates the reciprocal of the condition number of a triangular
/// matrix, in either the 1-norm or the infinity-norm.
void strcon_(char *norm, char *uplo, char *diag, f_int *n, f_float *a, f_int *lda, f_float *rcond, f_float *work, f_int *iwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void dtrcon_(char *norm, char *uplo, char *diag, f_int *n, f_double *a, f_int *lda, f_double *rcond, f_double *work, f_int *iwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void ctrcon_(char *norm, char *uplo, char *diag, f_int *n, f_cfloat *a, f_int *lda, f_float *rcond, f_cfloat *work, f_float *rwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);
void ztrcon_(char *norm, char *uplo, char *diag, f_int *n, f_cdouble *a, f_int *lda, f_double *rcond, f_cdouble *work, f_double *rwork, f_int *info, f_int norm_len, f_int uplo_len, f_int diag_len);

/// Computes some or all of the right and/or left eigenvectors of
/// an upper quasi-triangular matrix.
void strevc_(char *side, char *howmny, f_int *select, f_int *n, f_float *t, f_int *ldt, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_int *mm, f_int *m, f_float *work, f_int *info, f_int side_len, f_int howmny_len);
void dtrevc_(char *side, char *howmny, f_int *select, f_int *n, f_double *t, f_int *ldt, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_int *mm, f_int *m, f_double *work, f_int *info, f_int side_len, f_int howmny_len);
void ctrevc_(char *side, char *howmny, f_int *select, f_int *n, f_cfloat *t, f_int *ldt, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_int *mm, f_int *m, f_cfloat *work, f_float *rwork, f_int *info, f_int side_len, f_int howmny_len);
void ztrevc_(char *side, char *howmny, f_int *select, f_int *n, f_cdouble *t, f_int *ldt, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_int *mm, f_int *m, f_cdouble *work, f_double *rwork, f_int *info, f_int side_len, f_int howmny_len);

/// Reorders the Schur factorization of a matrix by an orthogonal
/// similarity transformation.
void strexc_(char *compq, f_int *n, f_float *t, f_int *ldt, f_float *q, f_int *ldq, f_int *ifst, f_int *ilst, f_float *work, f_int *info, f_int compq_len);
void dtrexc_(char *compq, f_int *n, f_double *t, f_int *ldt, f_double *q, f_int *ldq, f_int *ifst, f_int *ilst, f_double *work, f_int *info, f_int compq_len);
void ctrexc_(char *compq, f_int *n, f_cfloat *t, f_int *ldt, f_cfloat *q, f_int *ldq, f_int *ifst, f_int *ilst, f_int *info, f_int compq_len);
void ztrexc_(char *compq, f_int *n, f_cdouble *t, f_int *ldt, f_cdouble *q, f_int *ldq, f_int *ifst, f_int *ilst, f_int *info, f_int compq_len);

/// Provides forward and backward error bounds for the solution
/// of a triangular system of linear equations A X=B, A**T X=B or
/// A**H X=B.
void strrfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *x, f_int *ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtrrfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *x, f_int *ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctrrfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *x, f_int *ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztrrfs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *x, f_int *ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);

/// Reorders the Schur factorization of a matrix in order to find
/// an orthonormal basis of a right invariant subspace corresponding
/// to selected eigenvalues, and returns reciprocal condition numbers
/// (sensitivities) of the average of the cluster of eigenvalues
/// and of the invariant subspace.
void strsen_(char *job, char *compq, f_int *select, f_int *n, f_float *t, f_int *ldt, f_float *q, f_int *ldq, f_float *wr, f_float *wi, f_int *m, f_float *s, f_float *sep, f_float *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int job_len, f_int compq_len);
void dtrsen_(char *job, char *compq, f_int *select, f_int *n, f_double *t, f_int *ldt, f_double *q, f_int *ldq, f_double *wr, f_double *wi, f_int *m, f_double *s, f_double *sep, f_double *work, f_int *lwork, f_int *iwork, f_int *liwork, f_int *info, f_int job_len, f_int compq_len);
void ctrsen_(char *job, char *compq, f_int *select, f_int *n, f_cfloat *t, f_int *ldt, f_cfloat *q, f_int *ldq, f_cfloat *w, f_int *m, f_float *s, f_float *sep, f_cfloat *work, f_int *lwork, f_int *info, f_int job_len, f_int compq_len);
void ztrsen_(char *job, char *compq, f_int *select, f_int *n, f_cdouble *t, f_int *ldt, f_cdouble *q, f_int *ldq, f_cdouble *w, f_int *m, f_double *s, f_double *sep, f_cdouble *work, f_int *lwork, f_int *info, f_int job_len, f_int compq_len);

/// Estimates the reciprocal condition numbers (sensitivities)
/// of selected eigenvalues and eigenvectors of an upper
/// quasi-triangular matrix.
void strsna_(char *job, char *howmny, f_int *select, f_int *n, f_float *t, f_int *ldt, f_float *vl, f_int *ldvl, f_float *vr, f_int *ldvr, f_float *s, f_float *sep, f_int *mm, f_int *m, f_float *work, f_int *ldwork, f_int *iwork, f_int *info, f_int job_len, f_int howmny_len);
void dtrsna_(char *job, char *howmny, f_int *select, f_int *n, f_double *t, f_int *ldt, f_double *vl, f_int *ldvl, f_double *vr, f_int *ldvr, f_double *s, f_double *sep, f_int *mm, f_int *m, f_double *work, f_int *ldwork, f_int *iwork, f_int *info, f_int job_len, f_int howmny_len);
void ctrsna_(char *job, char *howmny, f_int *select, f_int *n, f_cfloat *t, f_int *ldt, f_cfloat *vl, f_int *ldvl, f_cfloat *vr, f_int *ldvr, f_float *s, f_float *sep, f_int *mm, f_int *m, f_cfloat *work, f_int *ldwork, f_float *rwork, f_int *info, f_int job_len, f_int howmny_len);
void ztrsna_(char *job, char *howmny, f_int *select, f_int *n, f_cdouble *t, f_int *ldt, f_cdouble *vl, f_int *ldvl, f_cdouble *vr, f_int *ldvr, f_double *s, f_double *sep, f_int *mm, f_int *m, f_cdouble *work, f_int *ldwork, f_double *rwork, f_int *info, f_int job_len, f_int howmny_len);

/// Solves the Sylvester matrix equation A X +/- X B=C where A
/// and B are upper quasi-triangular, and may be transposed.
void strsyl_(char *trana, char *tranb, f_int *isgn, f_int *m, f_int *n, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_float *c, f_int *ldc, f_float *scale, f_int *info, f_int trana_len, f_int tranb_len);
void dtrsyl_(char *trana, char *tranb, f_int *isgn, f_int *m, f_int *n, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_double *c, f_int *ldc, f_double *scale, f_int *info, f_int trana_len, f_int tranb_len);
void ctrsyl_(char *trana, char *tranb, f_int *isgn, f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_cfloat *c, f_int *ldc, f_float *scale, f_int *info, f_int trana_len, f_int tranb_len);
void ztrsyl_(char *trana, char *tranb, f_int *isgn, f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_cdouble *c, f_int *ldc, f_double *scale, f_int *info, f_int trana_len, f_int tranb_len);

/// Computes the inverse of a triangular matrix.
void strtri_(char *uplo, char *diag, f_int *n, f_float *a, f_int *lda, f_int *info, f_int uplo_len, f_int diag_len);
void dtrtri_(char *uplo, char *diag, f_int *n, f_double *a, f_int *lda, f_int *info, f_int uplo_len, f_int diag_len);
void ctrtri_(char *uplo, char *diag, f_int *n, f_cfloat *a, f_int *lda, f_int *info, f_int uplo_len, f_int diag_len);
void ztrtri_(char *uplo, char *diag, f_int *n, f_cdouble *a, f_int *lda, f_int *info, f_int uplo_len, f_int diag_len);

/// Solves a triangular system of linear equations AX=B,
/// A**T X=B or A**H X=B.
void strtrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_float *a, f_int *lda, f_float *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void dtrtrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_double *a, f_int *lda, f_double *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ctrtrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_cfloat *a, f_int *lda, f_cfloat *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);
void ztrtrs_(char *uplo, char *trans, char *diag, f_int *n, f_int *nrhs, f_cdouble *a, f_int *lda, f_cdouble *b, f_int *ldb, f_int *info, f_int uplo_len, f_int trans_len, f_int diag_len);

/// Computes an RQ factorization of an upper trapezoidal matrix.
void stzrqf_(f_int *m, f_int *n, f_float *a, f_int *lda, f_float *tau, f_int *info);
void dtzrqf_(f_int *m, f_int *n, f_double *a, f_int *lda, f_double *tau, f_int *info);
void ctzrqf_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *tau, f_int *info);
void ztzrqf_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *tau, f_int *info);

/// Computes an RZ factorization of an upper trapezoidal matrix
/// (blocked version of STZRQF).
void stzrzf_(f_int *m, f_int *n, f_float *a, f_int *lda, f_float *tau, f_float *work, f_int *lwork, f_int *info);
void dtzrzf_(f_int *m, f_int *n, f_double *a, f_int *lda, f_double *tau, f_double *work, f_int *lwork, f_int *info);
void ctzrzf_(f_int *m, f_int *n, f_cfloat *a, f_int *lda, f_cfloat *tau, f_cfloat *work, f_int *lwork, f_int *info);
void ztzrzf_(f_int *m, f_int *n, f_cdouble *a, f_int *lda, f_cdouble *tau, f_cdouble *work, f_int *lwork, f_int *info);


/// Multiplies a general matrix by the unitary
/// transformation matrix from a reduction to tridiagonal form
/// determined by CHPTRD.
void cupmtr_(char *side, char *uplo, char *trans, f_int *m, f_int *n, f_cfloat *ap, f_cfloat *tau, f_cfloat *c, f_int *ldc, f_cfloat *work, f_int *info, f_int side_len, f_int uplo_len, f_int trans_len);
void zupmtr_(char *side, char *uplo, char *trans, f_int *m, f_int *n, f_cdouble *ap, f_cdouble *tau, f_cdouble *c, f_int *ldc, f_cdouble *work, f_int *info, f_int side_len, f_int uplo_len, f_int trans_len);


//------------------------------------
//     ----- MISC routines -----
//------------------------------------

f_int ilaenv_(f_int *ispec, char *name, char *opts, f_int *n1, f_int *n2, f_int *n3, f_int *n4, f_int len_name, f_int len_opts);
void ilaenvset_(f_int *ispec, char *name, char *opts, f_int *n1, f_int *n2, f_int *n3, f_int *n4, f_int *nvalue, f_int *info, f_int len_name, f_int len_opts);

///
f_float slamch_(char *cmach, f_int cmach_len);
f_double dlamch_(char *cmach, f_int cmach_len);

///
lapack_float_ret_t second_();
f_double dsecnd_();


