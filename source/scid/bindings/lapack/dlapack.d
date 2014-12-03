/** LAPACK bindings for D.

    Authors:    William V. Baxter III (with slight modifications by
                Lars Tandle Kyllingstad).
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.bindings.lapack.dlapack;


public  import scid.bindings.blas.types;
import scid.bindings.lapack.lapack;
import scid.core.fortran;




/* LAPACK routines */

//--------------------------------------------------------
// ---- SIMPLE and DIVIDE AND CONQUER DRIVER routines ----
//---------------------------------------------------------

/// Solves a general system of linear equations AX=B.
void gesv(f_int n, f_int nrhs, f_float *a, f_int lda, f_int *ipiv, f_float *b, f_int ldb, ref f_int info) {
    sgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
}
void gesv(f_int n, f_int nrhs, f_double *a, f_int lda, f_int *ipiv, f_double *b, f_int ldb, ref f_int info) {
    dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
}
void gesv(f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    cgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
}
void gesv(f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);
}

/// Solves a general banded system of linear equations AX=B.
void gbsv(f_int n, f_int kl, f_int ku, f_int nrhs, f_float *ab, f_int ldab, f_int *ipiv, f_float *b, f_int ldb, ref f_int info) {
    sgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
}
void gbsv(f_int n, f_int kl, f_int ku, f_int nrhs, f_double *ab, f_int ldab, f_int *ipiv, f_double *b, f_int ldb, ref f_int info) {
    dgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
}
void gbsv(f_int n, f_int kl, f_int ku, f_int nrhs, f_cfloat *ab, f_int ldab, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    cgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
}
void gbsv(f_int n, f_int kl, f_int ku, f_int nrhs, f_cdouble *ab, f_int ldab, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info);
}

/// Solves a general tridiagonal system of linear equations AX=B.
void gtsv(f_int n, f_int nrhs, f_float *dl, f_float *d, f_float *du, f_float *b, f_int ldb, ref f_int info) {
    sgtsv_(&n, &nrhs, dl, d, du, b, &ldb, &info);
}
void gtsv(f_int n, f_int nrhs, f_double *dl, f_double *d, f_double *du, f_double *b, f_int ldb, ref f_int info) {
    dgtsv_(&n, &nrhs, dl, d, du, b, &ldb, &info);
}
void gtsv(f_int n, f_int nrhs, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *b, f_int ldb, ref f_int info) {
    cgtsv_(&n, &nrhs, dl, d, du, b, &ldb, &info);
}
void gtsv(f_int n, f_int nrhs, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *b, f_int ldb, ref f_int info) {
    zgtsv_(&n, &nrhs, dl, d, du, b, &ldb, &info);
}

/// Solves a symmetric positive definite system of linear
/// equations AX=B.
void posv(char uplo, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, ref f_int info) {
    sposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
}
void posv(char uplo, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, ref f_int info) {
    dposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
}
void posv(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, ref f_int info) {
    cposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
}
void posv(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, ref f_int info) {
    zposv_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
}

/// Solves a symmetric positive definite system of linear
/// equations AX=B, where A is held in packed storage.
void ppsv(char uplo, f_int n, f_int nrhs, f_float *ap, f_float *b, f_int ldb, ref f_int info) {
    sppsv_(&uplo, &n, &nrhs, ap, b, &ldb, &info, 1);
}
void ppsv(char uplo, f_int n, f_int nrhs, f_double *ap, f_double *b, f_int ldb, ref f_int info) {
    dppsv_(&uplo, &n, &nrhs, ap, b, &ldb, &info, 1);
}
void ppsv(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *b, f_int ldb, ref f_int info) {
    cppsv_(&uplo, &n, &nrhs, ap, b, &ldb, &info, 1);
}
void ppsv(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *b, f_int ldb, ref f_int info) {
    zppsv_(&uplo, &n, &nrhs, ap, b, &ldb, &info, 1);
}

/// Solves a symmetric positive definite banded system
/// of linear equations AX=B.
void pbsv(char uplo, f_int n, f_int kd, f_int nrhs, f_float *ab, f_int ldab, f_float *b, f_int ldb, ref f_int info) {
    spbsv_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1);
}
void pbsv(char uplo, f_int n, f_int kd, f_int nrhs, f_double *ab, f_int ldab, f_double *b, f_int ldb, ref f_int info) {
    dpbsv_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1);
}
void pbsv(char uplo, f_int n, f_int kd, f_int nrhs, f_cfloat *ab, f_int ldab, f_cfloat *b, f_int ldb, ref f_int info) {
    cpbsv_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1);
}
void pbsv(char uplo, f_int n, f_int kd, f_int nrhs, f_cdouble *ab, f_int ldab, f_cdouble *b, f_int ldb, ref f_int info) {
    zpbsv_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1);
}

/// Solves a symmetric positive definite tridiagonal system
/// of linear equations AX=B.
void ptsv(f_int n, f_int nrhs, f_float *d, f_float *e, f_float *b, f_int ldb, ref f_int info) {
    sptsv_(&n, &nrhs, d, e, b, &ldb, &info);
}
void ptsv(f_int n, f_int nrhs, f_double *d, f_double *e, f_double *b, f_int ldb, ref f_int info) {
    dptsv_(&n, &nrhs, d, e, b, &ldb, &info);
}
void ptsv(f_int n, f_int nrhs, f_float *d, f_cfloat *e, f_cfloat *b, f_int ldb, ref f_int info) {
    cptsv_(&n, &nrhs, d, e, b, &ldb, &info);
}
void ptsv(f_int n, f_int nrhs, f_double *d, f_cdouble *e, f_cdouble *b, f_int ldb, ref f_int info) {
    zptsv_(&n, &nrhs, d, e, b, &ldb, &info);
}


/// Solves a real symmetric indefinite system of linear equations AX=B.
void sysv(char uplo, f_int n, f_int nrhs, f_float *a, f_int lda, f_int *ipiv, f_float *b, f_int ldb, f_float *work, f_int lwork, ref f_int info) {
    ssysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info, 1);
}
void sysv(char uplo, f_int n, f_int nrhs, f_double *a, f_int lda, f_int *ipiv, f_double *b, f_int ldb, f_double *work, f_int lwork, ref f_int info) {
    dsysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info, 1);
}
void sysv(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *work, f_int lwork, ref f_int info) {
    csysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info, 1);
}
void sysv(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *work, f_int lwork, ref f_int info) {
    zsysv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info, 1);
}

/// Solves a complex Hermitian indefinite system of linear equations AX=B.
void hesv(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *work, f_int lwork, ref f_int info) {
    chesv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info, 1);
}
void hesv(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *work, f_int lwork, ref f_int info) {
    zhesv_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info, 1);
}

/// Solves a real symmetric indefinite system of linear equations AX=B,
/// where A is held in packed storage.
void spsv(char uplo, f_int n, f_int nrhs, f_float *ap, f_int *ipiv, f_float *b, f_int ldb, ref f_int info) {
    sspsv_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}
void spsv(char uplo, f_int n, f_int nrhs, f_double *ap, f_int *ipiv, f_double *b, f_int ldb, ref f_int info) {
    dspsv_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}
void spsv(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    cspsv_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}
void spsv(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zspsv_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}

/// Solves a complex Hermitian indefinite system of linear equations AX=B,
/// where A is held in packed storage.
void hpsv(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    chpsv_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}
void hpsv(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zhpsv_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}

/// Computes the least squares solution to an over-determined system
/// of linear equations, A X=B or A**H X=B,  or the minimum norm
/// solution of an under-determined system, where A is a general
/// rectangular matrix of full rank,  using a QR or LQ factorization
/// of A.
void gels(char trans, f_int m, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *work, f_int lwork, ref f_int info) {
    sgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info, 1);
}
void gels(char trans, f_int m, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *work, f_int lwork, ref f_int info) {
    dgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info, 1);
}
void gels(char trans, f_int m, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *work, f_int lwork, ref f_int info) {
    cgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info, 1);
}
void gels(char trans, f_int m, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *work, f_int lwork, ref f_int info) {
    zgels_(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info, 1);
}

/// Computes the least squares solution to an over-determined system
/// of linear equations, A X=B or A**H X=B,  or the minimum norm
/// solution of an under-determined system, using a divide and conquer
/// method, where A is a general rectangular matrix of full rank,
/// using a QR or LQ factorization of A.
void gelsd(f_int m, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *s, f_float rcond, out f_int rank, f_float *work, f_int lwork, f_int *iwork, ref f_int info) {
    sgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
}
void gelsd(f_int m, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *s, f_double rcond, out f_int rank, f_double *work, f_int lwork, f_int *iwork, ref f_int info) {
    dgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
}
void gelsd(f_int m, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_float *s, f_float rcond, out f_int rank, f_cfloat *work, f_int lwork, f_float *rwork, f_int *iwork, ref f_int info) {
    cgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);
}
void gelsd(f_int m, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_double *s, f_double rcond, out f_int rank, f_cdouble *work, f_int lwork, f_double *rwork, f_int *iwork, ref f_int info) {
    zgelsd_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, rwork, iwork, &info);
}

/// Solves the LSE (Constrained Linear Least Squares Problem) using
/// the GRQ (Generalized RQ) factorization
void gglse(f_int m, f_int n, f_int p, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *c, f_float *d, f_float *x, f_float *work, f_int lwork, ref f_int info) {
    sgglse_(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
}
void gglse(f_int m, f_int n, f_int p, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *c, f_double *d, f_double *x, f_double *work, f_int lwork, ref f_int info) {
    dgglse_(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
}
void gglse(f_int m, f_int n, f_int p, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *c, f_cfloat *d, f_cfloat *x, f_cfloat *work, f_int lwork, ref f_int info) {
    cgglse_(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
}
void gglse(f_int m, f_int n, f_int p, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *c, f_cdouble *d, f_cdouble *x, f_cdouble *work, f_int lwork, ref f_int info) {
    zgglse_(&m, &n, &p, a, &lda, b, &ldb, c, d, x, work, &lwork, &info);
}

/// Solves the GLM (Generalized Linear Regression Model) using
/// the GQR (Generalized QR) factorization
void ggglm(f_int n, f_int m, f_int p, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *d, f_float *x, f_float *y, f_float *work, f_int lwork, ref f_int info) {
    sggglm_(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
}
void ggglm(f_int n, f_int m, f_int p, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *d, f_double *x, f_double *y, f_double *work, f_int lwork, ref f_int info) {
    dggglm_(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
}
void ggglm(f_int n, f_int m, f_int p, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *d, f_cfloat *x, f_cfloat *y, f_cfloat *work, f_int lwork, ref f_int info) {
    cggglm_(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
}
void ggglm(f_int n, f_int m, f_int p, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *d, f_cdouble *x, f_cdouble *y, f_cdouble *work, f_int lwork, ref f_int info) {
    zggglm_(&n, &m, &p, a, &lda, b, &ldb, d, x, y, work, &lwork, &info);
}

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix.
void syev(char jobz, char uplo, f_int n, f_float *a, f_int lda, f_float *w, f_float *work, f_int lwork, ref f_int info) {
    ssyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info, 1, 1);
}
void syev(char jobz, char uplo, f_int n, f_double *a, f_int lda, f_double *w, f_double *work, f_int lwork, ref f_int info) {
    dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info, 1, 1);
}

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian matrix.
void heev(char jobz, char uplo, f_int n, f_cfloat *a, f_int lda, f_float *w, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info, 1, 1);
}
void heev(char jobz, char uplo, f_int n, f_cdouble *a, f_int lda, f_double *w, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zheev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &info, 1, 1);
}


/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix.  If eigenvectors are desired, it uses a divide
/// and conquer algorithm.
void syevd(char jobz, char uplo, f_int n, f_float *a, f_int lda, f_float *w, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    ssyevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void syevd(char jobz, char uplo, f_int n, f_double *a, f_int lda, f_double *w, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dsyevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian matrix.  If eigenvectors are desired, it uses a divide
/// and conquer algorithm.
void heevd(char jobz, char uplo, f_int n, f_cfloat *a, f_int lda, f_float *w, f_cfloat *work, f_int lwork, f_float *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    cheevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}
void heevd(char jobz, char uplo, f_int n, f_cdouble *a, f_int lda, f_double *w, f_cdouble *work, f_int lwork, f_double *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    zheevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix in packed storage.
void spev(char jobz, char uplo, f_int n, f_float *ap, f_float *w, f_float *z, f_int ldz, f_float *work, ref f_int info) {
    sspev_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info, 1, 1);
}
void spev(char jobz, char uplo, f_int n, f_double *ap, f_double *w, f_double *z, f_int ldz, f_double *work, ref f_int info) {
    dspev_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &info, 1, 1);
}

/// Computes selected eigenvalues, and optionally, eigenvectors of a complex
/// Hermitian matrix.  Eigenvalues are computed by the dqds
/// algorithm, and eigenvectors are computed from various "good" LDL^T
/// representations (also known as Relatively Robust Representations).
/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian matrix in packed storage.
void hpev(char jobz, char uplo, f_int n, f_cfloat *ap, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_float *rwork, ref f_int info) {
    chpev_(&jobz, &uplo, &n, ap, w, z, &ldz, work, rwork, &info, 1, 1);
}
void hpev(char jobz, char uplo, f_int n, f_cdouble *ap, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_double *rwork, ref f_int info) {
    zhpev_(&jobz, &uplo, &n, ap, w, z, &ldz, work, rwork, &info, 1, 1);
}

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix in packed storage.  If eigenvectors are desired,
/// it uses a divide and conquer algorithm.
void spevd(char jobz, char uplo, f_int n, f_float *ap, f_float *w, f_float *z, f_int ldz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    sspevd_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void spevd(char jobz, char uplo, f_int n, f_double *ap, f_double *w, f_double *z, f_int ldz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dspevd_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian matrix in packed storage.  If eigenvectors are desired, it
/// uses a divide and conquer algorithm.
void hpevd(char jobz, char uplo, f_int n, f_cfloat *ap, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, f_float *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    chpevd_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}
void hpevd(char jobz, char uplo, f_int n, f_cdouble *ap, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, f_double *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    zhpevd_(&jobz, &uplo, &n, ap, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric band matrix.
void sbev(char jobz, char uplo, f_int n, f_int kd, f_float *ab, f_int ldab, f_float *w, f_float *z, f_int ldz, f_float *work, ref f_int info) {
    ssbev_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info, 1, 1);
}
void sbev(char jobz, char uplo, f_int n, f_int kd, f_double *ab, f_int ldab, f_double *w, f_double *z, f_int ldz, f_double *work, ref f_int info) {
    dsbev_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &info, 1, 1);
}

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian band matrix.
void hbev(char jobz, char uplo, f_int n, f_int kd, f_cfloat *ab, f_int ldab, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_float *rwork, ref f_int info) {
    chbev_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, rwork, &info, 1, 1);
}
void hbev(char jobz, char uplo, f_int n, f_int kd, f_cdouble *ab, f_int ldab, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_double *rwork, ref f_int info) {
    zhbev_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, rwork, &info, 1, 1);
}

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric band matrix.  If eigenvectors are desired, it uses a
/// divide and conquer algorithm.
void sbevd(char jobz, char uplo, f_int n, f_int kd, f_float *ab, f_int ldab, f_float *w, f_float *z, f_int ldz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    ssbevd_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void sbevd(char jobz, char uplo, f_int n, f_int kd, f_double *ab, f_int ldab, f_double *w, f_double *z, f_int ldz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dsbevd_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all eigenvalues and, optionally, eigenvectors of a complex
/// Hermitian band matrix.  If eigenvectors are desired, it uses a divide
/// and conquer algorithm.
void hbevd(char jobz, char uplo, f_int n, f_int kd, f_cfloat *ab, f_int ldab, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, f_float *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    chbevd_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}
void hbevd(char jobz, char uplo, f_int n, f_int kd, f_cdouble *ab, f_int ldab, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, f_double *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    zhbevd_(&jobz, &uplo, &n, &kd, ab, &ldab, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric tridiagonal matrix.
void stev(char jobz, f_int n, f_float *d, f_float *e, f_float *z, f_int ldz, f_float *work, ref f_int info) {
    sstev_(&jobz, &n, d, e, z, &ldz, work, &info, 1);
}
void stev(char jobz, f_int n, f_double *d, f_double *e, f_double *z, f_int ldz, f_double *work, ref f_int info) {
    dstev_(&jobz, &n, d, e, z, &ldz, work, &info, 1);
}

/// Computes all eigenvalues, and optionally, eigenvectors of a real
/// symmetric tridiagonal matrix.  If eigenvectors are desired, it uses
/// a divide and conquer algorithm.
void stevd(char jobz, f_int n, f_float *d, f_float *e, f_float *z, f_int ldz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    sstevd_(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info, 1);
}
void stevd(char jobz, f_int n, f_double *d, f_double *e, f_double *z, f_int ldz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dstevd_(&jobz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info, 1);
}

/// Computes the eigenvalues and Schur factorization of a general
/// matrix, and orders the factorization so that selected eigenvalues
/// are at the top left of the Schur form.
void gees(char jobvs, char sort, FCB_SGEES_SELECT select, f_int n, f_float *a, f_int lda, f_int sdim, f_float *wr, f_float *wi, f_float *vs, f_int ldvs, f_float *work, f_int lwork, f_int bwork, ref f_int info) {
    sgees_(&jobvs, &sort, select, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, work, &lwork, &bwork, &info, 1, 1);
}
void gees(char jobvs, char sort, FCB_DGEES_SELECT select, f_int n, f_double *a, f_int lda, f_int sdim, f_double *wr, f_double *wi, f_double *vs, f_int ldvs, f_double *work, f_int lwork, f_int bwork, ref f_int info) {
    dgees_(&jobvs, &sort, select, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, work, &lwork, &bwork, &info, 1, 1);
}
void gees(char jobvs, char sort, FCB_CGEES_SELECT select, f_int n, f_cfloat *a, f_int lda, f_int sdim, f_cfloat *w, f_cfloat *vs, f_int ldvs, f_cfloat *work, f_int lwork, f_float *rwork, f_int bwork, ref f_int info) {
    cgees_(&jobvs, &sort, select, &n, a, &lda, &sdim, w, vs, &ldvs, work, &lwork, rwork, &bwork, &info, 1, 1);
}
void gees(char jobvs, char sort, FCB_ZGEES_SELECT select, f_int n, f_cdouble *a, f_int lda, f_int sdim, f_cdouble *w, f_cdouble *vs, f_int ldvs, f_cdouble *work, f_int lwork, f_double *rwork, f_int bwork, ref f_int info) {
    zgees_(&jobvs, &sort, select, &n, a, &lda, &sdim, w, vs, &ldvs, work, &lwork, rwork, &bwork, &info, 1, 1);
}

/// Computes the eigenvalues and left and right eigenvectors of
/// a general matrix.
void geev(char jobvl, char jobvr, f_int n, f_float *a, f_int lda, f_float *wr, f_float *wi, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_float *work, f_int lwork, ref f_int info) {
    sgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info, 1, 1);
}
void geev(char jobvl, char jobvr, f_int n, f_double *a, f_int lda, f_double *wr, f_double *wi, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_double *work, f_int lwork, ref f_int info) {
    dgeev_(&jobvl, &jobvr, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, work, &lwork, &info, 1, 1);
}
void geev(char jobvl, char jobvr, f_int n, f_cfloat *a, f_int lda, f_cfloat *w, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info, 1, 1);
}
void geev(char jobvl, char jobvr, f_int n, f_cdouble *a, f_int lda, f_cdouble *w, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zgeev_(&jobvl, &jobvr, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info, 1, 1);
}

/// Computes the singular value decomposition (SVD) of a general
/// rectangular matrix.
void gesvd(char jobu, char jobvt, f_int m, f_int n, f_float *a, f_int lda, f_float *s, f_float *u, f_int ldu, f_float *vt, f_int ldvt, f_float *work, f_int lwork, ref f_int info) {
    sgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info, 1, 1);
}
void gesvd(char jobu, char jobvt, f_int m, f_int n, f_double *a, f_int lda, f_double *s, f_double *u, f_int ldu, f_double *vt, f_int ldvt, f_double *work, f_int lwork, ref f_int info) {
    dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, &info, 1, 1);
}
void gesvd(char jobu, char jobvt, f_int m, f_int n, f_cfloat *a, f_int lda, f_float *s, f_cfloat *u, f_int ldu, f_cfloat *vt, f_int ldvt, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info, 1, 1);
}
void gesvd(char jobu, char jobvt, f_int m, f_int n, f_cdouble *a, f_int lda, f_double *s, f_cdouble *u, f_int ldu, f_cdouble *vt, f_int ldvt, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, &info, 1, 1);
}

/// Computes the singular value decomposition (SVD) of a general
/// rectangular matrix using divide-and-conquer.
void gesdd(char jobz, f_int m, f_int n, f_float *a, f_int lda, f_float *s, f_float *u, f_int ldu, f_float *vt, f_int ldvt, f_float *work, f_int lwork, f_int *iwork, ref f_int info) {
    sgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info, 1);
}
void gesdd(char jobz, f_int m, f_int n, f_double *a, f_int lda, f_double *s, f_double *u, f_int ldu, f_double *vt, f_int ldvt, f_double *work, f_int lwork, f_int *iwork, ref f_int info) {
    dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, &info, 1);
}
void gesdd(char jobz, f_int m, f_int n, f_cfloat *a, f_int lda, f_float *s, f_cfloat *u, f_int ldu, f_cfloat *vt, f_int ldvt, f_cfloat *work, f_int lwork, f_float *rwork, f_int *iwork, ref f_int info) {
    cgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info, 1);
}
void gesdd(char jobz, f_int m, f_int n, f_cdouble *a, f_int lda, f_double *s, f_cdouble *u, f_int ldu, f_cdouble *vt, f_int ldvt, f_cdouble *work, f_int lwork, f_double *rwork, f_int *iwork, ref f_int info) {
    zgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, rwork, iwork, &info, 1);
}

/// Computes all eigenvalues and the eigenvectors of  a generalized
/// symmetric-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
void sygv(f_int itype, char jobz, char uplo, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *w, f_float *work, f_int lwork, ref f_int info) {
    ssygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info, 1, 1);
}
void sygv(f_int itype, char jobz, char uplo, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *w, f_double *work, f_int lwork, ref f_int info) {
    dsygv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, &info, 1, 1);
}

/// Computes all eigenvalues and the eigenvectors of  a generalized
/// Hermitian-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
void hegv(f_int itype, char jobz, char uplo, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_float *w, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    chegv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &info, 1, 1);
}
void hegv(f_int itype, char jobz, char uplo, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_double *w, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zhegv_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &info, 1, 1);
}

/// Computes all eigenvalues and the eigenvectors of  a generalized
/// symmetric-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void sygvd(f_int itype, char jobz, char uplo, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *w, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    ssygvd_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void sygvd(f_int itype, char jobz, char uplo, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *w, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dsygvd_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, iwork, &liwork, &info, 1, 1);
}
/// Computes all eigenvalues and the eigenvectors of  a generalized
/// Hermitian-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void hegvd(f_int itype, char jobz, char uplo, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_float *w, f_cfloat *work, f_int lwork, f_float *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    chegvd_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}
void hegvd(f_int itype, char jobz, char uplo, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_double *w, f_cdouble *work, f_int lwork, f_double *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    zhegvd_(&itype, &jobz, &uplo, &n, a, &lda, b, &ldb, w, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all eigenvalues and eigenvectors of  a generalized
/// symmetric-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
void spgv(f_int itype, char jobz, char uplo, f_int n, f_float *ap, f_float *bp, f_float *w, f_float *z, f_int ldz, f_float *work, ref f_int info) {
    sspgv_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &info, 1, 1);
}
void spgv(f_int itype, char jobz, char uplo, f_int n, f_double *ap, f_double *bp, f_double *w, f_double *z, f_int ldz, f_double *work, ref f_int info) {
    dspgv_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &info, 1, 1);
}

/// Computes all eigenvalues and eigenvectors of  a generalized
/// Hermitian-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
void hpgv(f_int itype, char jobz, char uplo, f_int n, f_cfloat *ap, f_cfloat *bp, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_float *rwork, ref f_int info) {
    chpgv_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, rwork, &info, 1, 1);
}
void hpgv(f_int itype, char jobz, char uplo, f_int n, f_cdouble *ap, f_cdouble *bp, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_double *rwork, ref f_int info) {
    zhpgv_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, rwork, &info, 1, 1);
}

/// Computes all eigenvalues and eigenvectors of  a generalized
/// symmetric-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void spgvd(f_int itype, char jobz, char uplo, f_int n, f_float *ap, f_float *bp, f_float *w, f_float *z, f_int ldz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    sspgvd_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void spgvd(f_int itype, char jobz, char uplo, f_int n, f_double *ap, f_double *bp, f_double *w, f_double *z, f_int ldz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dspgvd_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all eigenvalues and eigenvectors of  a generalized
/// Hermitian-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void hpgvd(f_int itype, char jobz, char uplo, f_int n, f_cfloat *ap, f_cfloat *bp, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, f_float *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    chpgvd_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}
void hpgvd(f_int itype, char jobz, char uplo, f_int n, f_cdouble *ap, f_cdouble *bp, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, f_double *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    zhpgvd_(&itype, &jobz, &uplo, &n, ap, bp, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all the eigenvalues, and optionally, the eigenvectors
/// of a real generalized symmetric-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be symmetric
/// and banded, and B is also positive definite.
void sbgv(char jobz, char uplo, f_int n, f_int ka, f_int kb, f_float *ab, f_int ldab, f_float *bb, f_int ldbb, f_float *w, f_float *z, f_int ldz, f_float *work, ref f_int info) {
    ssbgv_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &info, 1, 1);
}
void sbgv(char jobz, char uplo, f_int n, f_int ka, f_int kb, f_double *ab, f_int ldab, f_double *bb, f_int ldbb, f_double *w, f_double *z, f_int ldz, f_double *work, ref f_int info) {
    dsbgv_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &info, 1, 1);
}

/// Computes all the eigenvalues, and optionally, the eigenvectors
/// of a complex generalized Hermitian-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be Hermitian
/// and banded, and B is also positive definite.
void hbgv(char jobz, char uplo, f_int n, f_int ka, f_int kb, f_cfloat *ab, f_int ldab, f_cfloat *bb, f_int ldbb, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_float *rwork, ref f_int info) {
    chbgv_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, rwork, &info, 1, 1);
}
void hbgv(char jobz, char uplo, f_int n, f_int ka, f_int kb, f_cdouble *ab, f_int ldab, f_cdouble *bb, f_int ldbb, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_double *rwork, ref f_int info) {
    zhbgv_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, rwork, &info, 1, 1);
}

/// Computes all the eigenvalues, and optionally, the eigenvectors
/// of a real generalized symmetric-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be symmetric
/// and banded, and B is also positive definite.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void sbgvd(char jobz, char uplo, f_int n, f_int ka, f_int kb, f_float *ab, f_int ldab, f_float *bb, f_int ldbb, f_float *w, f_float *z, f_int ldz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    ssbgvd_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void sbgvd(char jobz, char uplo, f_int n, f_int ka, f_int kb, f_double *ab, f_int ldab, f_double *bb, f_int ldbb, f_double *w, f_double *z, f_int ldz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dsbgvd_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, iwork, &liwork, &info, 1, 1);
}

/// Computes all the eigenvalues, and optionally, the eigenvectors
/// of a complex generalized Hermitian-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be Hermitian
/// and banded, and B is also positive definite.
/// If eigenvectors are desired, it uses a divide and conquer algorithm.
void hbgvd(char jobz, char uplo, f_int n, f_int ka, f_int kb, f_cfloat *ab, f_int ldab, f_cfloat *bb, f_int ldbb, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, f_float *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    chbgvd_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}
void hbgvd(char jobz, char uplo, f_int n, f_int ka, f_int kb, f_cdouble *ab, f_int ldab, f_cdouble *bb, f_int ldbb, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, f_double *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    zhbgvd_(&jobz, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, w, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1);
}

/// Computes the generalized eigenvalues, Schur form, and left and/or
/// right Schur vectors for a pair of nonsymmetric matrices
void gegs(char jobvsl, char jobvsr, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *vsl, f_int ldvsl, f_float *vsr, f_int ldvsr, f_float *work, f_int lwork, ref f_int info) {
    sgegs_(&jobvsl, &jobvsr, &n, a, &lda, b, &ldb, alphar, alphai, betav, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, &info, 1, 1);
}
void gegs(char jobvsl, char jobvsr, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *vsl, f_int ldvsl, f_double *vsr, f_int ldvsr, f_double *work, f_int lwork, ref f_int info) {
    dgegs_(&jobvsl, &jobvsr, &n, a, &lda, b, &ldb, alphar, alphai, betav, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, &info, 1, 1);
}
void gegs(char jobvsl, char jobvsr, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vsl, f_int ldvsl, f_cfloat *vsr, f_int ldvsr, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cgegs_(&jobvsl, &jobvsr, &n, a, &lda, b, &ldb, alphav, betav, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, &info, 1, 1);
}
void gegs(char jobvsl, char jobvsr, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vsl, f_int ldvsl, f_cdouble *vsr, f_int ldvsr, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zgegs_(&jobvsl, &jobvsr, &n, a, &lda, b, &ldb, alphav, betav, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, &info, 1, 1);
}

/// Computes the generalized eigenvalues, Schur form, and left and/or
/// right Schur vectors for a pair of nonsymmetric matrices
void gges(char jobvsl, char jobvsr, char sort, FCB_SGGES_SELCTG selctg, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_int sdim, f_float *alphar, f_float *alphai, f_float *betav, f_float *vsl, f_int ldvsl, f_float *vsr, f_int ldvsr, f_float *work, f_int lwork, f_int bwork, ref f_int info) {
    sgges_(&jobvsl, &jobvsr, &sort, selctg, &n, a, &lda, b, &ldb, &sdim, alphar, alphai, betav, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, &bwork, &info, 1, 1, 1);
}
void gges(char jobvsl, char jobvsr, char sort, FCB_DGGES_DELCTG delctg, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_int sdim, f_double *alphar, f_double *alphai, f_double *betav, f_double *vsl, f_int ldvsl, f_double *vsr, f_int ldvsr, f_double *work, f_int lwork, f_int bwork, ref f_int info) {
    dgges_(&jobvsl, &jobvsr, &sort, delctg, &n, a, &lda, b, &ldb, &sdim, alphar, alphai, betav, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, &bwork, &info, 1, 1, 1);
}
void gges(char jobvsl, char jobvsr, char sort, FCB_CGGES_SELCTG selctg, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_int sdim, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vsl, f_int ldvsl, f_cfloat *vsr, f_int ldvsr, f_cfloat *work, f_int lwork, f_float *rwork, f_int bwork, ref f_int info) {
    cgges_(&jobvsl, &jobvsr, &sort, selctg, &n, a, &lda, b, &ldb, &sdim, alphav, betav, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, &bwork, &info, 1, 1, 1);
}
void gges(char jobvsl, char jobvsr, char sort, FCB_ZGGES_DELCTG delctg, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_int sdim, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vsl, f_int ldvsl, f_cdouble *vsr, f_int ldvsr, f_cdouble *work, f_int lwork, f_double *rwork, f_int bwork, ref f_int info) {
    zgges_(&jobvsl, &jobvsr, &sort, delctg, &n, a, &lda, b, &ldb, &sdim, alphav, betav, vsl, &ldvsl, vsr, &ldvsr, work, &lwork, rwork, &bwork, &info, 1, 1, 1);
}

/// Computes the generalized eigenvalues, and left and/or right
/// generalized eigenvectors for a pair of nonsymmetric matrices
void gegv(char jobvl, char jobvr, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_float *work, f_int lwork, ref f_int info) {
    sgegv_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, betav, vl, &ldvl, vr, &ldvr, work, &lwork, &info, 1, 1);
}
void gegv(char jobvl, char jobvr, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_double *work, f_int lwork, ref f_int info) {
    dgegv_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, betav, vl, &ldvl, vr, &ldvr, work, &lwork, &info, 1, 1);
}
void gegv(char jobvl, char jobvr, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *alphar, f_cfloat *betav, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cgegv_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, betav, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info, 1, 1);
}
void gegv(char jobvl, char jobvr, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *alphar, f_cdouble *betav, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zgegv_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, betav, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info, 1, 1);
}

/// Computes the generalized eigenvalues, and left and/or right
/// generalized eigenvectors for a pair of nonsymmetric matrices
void ggev(char jobvl, char jobvr, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_float *work, f_int lwork, ref f_int info) {
    sggev_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, betav, vl, &ldvl, vr, &ldvr, work, &lwork, &info, 1, 1);
}
void ggev(char jobvl, char jobvr, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_double *work, f_int lwork, ref f_int info) {
    dggev_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphar, alphai, betav, vl, &ldvl, vr, &ldvr, work, &lwork, &info, 1, 1);
}
void ggev(char jobvl, char jobvr, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cggev_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphav, betav, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info, 1, 1);
}
void ggev(char jobvl, char jobvr, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zggev_(&jobvl, &jobvr, &n, a, &lda, b, &ldb, alphav, betav, vl, &ldvl, vr, &ldvr, work, &lwork, rwork, &info, 1, 1);
}

/// Computes the Generalized Singular Value Decomposition
void ggsvd(char jobu, char jobv, char jobq, f_int m, f_int n, f_int p, f_int k, f_int l, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *alphav, f_float *betav, f_float *u, f_int ldu, f_float *v, f_int ldv, f_float *q, f_int ldq, f_float *work, f_int *iwork, ref f_int info) {
    sggsvd_(&jobu, &jobv, &jobq, &m, &n, &p, &k, &l, a, &lda, b, &ldb, alphav, betav, u, &ldu, v, &ldv, q, &ldq, work, iwork, &info, 1, 1, 1);
}
void ggsvd(char jobu, char jobv, char jobq, f_int m, f_int n, f_int p, f_int k, f_int l, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *alphav, f_double *betav, f_double *u, f_int ldu, f_double *v, f_int ldv, f_double *q, f_int ldq, f_double *work, f_int *iwork, ref f_int info) {
    dggsvd_(&jobu, &jobv, &jobq, &m, &n, &p, &k, &l, a, &lda, b, &ldb, alphav, betav, u, &ldu, v, &ldv, q, &ldq, work, iwork, &info, 1, 1, 1);
}
void ggsvd(char jobu, char jobv, char jobq, f_int m, f_int n, f_int p, f_int k, f_int l, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_float *alphav, f_float *betav, f_cfloat *u, f_int ldu, f_cfloat *v, f_int ldv, f_cfloat *q, f_int ldq, f_cfloat *work, f_float *rwork, f_int *iwork, ref f_int info) {
    cggsvd_(&jobu, &jobv, &jobq, &m, &n, &p, &k, &l, a, &lda, b, &ldb, alphav, betav, u, &ldu, v, &ldv, q, &ldq, work, rwork, iwork, &info, 1, 1, 1);
}
void ggsvd(char jobu, char jobv, char jobq, f_int m, f_int n, f_int p, f_int k, f_int l, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_double *alphav, f_double *betav, f_cdouble *u, f_int ldu, f_cdouble *v, f_int ldv, f_cdouble *q, f_int ldq, f_cdouble *work, f_double *rwork, f_int *iwork, ref f_int info) {
    zggsvd_(&jobu, &jobv, &jobq, &m, &n, &p, &k, &l, a, &lda, b, &ldb, alphav, betav, u, &ldu, v, &ldv, q, &ldq, work, rwork, iwork, &info, 1, 1, 1);
}

//-----------------------------------------------------
//       ---- EXPERT and RRR DRIVER routines ----
//-----------------------------------------------------

/// Solves a general system of linear equations AX=B, A**T X=B
/// or A**H X=B, and provides an estimate of the condition number
/// and error bounds on the solution.
void gesvx(char fact, char trans, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *af, f_int ldaf, f_int *ipiv, char equed, f_float *r, f_float *c, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sgesvx_(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void gesvx(char fact, char trans, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *af, f_int ldaf, f_int *ipiv, char equed, f_double *r, f_double *c, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dgesvx_(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void gesvx(char fact, char trans, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *af, f_int ldaf, f_int *ipiv, char equed, f_float *r, f_float *c, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgesvx_(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}
void gesvx(char fact, char trans, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *af, f_int ldaf, f_int *ipiv, char equed, f_double *r, f_double *c, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgesvx_(&fact, &trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}

/// Solves a general banded system of linear equations AX=B,
/// A**T X=B or A**H X=B, and provides an estimate of the condition
/// number and error bounds on the solution.
void gbsvx(char fact, char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_float *ab, f_int ldab, f_float *afb, f_int ldafb, f_int *ipiv, char equed, f_float *r, f_float *c, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sgbsvx_(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void gbsvx(char fact, char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_double *ab, f_int ldab, f_double *afb, f_int ldafb, f_int *ipiv, char equed, f_double *r, f_double *c, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dgbsvx_(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void gbsvx(char fact, char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_cfloat *ab, f_int ldab, f_cfloat *afb, f_int ldafb, f_int *ipiv, char equed, f_float *r, f_float *c, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgbsvx_(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}
void gbsvx(char fact, char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_cdouble *ab, f_int ldab, f_cdouble *afb, f_int ldafb, f_int *ipiv, char equed, f_double *r, f_double *c, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgbsvx_(&fact, &trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, &equed, r, c, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}

/// Solves a general tridiagonal system of linear equations AX=B,
/// A**T X=B or A**H X=B, and provides an estimate of the condition
/// number  and error bounds on the solution.
void gtsvx(char fact, char trans, f_int n, f_int nrhs, f_float *dl, f_float *d, f_float *du, f_float *dlf, f_float *df, f_float *duf, f_float *du2, f_int *ipiv, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sgtsvx_(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1);
}
void gtsvx(char fact, char trans, f_int n, f_int nrhs, f_double *dl, f_double *d, f_double *du, f_double *dlf, f_double *df, f_double *duf, f_double *du2, f_int *ipiv, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dgtsvx_(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1);
}
void gtsvx(char fact, char trans, f_int n, f_int nrhs, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *dlf, f_cfloat *df, f_cfloat *duf, f_cfloat *du2, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgtsvx_(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1);
}
void gtsvx(char fact, char trans, f_int n, f_int nrhs, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *dlf, f_cdouble *df, f_cdouble *duf, f_cdouble *du2, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgtsvx_(&fact, &trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1);
}

/// Solves a symmetric positive definite system of linear
/// equations AX=B, and provides an estimate of the condition number
/// and error bounds on the solution.
void posvx(char fact, char uplo, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *af, f_int ldaf, char equed, f_float *s, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sposvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void posvx(char fact, char uplo, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *af, f_int ldaf, char equed, f_double *s, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dposvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void posvx(char fact, char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *af, f_int ldaf, char equed, f_float *s, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cposvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}
void posvx(char fact, char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *af, f_int ldaf, char equed, f_double *s, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zposvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}

/// Solves a symmetric positive definite system of linear
/// equations AX=B, where A is held in packed storage, and provides
/// an estimate of the condition number and error bounds on the
/// solution.
void ppsvx(char fact, char uplo, f_int n, f_int nrhs, f_float *ap, f_float *afp, char equed, f_float *s, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sppsvx_(&fact, &uplo, &n, &nrhs, ap, afp, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void ppsvx(char fact, char uplo, f_int n, f_int nrhs, f_double *ap, f_double *afp, char equed, f_double *s, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dppsvx_(&fact, &uplo, &n, &nrhs, ap, afp, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void ppsvx(char fact, char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *afp, char equed, f_float *s, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cppsvx_(&fact, &uplo, &n, &nrhs, ap, afp, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}
void ppsvx(char fact, char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *afp, char equed, f_double *s, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zppsvx_(&fact, &uplo, &n, &nrhs, ap, afp, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}

/// Solves a symmetric positive definite banded system
/// of linear equations AX=B, and provides an estimate of the condition
/// number and error bounds on the solution.
void pbsvx(char fact, char uplo, f_int n, f_int kd, f_int nrhs, f_float *ab, f_int ldab, f_float *afb, f_int ldafb, char equed, f_float *s, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    spbsvx_(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void pbsvx(char fact, char uplo, f_int n, f_int kd, f_int nrhs, f_double *ab, f_int ldab, f_double *afb, f_int ldafb, char equed, f_double *s, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dpbsvx_(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void pbsvx(char fact, char uplo, f_int n, f_int kd, f_int nrhs, f_cfloat *ab, f_int ldab, f_cfloat *afb, f_int ldafb, char equed, f_float *s, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cpbsvx_(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}
void pbsvx(char fact, char uplo, f_int n, f_int kd, f_int nrhs, f_cdouble *ab, f_int ldab, f_cdouble *afb, f_int ldafb, char equed, f_double *s, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zpbsvx_(&fact, &uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, &equed, s, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1, 1);
}

/// Solves a symmetric positive definite tridiagonal
/// system of linear equations AX=B, and provides an estimate of
/// the condition number and error bounds on the solution.
void ptsvx(char fact, f_int n, f_int nrhs, f_float *d, f_float *e, f_float *df, f_float *ef, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, ref f_int info) {
    sptsvx_(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &info, 1);
}
void ptsvx(char fact, f_int n, f_int nrhs, f_double *d, f_double *e, f_double *df, f_double *ef, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, ref f_int info) {
    dptsvx_(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &info, 1);
}
void ptsvx(char fact, f_int n, f_int nrhs, f_float *d, f_cfloat *e, f_float *df, f_cfloat *ef, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cptsvx_(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1);
}
void ptsvx(char fact, f_int n, f_int nrhs, f_double *d, f_cdouble *e, f_double *df, f_cdouble *ef, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zptsvx_(&fact, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1);
}

/// Solves a real symmetric
/// indefinite system  of linear equations AX=B, and provides an
/// estimate of the condition number and error bounds on the solution.
void sysvx(char fact, char uplo, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *af, f_int ldaf, f_int *ipiv, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, f_int lwork, f_int *iwork, ref f_int info) {
    ssysvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &lwork, iwork, &info, 1, 1);
}
void sysvx(char fact, char uplo, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *af, f_int ldaf, f_int *ipiv, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, f_int lwork, f_int *iwork, ref f_int info) {
    dsysvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &lwork, iwork, &info, 1, 1);
}
void sysvx(char fact, char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *af, f_int ldaf, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    csysvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &lwork, rwork, &info, 1, 1);
}
void sysvx(char fact, char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *af, f_int ldaf, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zsysvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &lwork, rwork, &info, 1, 1);
}

/// Solves a complex Hermitian
/// indefinite system  of linear equations AX=B, and provides an
/// estimate of the condition number and error bounds on the solution.
void hesvx(char fact, char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *af, f_int ldaf, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    chesvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &lwork, rwork, &info, 1, 1);
}
void hesvx(char fact, char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *af, f_int ldaf, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zhesvx_(&fact, &uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, &lwork, rwork, &info, 1, 1);
}

/// Solves a real symmetric
/// indefinite system of linear equations AX=B, where A is held
/// in packed storage, and provides an estimate of the condition
/// number and error bounds on the solution.
void spsvx(char fact, char uplo, f_int n, f_int nrhs, f_float *ap, f_float *afp, f_int *ipiv, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sspsvx_(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1);
}
void spsvx(char fact, char uplo, f_int n, f_int nrhs, f_double *ap, f_double *afp, f_int *ipiv, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dspsvx_(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, iwork, &info, 1, 1);
}
void spsvx(char fact, char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *afp, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cspsvx_(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1);
}
void spsvx(char fact, char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *afp, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zspsvx_(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1);
}

/// Solves a complex Hermitian
/// indefinite system of linear equations AX=B, where A is held
/// in packed storage, and provides an estimate of the condition
/// number and error bounds on the solution.
void hpsvx(char fact, char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *afp, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float rcond, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    chpsvx_(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1);
}
void hpsvx(char fact, char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *afp, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double rcond, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zhpsvx_(&fact, &uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, &rcond, ferr, berr, work, rwork, &info, 1, 1);
}

/// Computes the minimum norm least squares solution to an over-
/// or under-determined system of linear equations A X=B, using a
/// complete orthogonal factorization of A.
void gelsx(f_int m, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, f_int jpvt, f_float rcond, out f_int rank, f_float *work, ref f_int info) {
    sgelsx_(&m, &n, &nrhs, a, &lda, b, &ldb, &jpvt, &rcond, &rank, work, &info);
}
void gelsx(f_int m, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, f_int jpvt, f_double rcond, out f_int rank, f_double *work, ref f_int info) {
    dgelsx_(&m, &n, &nrhs, a, &lda, b, &ldb, &jpvt, &rcond, &rank, work, &info);
}
void gelsx(f_int m, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_int jpvt, f_float rcond, out f_int rank, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgelsx_(&m, &n, &nrhs, a, &lda, b, &ldb, &jpvt, &rcond, &rank, work, rwork, &info);
}
void gelsx(f_int m, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_int jpvt, f_double rcond, out f_int rank, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgelsx_(&m, &n, &nrhs, a, &lda, b, &ldb, &jpvt, &rcond, &rank, work, rwork, &info);
}

/// Computes the minimum norm least squares solution to an over-
/// or under-determined system of linear equations A X=B, using a
/// complete orthogonal factorization of A.
void gelsy(f_int m, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, f_int jpvt, f_float rcond, out f_int rank, f_float *work, f_int lwork, ref f_int info) {
    sgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, &jpvt, &rcond, &rank, work, &lwork, &info);
}
void gelsy(f_int m, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, f_int jpvt, f_double rcond, out f_int rank, f_double *work, f_int lwork, ref f_int info) {
    dgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, &jpvt, &rcond, &rank, work, &lwork, &info);
}
void gelsy(f_int m, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_int jpvt, f_float rcond, out f_int rank, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, &jpvt, &rcond, &rank, work, &lwork, rwork, &info);
}
void gelsy(f_int m, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_int jpvt, f_double rcond, out f_int rank, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zgelsy_(&m, &n, &nrhs, a, &lda, b, &ldb, &jpvt, &rcond, &rank, work, &lwork, rwork, &info);
}

/// Computes the minimum norm least squares solution to an over-
/// or under-determined system of linear equations A X=B,  using
/// the singular value decomposition of A.
void gelss(f_int m, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *s, f_float rcond, out f_int rank, f_float *work, f_int lwork, ref f_int info) {
    sgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info);
}
void gelss(f_int m, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *s, f_double rcond, out f_int rank, f_double *work, f_int lwork, ref f_int info) {
    dgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, &info);
}
void gelss(f_int m, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_float *s, f_float rcond, out f_int rank, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, rwork, &info);
}
void gelss(f_int m, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_double *s, f_double rcond, out f_int rank, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zgelss_(&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, rwork, &info);
}

/// Computes selected eigenvalues and eigenvectors of a symmetric matrix.
void syevx(char jobz, char range, char uplo, f_int n, f_float *a, f_int lda, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_float *work, f_int lwork, f_int *iwork, f_int ifail, ref f_int info) {
    ssyevx_(&jobz, &range, &uplo, &n, a, &lda, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, &lwork, iwork, &ifail, &info, 1, 1, 1);
}
void syevx(char jobz, char range, char uplo, f_int n, f_double *a, f_int lda, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_double *work, f_int lwork, f_int *iwork, f_int ifail, ref f_int info) {
    dsyevx_(&jobz, &range, &uplo, &n, a, &lda, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, &lwork, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues and eigenvectors of a Hermitian matrix.
void heevx(char jobz, char range, char uplo, f_int n, f_cfloat *a, f_int lda, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, f_float *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    cheevx_(&jobz, &range, &uplo, &n, a, &lda, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, &lwork, rwork, iwork, &ifail, &info, 1, 1, 1);
}
void heevx(char jobz, char range, char uplo, f_int n, f_cdouble *a, f_int lda, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, f_double *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    zheevx_(&jobz, &range, &uplo, &n, a, &lda, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, &lwork, rwork, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues, and optionally, eigenvectors of a real
/// symmetric matrix.  Eigenvalues are computed by the dqds
/// algorithm, and eigenvectors are computed from various "good" LDL^T
/// representations (also known as Relatively Robust Representations).
void syevr(char jobz, char range, char uplo, f_int n, f_float *a, f_int lda, f_float vl, f_float vu, f_int il, f_int iu, f_float abstol, out f_int m, f_float *w, f_float *z, f_int ldz, f_int *isuppz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    ssyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info, 1, 1, 1);
}
void syevr(char jobz, char range, char uplo, f_int n, f_double *a, f_int lda, f_double vl, f_double vu, f_int il, f_int iu, f_double abstol, out f_int m, f_double *w, f_double *z, f_int ldz, f_int *isuppz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dsyevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork, &info, 1, 1, 1);
}

/// Computes selected eigenvalues, and optionally, eigenvectors of a complex
/// Hermitian matrix.  Eigenvalues are computed by the dqds
/// algorithm, and eigenvectors are computed from various "good" LDL^T
/// representations (also known as Relatively Robust Representations).
void heevr(char jobz, char range, char uplo, f_int n, f_cfloat *a, f_int lda, f_float vl, f_float vu, f_int il, f_int iu, f_float abstol, out f_int m, f_float *w, f_cfloat *z, f_int ldz, f_int *isuppz, f_cfloat *work, f_int lwork, f_float *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    cheevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz,
        isuppz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1, 1);
}
void heevr(char jobz, char range, char uplo, f_int n, f_cdouble *a, f_int lda, f_double vl, f_double vu, f_int il, f_int iu, f_double abstol, out f_int m, f_double *w, f_cdouble *z, f_int ldz, f_int* isuppz, f_cdouble *work, f_int lwork, f_double *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    zheevr_(&jobz, &range, &uplo, &n, a, &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, isuppz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1, 1, 1);
}


/// Computes selected eigenvalues, and optionally, the eigenvectors of
/// a generalized symmetric-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
void sygvx(f_int itype, char jobz, char range, char uplo, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_float *work, f_int lwork, f_int *iwork, f_int ifail, ref f_int info) {
    ssygvx_(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, &lwork, iwork, &ifail, &info, 1, 1, 1);
}
void sygvx(f_int itype, char jobz, char range, char uplo, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_double *work, f_int lwork, f_int *iwork, f_int ifail, ref f_int info) {
    dsygvx_(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, &lwork, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues, and optionally, the eigenvectors of
/// a generalized Hermitian-definite generalized eigenproblem,
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x.
void hegvx(f_int itype, char jobz, char range, char uplo, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, f_float *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    chegvx_(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, &lwork, rwork, iwork, &ifail, &info, 1, 1, 1);
}
void hegvx(f_int itype, char jobz, char range, char uplo, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, f_double *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    zhegvx_(&itype, &jobz, &range, &uplo, &n, a, &lda, b, &ldb, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, &lwork, rwork, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues and eigenvectors of a
/// symmetric matrix in packed storage.
void spevx(char jobz, char range, char uplo, f_int n, f_float *ap, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_float *work, f_int *iwork, f_int ifail, ref f_int info) {
    sspevx_(&jobz, &range, &uplo, &n, ap, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1, 1);
}
void spevx(char jobz, char range, char uplo, f_int n, f_double *ap, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_double *work, f_int *iwork, f_int ifail, ref f_int info) {
    dspevx_(&jobz, &range, &uplo, &n, ap, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues and eigenvectors of a
/// Hermitian matrix in packed storage.
void hpevx(char jobz, char range, char uplo, f_int n, f_cfloat *ap, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_float *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    chpevx_(&jobz, &range, &uplo, &n, ap, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, rwork, iwork, &ifail, &info, 1, 1, 1);
}
void hpevx(char jobz, char range, char uplo, f_int n, f_cdouble *ap, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_double *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    zhpevx_(&jobz, &range, &uplo, &n, ap, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, rwork, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues, and optionally, eigenvectors of
/// a generalized symmetric-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
void spgvx(f_int itype, char jobz, char range, char uplo, f_int n, f_float *ap, f_float *bp, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_float *work, f_int *iwork, f_int ifail, ref f_int info) {
    sspgvx_(&itype, &jobz, &range, &uplo, &n, ap, bp, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1, 1);
}
void spgvx(f_int itype, char jobz, char range, char uplo, f_int n, f_double *ap, f_double *bp, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_double *work, f_int *iwork, f_int ifail, ref f_int info) {
    dspgvx_(&itype, &jobz, &range, &uplo, &n, ap, bp, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues, and optionally, the eigenvectors of
/// a generalized Hermitian-definite generalized eigenproblem,  Ax= lambda
/// Bx,  ABx= lambda x,  or BAx= lambda x, where A and B are in packed
/// storage.
void hpgvx(f_int itype, char jobz, char range, char uplo, f_int n, f_cfloat *ap, f_cfloat *bp, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_float *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    chpgvx_(&itype, &jobz, &range, &uplo, &n, ap, bp, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, rwork, iwork, &ifail, &info, 1, 1, 1);
}
void hpgvx(f_int itype, char jobz, char range, char uplo, f_int n, f_cdouble *ap, f_cdouble *bp, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_double *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    zhpgvx_(&itype, &jobz, &range, &uplo, &n, ap, bp, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, rwork, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues and eigenvectors of a
/// symmetric band matrix.
void sbevx(char jobz, char range, char uplo, f_int n, f_int kd, f_float *ab, f_int ldab, f_float *q, f_int ldq, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_float *work, f_int *iwork, f_int ifail, ref f_int info) {
    ssbevx_(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1, 1);
}
void sbevx(char jobz, char range, char uplo, f_int n, f_int kd, f_double *ab, f_int ldab, f_double *q, f_int ldq, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_double *work, f_int *iwork, f_int ifail, ref f_int info) {
    dsbevx_(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues and eigenvectors of a
/// Hermitian band matrix.
void hbevx(char jobz, char range, char uplo, f_int n, f_int kd, f_cfloat *ab, f_int ldab, f_cfloat *q, f_int ldq, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_float *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    chbevx_(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, rwork, iwork, &ifail, &info, 1, 1, 1);
}
void hbevx(char jobz, char range, char uplo, f_int n, f_int kd, f_cdouble *ab, f_int ldab, f_cdouble *q, f_int ldq, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_double *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    zhbevx_(&jobz, &range, &uplo, &n, &kd, ab, &ldab, q, &ldq, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, rwork, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues, and optionally, the eigenvectors
/// of a real generalized symmetric-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be symmetric
/// and banded, and B is also positive definite.
void sbgvx(char jobz, char range, char uplo, f_int n, f_int ka, f_int kb, f_float *ab, f_int ldab, f_float *bb, f_int ldbb, f_float *q, f_int ldq, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_float *work, f_int *iwork, f_int ifail, ref f_int info) {
    ssbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1, 1);
}
void sbgvx(char jobz, char range, char uplo, f_int n, f_int ka, f_int kb, f_double *ab, f_int ldab, f_double *bb, f_int ldbb, f_double *q, f_int ldq, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_double *work, f_int *iwork, f_int ifail, ref f_int info) {
    dsbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues, and optionally, the eigenvectors
/// of a complex generalized Hermitian-definite banded eigenproblem, of
/// the form A*x=(lambda)*B*x.  A and B are assumed to be Hermitian
/// and banded, and B is also positive definite.
void hbgvx(char jobz, char range, char uplo, f_int n, f_int ka, f_int kb, f_cfloat *ab, f_int ldab, f_cfloat *bb, f_int ldbb, f_cfloat *q, f_int ldq, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_float *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    chbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, rwork, iwork, &ifail, &info, 1, 1, 1);
}
void hbgvx(char jobz, char range, char uplo, f_int n, f_int ka, f_int kb, f_cdouble *ab, f_int ldab, f_cdouble *bb, f_int ldbb, f_cdouble *q, f_int ldq, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_double *rwork, f_int *iwork, f_int ifail, ref f_int info) {
    zhbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, q, &ldq, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, rwork, iwork, &ifail, &info, 1, 1, 1);
}

/// Computes selected eigenvalues and eigenvectors of a real
/// symmetric tridiagonal matrix.
void stevx(char jobz, char range, f_int n, f_float *d, f_float *e, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_float *work, f_int *iwork, f_int ifail, ref f_int info) {
    sstevx_(&jobz, &range, &n, d, e, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1);
}
void stevx(char jobz, char range, f_int n, f_double *d, f_double *e, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_double *work, f_int *iwork, f_int ifail, ref f_int info) {
    dstevx_(&jobz, &range, &n, d, e, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, work, iwork, &ifail, &info, 1, 1);
}

/// Computes selected eigenvalues, and optionally, eigenvectors of a real
/// symmetric tridiagonal matrix.  Eigenvalues are computed by the dqds
/// algorithm, and eigenvectors are computed from various "good" LDL^T
/// representations (also known as Relatively Robust Representations).
void stevr(char jobz, char range, f_int n, f_float *d, f_float *e, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_int isuppz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    sstevr_(&jobz, &range, &n, d, e, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, &isuppz, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void stevr(char jobz, char range, f_int n, f_double *d, f_double *e, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_int isuppz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dstevr_(&jobz, &range, &n, d, e, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, &isuppz, work, &lwork, iwork, &liwork, &info, 1, 1);
}

/// Computes the eigenvalues and Schur factorization of a general
/// matrix, orders the factorization so that selected eigenvalues
/// are at the top left of the Schur form, and computes reciprocal
/// condition numbers for the average of the selected eigenvalues,
/// and for the associated right invariant subspace.
void geesx(char jobvs, char sort, FCB_SGEESX_SELECT select, char sense, f_int n, f_float *a, f_int lda, f_int sdim, f_float *wr, f_float *wi, f_float *vs, f_int ldvs, f_float *rconde, f_float *rcondv, f_float *work, f_int lwork, f_int *iwork, f_int liwork, f_int bwork, ref f_int info) {
    sgeesx_(&jobvs, &sort, select, &sense, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, rconde, rcondv, work, &lwork, iwork, &liwork, &bwork, &info, 1, 1, 1);
}
void geesx(char jobvs, char sort, FCB_DGEESX_SELECT select, char sense, f_int n, f_double *a, f_int lda, f_int sdim, f_double *wr, f_double *wi, f_double *vs, f_int ldvs, f_double *rconde, f_double *rcondv, f_double *work, f_int lwork, f_int *iwork, f_int liwork, f_int bwork, ref f_int info) {
    dgeesx_(&jobvs, &sort, select, &sense, &n, a, &lda, &sdim, wr, wi, vs, &ldvs, rconde, rcondv, work, &lwork, iwork, &liwork, &bwork, &info, 1, 1, 1);
}
void geesx(char jobvs, char sort, FCB_CGEESX_SELECT select, char sense, f_int n, f_cfloat *a, f_int lda, f_int sdim, f_cfloat *w, f_cfloat *vs, f_int ldvs, f_float *rconde, f_float *rcondv, f_cfloat *work, f_int lwork, f_float *rwork, f_int bwork, ref f_int info) {
    cgeesx_(&jobvs, &sort, select, &sense, &n, a, &lda, &sdim, w, vs, &ldvs, rconde, rcondv, work, &lwork, rwork, &bwork, &info, 1, 1, 1);
}
void geesx(char jobvs, char sort, FCB_ZGEESX_SELECT select, char sense, f_int n, f_cdouble *a, f_int lda, f_int sdim, f_cdouble *w, f_cdouble *vs, f_int ldvs, f_double *rconde, f_double *rcondv, f_cdouble *work, f_int lwork, f_double *rwork, f_int bwork, ref f_int info) {
    zgeesx_(&jobvs, &sort, select, &sense, &n, a, &lda, &sdim, w, vs, &ldvs, rconde, rcondv, work, &lwork, rwork, &bwork, &info, 1, 1, 1);
}

/// Computes the generalized eigenvalues, the real Schur form, and,
/// optionally, the left and/or right matrices of Schur vectors.
void ggesx(char jobvsl, char jobvsr, char sort, FCB_SGGESX_SELCTG selctg, char sense, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_int sdim, f_float *alphar, f_float *alphai, f_float *betav, f_float *vsl, f_int ldvsl, f_float *vsr, f_int ldvsr, f_float *rconde, f_float *rcondv, f_float *work, f_int lwork, f_int *iwork, f_int liwork, f_int bwork, ref f_int info) {
    sggesx_(&jobvsl, &jobvsr, &sort, selctg, &sense, &n, a, &lda, b, &ldb, &sdim, alphar, alphai, betav, vsl, &ldvsl, vsr, &ldvsr, rconde, rcondv, work, &lwork, iwork, &liwork, &bwork, &info, 1, 1, 1, 1);
}
void ggesx(char jobvsl, char jobvsr, char sort, FCB_DGGESX_DELCTG delctg, char sense, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_int sdim, f_double *alphar, f_double *alphai, f_double *betav, f_double *vsl, f_int ldvsl, f_double *vsr, f_int ldvsr, f_double *rconde, f_double *rcondv, f_double *work, f_int lwork, f_int *iwork, f_int liwork, f_int bwork, ref f_int info) {
    dggesx_(&jobvsl, &jobvsr, &sort, delctg, &sense, &n, a, &lda, b, &ldb, &sdim, alphar, alphai, betav, vsl, &ldvsl, vsr, &ldvsr, rconde, rcondv, work, &lwork, iwork, &liwork, &bwork, &info, 1, 1, 1, 1);
}
void ggesx(char jobvsl, char jobvsr, char sort, FCB_CGGESX_SELCTG selctg, char sense, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_int sdim, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vsl, f_int ldvsl, f_cfloat *vsr, f_int ldvsr, f_float *rconde, f_float *rcondv, f_cfloat *work, f_int lwork, f_float *rwork, f_int *iwork, f_int liwork, f_int bwork, ref f_int info) {
    cggesx_(&jobvsl, &jobvsr, &sort, selctg, &sense, &n, a, &lda, b, &ldb, &sdim, alphav, betav, vsl, &ldvsl, vsr, &ldvsr, rconde, rcondv, work, &lwork, rwork, iwork, &liwork, &bwork, &info, 1, 1, 1, 1);
}
void ggesx(char jobvsl, char jobvsr, char sort, FCB_ZGGESX_DELCTG delctg, char sense, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_int sdim, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vsl, f_int ldvsl, f_cdouble *vsr, f_int ldvsr, f_double *rconde, f_double *rcondv, f_cdouble *work, f_int lwork, f_double *rwork, f_int *iwork, f_int liwork, f_int bwork, ref f_int info) {
    zggesx_(&jobvsl, &jobvsr, &sort, delctg, &sense, &n, a, &lda, b, &ldb, &sdim, alphav, betav, vsl, &ldvsl, vsr, &ldvsr, rconde, rcondv, work, &lwork, rwork, iwork, &liwork, &bwork, &info, 1, 1, 1, 1);
}

/// Computes the eigenvalues and left and right eigenvectors of
/// a general matrix,  with preliminary balancing of the matrix,
/// and computes reciprocal condition numbers for the eigenvalues
/// and right eigenvectors.
void geevx(char balanc, char jobvl, char jobvr, char sense, f_int n, f_float *a, f_int lda, f_float *wr, f_float *wi, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_int ilo, f_int ihi, f_float *scale, f_float *abnrm, f_float *rconde, f_float *rcondv, f_float *work, f_int lwork, f_int *iwork, ref f_int info) {
    sgeevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, abnrm, rconde, rcondv, work, &lwork, iwork, &info, 1, 1, 1, 1);
}
void geevx(char balanc, char jobvl, char jobvr, char sense, f_int n, f_double *a, f_int lda, f_double *wr, f_double *wi, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_int ilo, f_int ihi, f_double *scale, f_double *abnrm, f_double *rconde, f_double *rcondv, f_double *work, f_int lwork, f_int *iwork, ref f_int info) {
    dgeevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, wr, wi, vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, abnrm, rconde, rcondv, work, &lwork, iwork, &info, 1, 1, 1, 1);
}
void geevx(char balanc, char jobvl, char jobvr, char sense, f_int n, f_cfloat *a, f_int lda, f_cfloat *w, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_int ilo, f_int ihi, f_float *scale, f_float *abnrm, f_float *rconde, f_float *rcondv, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cgeevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, abnrm, rconde, rcondv, work, &lwork, rwork, &info, 1, 1, 1, 1);
}
void geevx(char balanc, char jobvl, char jobvr, char sense, f_int n, f_cdouble *a, f_int lda, f_cdouble *w, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_int ilo, f_int ihi, f_double *scale, f_double *abnrm, f_double *rconde, f_double *rcondv, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zgeevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, w, vl, &ldvl, vr, &ldvr, &ilo, &ihi, scale, abnrm, rconde, rcondv, work, &lwork, rwork, &info, 1, 1, 1, 1);
}

/// Computes the generalized eigenvalues, and optionally, the left
/// and/or right generalized eigenvectors.
void ggevx(char balanc, char jobvl, char jobvr, char sense, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_int ilo, f_int ihi, f_float *lscale, f_float *rscale, f_float *abnrm, f_float *bbnrm, f_float *rconde, f_float *rcondv, f_float *work, f_int lwork, f_int *iwork, f_int bwork, ref f_int info) {
    sggevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphar, alphai, betav, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, iwork, &bwork, &info, 1, 1, 1, 1);
}
void ggevx(char balanc, char jobvl, char jobvr, char sense, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_int ilo, f_int ihi, f_double *lscale, f_double *rscale, f_double *abnrm, f_double *bbnrm, f_double *rconde, f_double *rcondv, f_double *work, f_int lwork, f_int *iwork, f_int bwork, ref f_int info) {
    dggevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphar, alphai, betav, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, iwork, &bwork, &info, 1, 1, 1, 1);
}
void ggevx(char balanc, char jobvl, char jobvr, char sense, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_int ilo, f_int ihi, f_float *lscale, f_float *rscale, f_float *abnrm, f_float *bbnrm, f_float *rconde, f_float *rcondv, f_cfloat *work, f_int lwork, f_float *rwork, f_int *iwork, f_int bwork, ref f_int info) {
    cggevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphav, betav, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, &bwork, &info, 1, 1, 1, 1);
}
void ggevx(char balanc, char jobvl, char jobvr, char sense, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_int ilo, f_int ihi, f_double *lscale, f_double *rscale, f_double *abnrm, f_double *bbnrm, f_double *rconde, f_double *rcondv, f_cdouble *work, f_int lwork, f_double *rwork, f_int *iwork, f_int bwork, ref f_int info) {
    zggevx_(&balanc, &jobvl, &jobvr, &sense, &n, a, &lda, b, &ldb, alphav, betav, vl, &ldvl, vr, &ldvr, &ilo, &ihi, lscale, rscale, abnrm, bbnrm, rconde, rcondv, work, &lwork, rwork, iwork, &bwork, &info, 1, 1, 1, 1);
}



//----------------------------------------
//    ---- COMPUTATIONAL routines ----
//----------------------------------------


/// Computes the singular value decomposition (SVD) of a real bidiagonal
/// matrix, using a divide and conquer method.
void bdsdc(char uplo, char compq, f_int n, f_float *d, f_float *e, f_float *u, f_int ldu, f_float *vt, f_int ldvt, f_float *q, f_int iq, f_float *work, f_int *iwork, ref f_int info) {
    sbdsdc_(&uplo, &compq, &n, d, e, u, &ldu, vt, &ldvt, q, &iq, work, iwork, &info, 1, 1);
}
void bdsdc(char uplo, char compq, f_int n, f_double *d, f_double *e, f_double *u, f_int ldu, f_double *vt, f_int ldvt, f_double *q, f_int iq, f_double *work, f_int *iwork, ref f_int info) {
    dbdsdc_(&uplo, &compq, &n, d, e, u, &ldu, vt, &ldvt, q, &iq, work, iwork, &info, 1, 1);
}

/// Computes the singular value decomposition (SVD) of a real bidiagonal
/// matrix, using the bidiagonal QR algorithm.
void bdsqr(char uplo, f_int n, f_int ncvt, f_int nru, f_int ncc, f_float *d, f_float *e, f_float *vt, f_int ldvt, f_float *u, f_int ldu, f_float *c, f_int ldc, f_float *work, ref f_int info) {
    sbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info, 1);
}
void bdsqr(char uplo, f_int n, f_int ncvt, f_int nru, f_int ncc, f_double *d, f_double *e, f_double *vt, f_int ldvt, f_double *u, f_int ldu, f_double *c, f_int ldc, f_double *work, ref f_int info) {
    dbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, work, &info, 1);
}
void bdsqr(char uplo, f_int n, f_int ncvt, f_int nru, f_int ncc, f_float *d, f_float *e, f_cfloat *vt, f_int ldvt, f_cfloat *u, f_int ldu, f_cfloat *c, f_int ldc, f_float *rwork, ref f_int info) {
    cbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, rwork, &info, 1);
}
void bdsqr(char uplo, f_int n, f_int ncvt, f_int nru, f_int ncc, f_double *d, f_double *e, f_cdouble *vt, f_int ldvt, f_cdouble *u, f_int ldu, f_cdouble *c, f_int ldc, f_double *rwork, ref f_int info) {
    zbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, d, e, vt, &ldvt, u, &ldu, c, &ldc, rwork, &info, 1);
}

/// Computes the reciprocal condition numbers for the eigenvectors of a
/// real symmetric or complex Hermitian matrix or for the left or right
/// singular vectors of a general matrix.
void disna(char job, f_int m, f_int n, f_float *d, f_float *sep, ref f_int info) {
    sdisna_(&job, &m, &n, d, sep, &info, 1);
}
void disna(char job, f_int m, f_int n, f_double *d, f_double *sep, ref f_int info) {
    ddisna_(&job, &m, &n, d, sep, &info, 1);
}

/// Reduces a general band matrix to real upper bidiagonal form
/// by an orthogonal transformation.
void gbbrd(char vect, f_int m, f_int n, f_int ncc, f_int kl, f_int ku, f_float *ab, f_int ldab, f_float *d, f_float *e, f_float *q, f_int ldq, f_float *pt, f_int ldpt, f_float *c, f_int ldc, f_float *work, ref f_int info) {
    sgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info, 1);
}
void gbbrd(char vect, f_int m, f_int n, f_int ncc, f_int kl, f_int ku, f_double *ab, f_int ldab, f_double *d, f_double *e, f_double *q, f_int ldq, f_double *pt, f_int ldpt, f_double *c, f_int ldc, f_double *work, ref f_int info) {
    dgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, &info, 1);
}
void gbbrd(char vect, f_int m, f_int n, f_int ncc, f_int kl, f_int ku, f_cfloat *ab, f_int ldab, f_float *d, f_float *e, f_cfloat *q, f_int ldq, f_cfloat *pt, f_int ldpt, f_cfloat *c, f_int ldc, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, rwork, &info, 1);
}
void gbbrd(char vect, f_int m, f_int n, f_int ncc, f_int kl, f_int ku, f_cdouble *ab, f_int ldab, f_double *d, f_double *e, f_cdouble *q, f_int ldq, f_cdouble *pt, f_int ldpt, f_cdouble *c, f_int ldc, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgbbrd_(&vect, &m, &n, &ncc, &kl, &ku, ab, &ldab, d, e, q, &ldq, pt, &ldpt, c, &ldc, work, rwork, &info, 1);
}

/// Estimates the reciprocal of the condition number of a general
/// band matrix, in either the 1-norm or the infinity-norm, using
/// the LU factorization computed by SGBTRF.
void gbcon(char norm, f_int n, f_int kl, f_int ku, f_float *ab, f_int ldab, f_int *ipiv, f_float *anorm, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    sgbcon_(&norm, &n, &kl, &ku, ab, &ldab, ipiv, anorm, &rcond, work, iwork, &info, 1);
}
void gbcon(char norm, f_int n, f_int kl, f_int ku, f_double *ab, f_int ldab, f_int *ipiv, f_double *anorm, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dgbcon_(&norm, &n, &kl, &ku, ab, &ldab, ipiv, anorm, &rcond, work, iwork, &info, 1);
}
void gbcon(char norm, f_int n, f_int kl, f_int ku, f_cfloat *ab, f_int ldab, f_int *ipiv, f_float *anorm, f_float rcond, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgbcon_(&norm, &n, &kl, &ku, ab, &ldab, ipiv, anorm, &rcond, work, rwork, &info, 1);
}
void gbcon(char norm, f_int n, f_int kl, f_int ku, f_cdouble *ab, f_int ldab, f_int *ipiv, f_double *anorm, f_double rcond, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgbcon_(&norm, &n, &kl, &ku, ab, &ldab, ipiv, anorm, &rcond, work, rwork, &info, 1);
}

/// Computes row and column scalings to equilibrate a general band
/// matrix and reduce its condition number.
void gbequ(f_int m, f_int n, f_int kl, f_int ku, f_float *ab, f_int ldab, f_float *r, f_float *c, f_float *rowcnd, f_float *colcnd, f_float *amax, ref f_int info) {
    sgbequ_(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
}
void gbequ(f_int m, f_int n, f_int kl, f_int ku, f_double *ab, f_int ldab, f_double *r, f_double *c, f_double *rowcnd, f_double *colcnd, f_double *amax, ref f_int info) {
    dgbequ_(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
}
void gbequ(f_int m, f_int n, f_int kl, f_int ku, f_cfloat *ab, f_int ldab, f_float *r, f_float *c, f_float *rowcnd, f_float *colcnd, f_float *amax, ref f_int info) {
    cgbequ_(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
}
void gbequ(f_int m, f_int n, f_int kl, f_int ku, f_cdouble *ab, f_int ldab, f_double *r, f_double *c, f_double *rowcnd, f_double *colcnd, f_double *amax, ref f_int info) {
    zgbequ_(&m, &n, &kl, &ku, ab, &ldab, r, c, rowcnd, colcnd, amax, &info);
}

/// Improves the computed solution to a general banded system of
/// linear equations AX=B, A**T X=B or A**H X=B, and provides forward
/// and backward error bounds for the solution.
void gbrfs(char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_float *ab, f_int ldab, f_float *afb, f_int ldafb, f_int *ipiv, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sgbrfs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void gbrfs(char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_double *ab, f_int ldab, f_double *afb, f_int ldafb, f_int *ipiv, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dgbrfs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void gbrfs(char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_cfloat *ab, f_int ldab, f_cfloat *afb, f_int ldafb, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgbrfs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void gbrfs(char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_cdouble *ab, f_int ldab, f_cdouble *afb, f_int ldafb, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgbrfs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, afb, &ldafb, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Computes an LU factorization of a general band matrix, using
/// partial pivoting with row interchanges.
void gbtrf(f_int m, f_int n, f_int kl, f_int ku, f_float *ab, f_int ldab, f_int *ipiv, ref f_int info) {
    sgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
}
void gbtrf(f_int m, f_int n, f_int kl, f_int ku, f_double *ab, f_int ldab, f_int *ipiv, ref f_int info) {
    dgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
}
void gbtrf(f_int m, f_int n, f_int kl, f_int ku, f_cfloat *ab, f_int ldab, f_int *ipiv, ref f_int info) {
    cgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
}
void gbtrf(f_int m, f_int n, f_int kl, f_int ku, f_cdouble *ab, f_int ldab, f_int *ipiv, ref f_int info) {
    zgbtrf_(&m, &n, &kl, &ku, ab, &ldab, ipiv, &info);
}

/// Solves a general banded system of linear equations AX=B,
/// A**T X=B or A**H X=B, using the LU factorization computed
/// by SGBTRF.
void gbtrs(char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_float *ab, f_int ldab, f_int *ipiv, f_float *b, f_int ldb, ref f_int info) {
    sgbtrs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info, 1);
}
void gbtrs(char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_double *ab, f_int ldab, f_int *ipiv, f_double *b, f_int ldb, ref f_int info) {
    dgbtrs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info, 1);
}
void gbtrs(char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_cfloat *ab, f_int ldab, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    cgbtrs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info, 1);
}
void gbtrs(char trans, f_int n, f_int kl, f_int ku, f_int nrhs, f_cdouble *ab, f_int ldab, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zgbtrs_(&trans, &n, &kl, &ku, &nrhs, ab, &ldab, ipiv, b, &ldb, &info, 1);
}

/// Transforms eigenvectors of a balanced matrix to those of the
/// original matrix supplied to SGEBAL.
void gebak(char job, char side, f_int n, f_int ilo, f_int ihi, f_float *scale, f_int m, f_float *v, f_int ldv, ref f_int info) {
    sgebak_(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info, 1, 1);
}
void gebak(char job, char side, f_int n, f_int ilo, f_int ihi, f_double *scale, f_int m, f_double *v, f_int ldv, ref f_int info) {
    dgebak_(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info, 1, 1);
}
void gebak(char job, char side, f_int n, f_int ilo, f_int ihi, f_float *scale, f_int m, f_cfloat *v, f_int ldv, ref f_int info) {
    cgebak_(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info, 1, 1);
}
void gebak(char job, char side, f_int n, f_int ilo, f_int ihi, f_double *scale, f_int m, f_cdouble *v, f_int ldv, ref f_int info) {
    zgebak_(&job, &side, &n, &ilo, &ihi, scale, &m, v, &ldv, &info, 1, 1);
}

/// Balances a general matrix in order to improve the accuracy
/// of computed eigenvalues.
void gebal(char job, f_int n, f_float *a, f_int lda, f_int ilo, f_int ihi, f_float *scale, ref f_int info) {
    sgebal_(&job, &n, a, &lda, &ilo, &ihi, scale, &info, 1);
}
void gebal(char job, f_int n, f_double *a, f_int lda, f_int ilo, f_int ihi, f_double *scale, ref f_int info) {
    dgebal_(&job, &n, a, &lda, &ilo, &ihi, scale, &info, 1);
}
void gebal(char job, f_int n, f_cfloat *a, f_int lda, f_int ilo, f_int ihi, f_float *scale, ref f_int info) {
    cgebal_(&job, &n, a, &lda, &ilo, &ihi, scale, &info, 1);
}
void gebal(char job, f_int n, f_cdouble *a, f_int lda, f_int ilo, f_int ihi, f_double *scale, ref f_int info) {
    zgebal_(&job, &n, a, &lda, &ilo, &ihi, scale, &info, 1);
}

/// Reduces a general rectangular matrix to real bidiagonal form
/// by an orthogonal transformation.
void gebrd(f_int m, f_int n, f_float *a, f_int lda, f_float *d, f_float *e, f_float *tauq, f_float *taup, f_float *work, f_int lwork, ref f_int info) {
    sgebrd_(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
}
void gebrd(f_int m, f_int n, f_double *a, f_int lda, f_double *d, f_double *e, f_double *tauq, f_double *taup, f_double *work, f_int lwork, ref f_int info) {
    dgebrd_(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
}
void gebrd(f_int m, f_int n, f_cfloat *a, f_int lda, f_float *d, f_float *e, f_cfloat *tauq, f_cfloat *taup, f_cfloat *work, f_int lwork, ref f_int info) {
    cgebrd_(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
}
void gebrd(f_int m, f_int n, f_cdouble *a, f_int lda, f_double *d, f_double *e, f_cdouble *tauq, f_cdouble *taup, f_cdouble *work, f_int lwork, ref f_int info) {
    zgebrd_(&m, &n, a, &lda, d, e, tauq, taup, work, &lwork, &info);
}

/// Estimates the reciprocal of the condition number of a general
/// matrix, in either the 1-norm or the infinity-norm, using the
/// LU factorization computed by SGETRF.
void gecon(char norm, f_int n, f_float *a, f_int lda, f_float *anorm, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    sgecon_(&norm, &n, a, &lda, anorm, &rcond, work, iwork, &info, 1);
}
void gecon(char norm, f_int n, f_double *a, f_int lda, f_double *anorm, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dgecon_(&norm, &n, a, &lda, anorm, &rcond, work, iwork, &info, 1);
}
void gecon(char norm, f_int n, f_cfloat *a, f_int lda, f_float *anorm, f_float rcond, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgecon_(&norm, &n, a, &lda, anorm, &rcond, work, rwork, &info, 1);
}
void gecon(char norm, f_int n, f_cdouble *a, f_int lda, f_double *anorm, f_double rcond, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgecon_(&norm, &n, a, &lda, anorm, &rcond, work, rwork, &info, 1);
}

/// Computes row and column scalings to equilibrate a general
/// rectangular matrix and reduce its condition number.
void geequ(f_int m, f_int n, f_float *a, f_int lda, f_float *r, f_float *c, f_float *rowcnd, f_float *colcnd, f_float *amax, ref f_int info) {
    sgeequ_(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
}
void geequ(f_int m, f_int n, f_double *a, f_int lda, f_double *r, f_double *c, f_double *rowcnd, f_double *colcnd, f_double *amax, ref f_int info) {
    dgeequ_(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
}
void geequ(f_int m, f_int n, f_cfloat *a, f_int lda, f_float *r, f_float *c, f_float *rowcnd, f_float *colcnd, f_float *amax, ref f_int info) {
    cgeequ_(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
}
void geequ(f_int m, f_int n, f_cdouble *a, f_int lda, f_double *r, f_double *c, f_double *rowcnd, f_double *colcnd, f_double *amax, ref f_int info) {
    zgeequ_(&m, &n, a, &lda, r, c, rowcnd, colcnd, amax, &info);
}

/// Reduces a general matrix to upper Hessenberg form by an
/// orthogonal similarity transformation.
void gehrd(f_int n, f_int ilo, f_int ihi, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sgehrd_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
}
void gehrd(f_int n, f_int ilo, f_int ihi, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dgehrd_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
}
void gehrd(f_int n, f_int ilo, f_int ihi, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cgehrd_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
}
void gehrd(f_int n, f_int ilo, f_int ihi, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zgehrd_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
}

/// Computes an LQ factorization of a general rectangular matrix.
void gelqf(f_int m, f_int n, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sgelqf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void gelqf(f_int m, f_int n, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dgelqf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void gelqf(f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cgelqf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void gelqf(f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zgelqf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}

/// Computes a QL factorization of a general rectangular matrix.
void geqlf(f_int m, f_int n, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sgeqlf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void geqlf(f_int m, f_int n, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dgeqlf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void geqlf(f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cgeqlf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void geqlf(f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zgeqlf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}

/// Computes a QR factorization with column pivoting of a general
/// rectangular matrix using Level 3 BLAS.
void geqp3(f_int m, f_int n, f_float *a, f_int lda, f_int jpvt, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sgeqp3_(&m, &n, a, &lda, &jpvt, tau, work, &lwork, &info);
}
void geqp3(f_int m, f_int n, f_double *a, f_int lda, f_int jpvt, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dgeqp3_(&m, &n, a, &lda, &jpvt, tau, work, &lwork, &info);
}
void geqp3(f_int m, f_int n, f_cfloat *a, f_int lda, f_int jpvt, f_cfloat *tau, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    cgeqp3_(&m, &n, a, &lda, &jpvt, tau, work, &lwork, rwork, &info);
}
void geqp3(f_int m, f_int n, f_cdouble *a, f_int lda, f_int jpvt, f_cdouble *tau, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zgeqp3_(&m, &n, a, &lda, &jpvt, tau, work, &lwork, rwork, &info);
}

/// Computes a QR factorization with column pivoting of a general
/// rectangular matrix.
void geqpf(f_int m, f_int n, f_float *a, f_int lda, f_int jpvt, f_float *tau, f_float *work, ref f_int info) {
    sgeqpf_(&m, &n, a, &lda, &jpvt, tau, work, &info);
}
void geqpf(f_int m, f_int n, f_double *a, f_int lda, f_int jpvt, f_double *tau, f_double *work, ref f_int info) {
    dgeqpf_(&m, &n, a, &lda, &jpvt, tau, work, &info);
}
void geqpf(f_int m, f_int n, f_cfloat *a, f_int lda, f_int jpvt, f_cfloat *tau, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgeqpf_(&m, &n, a, &lda, &jpvt, tau, work, rwork, &info);
}
void geqpf(f_int m, f_int n, f_cdouble *a, f_int lda, f_int jpvt, f_cdouble *tau, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgeqpf_(&m, &n, a, &lda, &jpvt, tau, work, rwork, &info);
}

/// Computes a QR factorization of a general rectangular matrix.
void geqrf(f_int m, f_int n, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void geqrf(f_int m, f_int n, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void geqrf(f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void geqrf(f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zgeqrf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}

/// Improves the computed solution to a general system of linear
/// equations AX=B, A**T X=B or A**H X=B, and provides forward and
/// backward error bounds for the solution.
void gerfs(char trans, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *af, f_int ldaf, f_int *ipiv, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sgerfs_(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void gerfs(char trans, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *af, f_int ldaf, f_int *ipiv, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dgerfs_(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void gerfs(char trans, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *af, f_int ldaf, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgerfs_(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void gerfs(char trans, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *af, f_int ldaf, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgerfs_(&trans, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Computes an RQ factorization of a general rectangular matrix.
void gerqf(f_int m, f_int n, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sgerqf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void gerqf(f_int m, f_int n, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dgerqf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void gerqf(f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cgerqf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void gerqf(f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zgerqf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}

/// Computes an LU factorization of a general matrix, using partial
/// pivoting with row interchanges.
void getrf(f_int m, f_int n, f_float *a, f_int lda, f_int *ipiv, ref f_int info) {
    sgetrf_(&m, &n, a, &lda, ipiv, &info);
}
void getrf(f_int m, f_int n, f_double *a, f_int lda, f_int *ipiv, ref f_int info) {
    dgetrf_(&m, &n, a, &lda, ipiv, &info);
}
void getrf(f_int m, f_int n, f_cfloat *a, f_int lda, f_int *ipiv, ref f_int info) {
    cgetrf_(&m, &n, a, &lda, ipiv, &info);
}
void getrf(f_int m, f_int n, f_cdouble *a, f_int lda, f_int *ipiv, ref f_int info) {
    zgetrf_(&m, &n, a, &lda, ipiv, &info);
}

/// Computes the inverse of a general matrix, using the LU factorization
/// computed by SGETRF.
void getri(f_int n, f_float *a, f_int lda, f_int *ipiv, f_float *work, f_int lwork, ref f_int info) {
    sgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
}
void getri(f_int n, f_double *a, f_int lda, f_int *ipiv, f_double *work, f_int lwork, ref f_int info) {
    dgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
}
void getri(f_int n, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *work, f_int lwork, ref f_int info) {
    cgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
}
void getri(f_int n, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *work, f_int lwork, ref f_int info) {
    zgetri_(&n, a, &lda, ipiv, work, &lwork, &info);
}

/// Solves a general system of linear equations AX=B, A**T X=B
/// or A**H X=B, using the LU factorization computed by SGETRF.
void getrs(char trans, f_int n, f_int nrhs, f_float *a, f_int lda, f_int *ipiv, f_float *b, f_int ldb, ref f_int info) {
    sgetrs_(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}
void getrs(char trans, f_int n, f_int nrhs, f_double *a, f_int lda, f_int *ipiv, f_double *b, f_int ldb, ref f_int info) {
    dgetrs_(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}
void getrs(char trans, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    cgetrs_(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}
void getrs(char trans, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zgetrs_(&trans, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}

/// Forms the right or left eigenvectors of the generalized eigenvalue
/// problem by backward transformation on the computed eigenvectors of
/// the balanced pair of matrices output by SGGBAL.
void ggbak(char job, char side, f_int n, f_int ilo, f_int ihi, f_float *lscale, f_float *rscale, f_int m, f_float *v, f_int ldv, ref f_int info) {
    sggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info, 1, 1);
}
void ggbak(char job, char side, f_int n, f_int ilo, f_int ihi, f_double *lscale, f_double *rscale, f_int m, f_double *v, f_int ldv, ref f_int info) {
    dggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info, 1, 1);
}
void ggbak(char job, char side, f_int n, f_int ilo, f_int ihi, f_float *lscale, f_float *rscale, f_int m, f_cfloat *v, f_int ldv, ref f_int info) {
    cggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info, 1, 1);
}
void ggbak(char job, char side, f_int n, f_int ilo, f_int ihi, f_double *lscale, f_double *rscale, f_int m, f_cdouble *v, f_int ldv, ref f_int info) {
    zggbak_(&job, &side, &n, &ilo, &ihi, lscale, rscale, &m, v, &ldv, &info, 1, 1);
}

/// Balances a pair of general real matrices for the generalized
/// eigenvalue problem A x = lambda B x.
void ggbal(char job, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_int ilo, f_int ihi, f_float *lscale, f_float *rscale, f_float *work, ref f_int info) {
    sggbal_(&job, &n, a, &lda, b, &ldb, &ilo, &ihi, lscale, rscale, work, &info, 1);
}
void ggbal(char job, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_int ilo, f_int ihi, f_double *lscale, f_double *rscale, f_double *work, ref f_int info) {
    dggbal_(&job, &n, a, &lda, b, &ldb, &ilo, &ihi, lscale, rscale, work, &info, 1);
}
void ggbal(char job, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_int ilo, f_int ihi, f_float *lscale, f_float *rscale, f_float *work, ref f_int info) {
    cggbal_(&job, &n, a, &lda, b, &ldb, &ilo, &ihi, lscale, rscale, work, &info, 1);
}
void ggbal(char job, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_int ilo, f_int ihi, f_double *lscale, f_double *rscale, f_double *work, ref f_int info) {
    zggbal_(&job, &n, a, &lda, b, &ldb, &ilo, &ihi, lscale, rscale, work, &info, 1);
}

/// Reduces a pair of real matrices to generalized upper
/// Hessenberg form using orthogonal transformations 
void gghrd(char compq, char compz, f_int n, f_int ilo, f_int ihi, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *q, f_int ldq, f_float *z, f_int ldz, ref f_int info) {
    sgghrd_(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info, 1, 1);
}
void gghrd(char compq, char compz, f_int n, f_int ilo, f_int ihi, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *q, f_int ldq, f_double *z, f_int ldz, ref f_int info) {
    dgghrd_(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info, 1, 1);
}
void gghrd(char compq, char compz, f_int n, f_int ilo, f_int ihi, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *q, f_int ldq, f_cfloat *z, f_int ldz, ref f_int info) {
    cgghrd_(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info, 1, 1);
}
void gghrd(char compq, char compz, f_int n, f_int ilo, f_int ihi, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *q, f_int ldq, f_cdouble *z, f_int ldz, ref f_int info) {
    zgghrd_(&compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, q, &ldq, z, &ldz, &info, 1, 1);
}

/// Computes a generalized QR factorization of a pair of matrices. 
void ggqrf(f_int n, f_int m, f_int p, f_float *a, f_int lda, f_float *taua, f_float *b, f_int ldb, f_float *taub, f_float *work, f_int lwork, ref f_int info) {
    sggqrf_(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
}
void ggqrf(f_int n, f_int m, f_int p, f_double *a, f_int lda, f_double *taua, f_double *b, f_int ldb, f_double *taub, f_double *work, f_int lwork, ref f_int info) {
    dggqrf_(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
}
void ggqrf(f_int n, f_int m, f_int p, f_cfloat *a, f_int lda, f_cfloat *taua, f_cfloat *b, f_int ldb, f_cfloat *taub, f_cfloat *work, f_int lwork, ref f_int info) {
    cggqrf_(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
}
void ggqrf(f_int n, f_int m, f_int p, f_cdouble *a, f_int lda, f_cdouble *taua, f_cdouble *b, f_int ldb, f_cdouble *taub, f_cdouble *work, f_int lwork, ref f_int info) {
    zggqrf_(&n, &m, &p, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
}

/// Computes a generalized RQ factorization of a pair of matrices.
void ggrqf(f_int m, f_int p, f_int n, f_float *a, f_int lda, f_float *taua, f_float *b, f_int ldb, f_float *taub, f_float *work, f_int lwork, ref f_int info) {
    sggrqf_(&m, &p, &n, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
}
void ggrqf(f_int m, f_int p, f_int n, f_double *a, f_int lda, f_double *taua, f_double *b, f_int ldb, f_double *taub, f_double *work, f_int lwork, ref f_int info) {
    dggrqf_(&m, &p, &n, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
}
void ggrqf(f_int m, f_int p, f_int n, f_cfloat *a, f_int lda, f_cfloat *taua, f_cfloat *b, f_int ldb, f_cfloat *taub, f_cfloat *work, f_int lwork, ref f_int info) {
    cggrqf_(&m, &p, &n, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
}
void ggrqf(f_int m, f_int p, f_int n, f_cdouble *a, f_int lda, f_cdouble *taua, f_cdouble *b, f_int ldb, f_cdouble *taub, f_cdouble *work, f_int lwork, ref f_int info) {
    zggrqf_(&m, &p, &n, a, &lda, taua, b, &ldb, taub, work, &lwork, &info);
}

/// Computes orthogonal matrices as a preprocessing step
/// for computing the generalized singular value decomposition
void ggsvp(char jobu, char jobv, char jobq, f_int m, f_int p, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *tola, f_float *tolb, f_int k, f_int l, f_float *u, f_int ldu, f_float *v, f_int ldv, f_float *q, f_int ldq, f_int *iwork, f_float *tau, f_float *work, ref f_int info) {
    sggsvp_(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, tola, tolb, &k, &l, u, &ldu, v, &ldv, q, &ldq, iwork, tau, work, &info, 1, 1, 1);
}
void ggsvp(char jobu, char jobv, char jobq, f_int m, f_int p, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *tola, f_double *tolb, f_int k, f_int l, f_double *u, f_int ldu, f_double *v, f_int ldv, f_double *q, f_int ldq, f_int *iwork, f_double *tau, f_double *work, ref f_int info) {
    dggsvp_(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, tola, tolb, &k, &l, u, &ldu, v, &ldv, q, &ldq, iwork, tau, work, &info, 1, 1, 1);
}
void ggsvp(char jobu, char jobv, char jobq, f_int m, f_int p, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_float *tola, f_float *tolb, f_int k, f_int l, f_cfloat *u, f_int ldu, f_cfloat *v, f_int ldv, f_cfloat *q, f_int ldq, f_int *iwork, f_float *rwork, f_cfloat *tau, f_cfloat *work, ref f_int info) {
    cggsvp_(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, tola, tolb, &k, &l, u, &ldu, v, &ldv, q, &ldq, iwork, rwork, tau, work, &info, 1, 1, 1);
}
void ggsvp(char jobu, char jobv, char jobq, f_int m, f_int p, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_double *tola, f_double *tolb, f_int k, f_int l, f_cdouble *u, f_int ldu, f_cdouble *v, f_int ldv, f_cdouble *q, f_int ldq, f_int *iwork, f_double *rwork, f_cdouble *tau, f_cdouble *work, ref f_int info) {
    zggsvp_(&jobu, &jobv, &jobq, &m, &p, &n, a, &lda, b, &ldb, tola, tolb, &k, &l, u, &ldu, v, &ldv, q, &ldq, iwork, rwork, tau, work, &info, 1, 1, 1);
}

/// Estimates the reciprocal of the condition number of a general
/// tridiagonal matrix, in either the 1-norm or the infinity-norm,
/// using the LU factorization computed by SGTTRF.
void gtcon(char norm, f_int n, f_float *dl, f_float *d, f_float *du, f_float *du2, f_int *ipiv, f_float *anorm, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    sgtcon_(&norm, &n, dl, d, du, du2, ipiv, anorm, &rcond, work, iwork, &info, 1);
}
void gtcon(char norm, f_int n, f_double *dl, f_double *d, f_double *du, f_double *du2, f_int *ipiv, f_double *anorm, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dgtcon_(&norm, &n, dl, d, du, du2, ipiv, anorm, &rcond, work, iwork, &info, 1);
}
void gtcon(char norm, f_int n, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *du2, f_int *ipiv, f_float *anorm, f_float rcond, f_cfloat *work, ref f_int info) {
    cgtcon_(&norm, &n, dl, d, du, du2, ipiv, anorm, &rcond, work, &info, 1);
}
void gtcon(char norm, f_int n, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *du2, f_int *ipiv, f_double *anorm, f_double rcond, f_cdouble *work, ref f_int info) {
    zgtcon_(&norm, &n, dl, d, du, du2, ipiv, anorm, &rcond, work, &info, 1);
}

/// Improves the computed solution to a general tridiagonal system
/// of linear equations AX=B, A**T X=B or A**H X=B, and provides
/// forward and backward error bounds for the solution.
void gtrfs(char trans, f_int n, f_int nrhs, f_float *dl, f_float *d, f_float *du, f_float *dlf, f_float *df, f_float *duf, f_float *du2, f_int *ipiv, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sgtrfs_(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void gtrfs(char trans, f_int n, f_int nrhs, f_double *dl, f_double *d, f_double *du, f_double *dlf, f_double *df, f_double *duf, f_double *du2, f_int *ipiv, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dgtrfs_(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void gtrfs(char trans, f_int n, f_int nrhs, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *dlf, f_cfloat *df, f_cfloat *duf, f_cfloat *du2, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cgtrfs_(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void gtrfs(char trans, f_int n, f_int nrhs, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *dlf, f_cdouble *df, f_cdouble *duf, f_cdouble *du2, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zgtrfs_(&trans, &n, &nrhs, dl, d, du, dlf, df, duf, du2, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Computes an LU factorization of a general tridiagonal matrix,
/// using partial pivoting with row interchanges.
void gttrf(f_int n, f_float *dl, f_float *d, f_float *du, f_float *du2, f_int *ipiv, ref f_int info) {
    sgttrf_(&n, dl, d, du, du2, ipiv, &info);
}
void gttrf(f_int n, f_double *dl, f_double *d, f_double *du, f_double *du2, f_int *ipiv, ref f_int info) {
    dgttrf_(&n, dl, d, du, du2, ipiv, &info);
}
void gttrf(f_int n, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *du2, f_int *ipiv, ref f_int info) {
    cgttrf_(&n, dl, d, du, du2, ipiv, &info);
}
void gttrf(f_int n, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *du2, f_int *ipiv, ref f_int info) {
    zgttrf_(&n, dl, d, du, du2, ipiv, &info);
}

/// Solves a general tridiagonal system of linear equations AX=B,
/// A**T X=B or A**H X=B, using the LU factorization computed by
/// SGTTRF.
void gttrs(char trans, f_int n, f_int nrhs, f_float *dl, f_float *d, f_float *du, f_float *du2, f_int *ipiv, f_float *b, f_int ldb, ref f_int info) {
    sgttrs_(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info, 1);
}
void gttrs(char trans, f_int n, f_int nrhs, f_double *dl, f_double *d, f_double *du, f_double *du2, f_int *ipiv, f_double *b, f_int ldb, ref f_int info) {
    dgttrs_(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info, 1);
}
void gttrs(char trans, f_int n, f_int nrhs, f_cfloat *dl, f_cfloat *d, f_cfloat *du, f_cfloat *du2, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    cgttrs_(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info, 1);
}
void gttrs(char trans, f_int n, f_int nrhs, f_cdouble *dl, f_cdouble *d, f_cdouble *du, f_cdouble *du2, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zgttrs_(&trans, &n, &nrhs, dl, d, du, du2, ipiv, b, &ldb, &info, 1);
}

/// Implements a single-/f_double-shift version of the QZ method for
/// finding the generalized eigenvalues of the equation 
/// det(A - w(i) B) = 0
void hgeqz(char job, char compq, char compz, f_int n, f_int ilo, f_int ihi, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *q, f_int ldq, f_float *z, f_int ldz, f_float *work, f_int lwork, ref f_int info) {
    shgeqz_(&job, &compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, alphar, alphai, betav, q, &ldq, z, &ldz, work, &lwork, &info, 1, 1, 1);
}
void hgeqz(char job, char compq, char compz, f_int n, f_int ilo, f_int ihi, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *q, f_int ldq, f_double *z, f_int ldz, f_double *work, f_int lwork, ref f_int info) {
    dhgeqz_(&job, &compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, alphar, alphai, betav, q, &ldq, z, &ldz, work, &lwork, &info, 1, 1, 1);
}
void hgeqz(char job, char compq, char compz, f_int n, f_int ilo, f_int ihi, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *q, f_int ldq, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, f_float *rwork, ref f_int info) {
    chgeqz_(&job, &compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, alphav, betav, q, &ldq, z, &ldz, work, &lwork, rwork, &info, 1, 1, 1);
}
void hgeqz(char job, char compq, char compz, f_int n, f_int ilo, f_int ihi, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *q, f_int ldq, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, f_double *rwork, ref f_int info) {
    zhgeqz_(&job, &compq, &compz, &n, &ilo, &ihi, a, &lda, b, &ldb, alphav, betav, q, &ldq, z, &ldz, work, &lwork, rwork, &info, 1, 1, 1);
}

/// Computes specified right and/or left eigenvectors of an upper
/// Hessenberg matrix by inverse iteration.
void hsein(char side, char eigsrc, char initv, f_int select, f_int n, f_float *h, f_int ldh, f_float *wr, f_float *wi, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_int mm, f_int m, f_float *work, f_int ifaill, f_int ifailr, ref f_int info) {
    shsein_(&side, &eigsrc, &initv, &select, &n, h, &ldh, wr, wi, vl, &ldvl, vr, &ldvr, &mm, &m, work, &ifaill, &ifailr, &info, 1, 1, 1);
}
void hsein(char side, char eigsrc, char initv, f_int select, f_int n, f_double *h, f_int ldh, f_double *wr, f_double *wi, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_int mm, f_int m, f_double *work, f_int ifaill, f_int ifailr, ref f_int info) {
    dhsein_(&side, &eigsrc, &initv, &select, &n, h, &ldh, wr, wi, vl, &ldvl, vr, &ldvr, &mm, &m, work, &ifaill, &ifailr, &info, 1, 1, 1);
}
void hsein(char side, char eigsrc, char initv, f_int select, f_int n, f_cfloat *h, f_int ldh, f_cfloat *w, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_int mm, f_int m, f_cfloat *work, f_float *rwork, f_int ifaill, f_int ifailr, ref f_int info) {
    chsein_(&side, &eigsrc, &initv, &select, &n, h, &ldh, w, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &ifaill, &ifailr, &info, 1, 1, 1);
}
void hsein(char side, char eigsrc, char initv, f_int select, f_int n, f_cdouble *h, f_int ldh, f_cdouble *w, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_int mm, f_int m, f_cdouble *work, f_double *rwork, f_int ifaill, f_int ifailr, ref f_int info) {
    zhsein_(&side, &eigsrc, &initv, &select, &n, h, &ldh, w, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &ifaill, &ifailr, &info, 1, 1, 1);
}

/// Computes the eigenvalues and Schur factorization of an upper
/// Hessenberg matrix, using the multishift QR algorithm.
void hseqr(char job, char compz, f_int n, f_int ilo, f_int ihi, f_float *h, f_int ldh, f_float *wr, f_float *wi, f_float *z, f_int ldz, f_float *work, f_int lwork, ref f_int info) {
    shseqr_(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, work, &lwork, &info, 1, 1);
}
void hseqr(char job, char compz, f_int n, f_int ilo, f_int ihi, f_double *h, f_int ldh, f_double *wr, f_double *wi, f_double *z, f_int ldz, f_double *work, f_int lwork, ref f_int info) {
    dhseqr_(&job, &compz, &n, &ilo, &ihi, h, &ldh, wr, wi, z, &ldz, work, &lwork, &info, 1, 1);
}
void hseqr(char job, char compz, f_int n, f_int ilo, f_int ihi, f_cfloat *h, f_int ldh, f_cfloat *w, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, ref f_int info) {
    chseqr_(&job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, work, &lwork, &info, 1, 1);
}
void hseqr(char job, char compz, f_int n, f_int ilo, f_int ihi, f_cdouble *h, f_int ldh, f_cdouble *w, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, ref f_int info) {
    zhseqr_(&job, &compz, &n, &ilo, &ihi, h, &ldh, w, z, &ldz, work, &lwork, &info, 1, 1);
}

/// Generates the orthogonal transformation matrix from
/// a reduction to tridiagonal form determined by SSPTRD.
void opgtr(char uplo, f_int n, f_float *ap, f_float *tau, f_float *q, f_int ldq, f_float *work, ref f_int info) {
    sopgtr_(&uplo, &n, ap, tau, q, &ldq, work, &info, 1);
}
void opgtr(char uplo, f_int n, f_double *ap, f_double *tau, f_double *q, f_int ldq, f_double *work, ref f_int info) {
    dopgtr_(&uplo, &n, ap, tau, q, &ldq, work, &info, 1);
}

/// Generates the unitary transformation matrix from
/// a reduction to tridiagonal form determined by CHPTRD.
void upgtr(char uplo, f_int n, f_cfloat *ap, f_cfloat *tau, f_cfloat *q, f_int ldq, f_cfloat *work, ref f_int info) {
    cupgtr_(&uplo, &n, ap, tau, q, &ldq, work, &info, 1);
}
void upgtr(char uplo, f_int n, f_cdouble *ap, f_cdouble *tau, f_cdouble *q, f_int ldq, f_cdouble *work, ref f_int info) {
    zupgtr_(&uplo, &n, ap, tau, q, &ldq, work, &info, 1);
}


/// Multiplies a general matrix by the orthogonal
/// transformation matrix from a reduction to tridiagonal form
/// determined by SSPTRD.
void opmtr(char side, char uplo, char trans, f_int m, f_int n, f_float *ap, f_float *tau, f_float *c, f_int ldc, f_float *work, ref f_int info) {
    sopmtr_(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info, 1, 1, 1);
}
void opmtr(char side, char uplo, char trans, f_int m, f_int n, f_double *ap, f_double *tau, f_double *c, f_int ldc, f_double *work, ref f_int info) {
    dopmtr_(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info, 1, 1, 1);
}

/// Generates the orthogonal transformation matrices from
/// a reduction to bidiagonal form determined by SGEBRD.
void orgbr(char vect, f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sorgbr_(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info, 1);
}
void orgbr(char vect, f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dorgbr_(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info, 1);
}

/// Generates the unitary transformation matrices from
/// a reduction to bidiagonal form determined by CGEBRD.
void ungbr(char vect, f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cungbr_(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info, 1);
}
void ungbr(char vect, f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zungbr_(&vect, &m, &n, &k, a, &lda, tau, work, &lwork, &info, 1);
}

/// Generates the orthogonal transformation matrix from
/// a reduction to Hessenberg form determined by SGEHRD.
void orghr(f_int n, f_int ilo, f_int ihi, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sorghr_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
}
void orghr(f_int n, f_int ilo, f_int ihi, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dorghr_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
}

/// Generates the unitary transformation matrix from
/// a reduction to Hessenberg form determined by CGEHRD.
void unghr(f_int n, f_int ilo, f_int ihi, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cunghr_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
}
void unghr(f_int n, f_int ilo, f_int ihi, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zunghr_(&n, &ilo, &ihi, a, &lda, tau, work, &lwork, &info);
}

/// Generates all or part of the orthogonal matrix Q from
/// an LQ factorization determined by SGELQF.
void orglq(f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sorglq_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}
void orglq(f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dorglq_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

/// Generates all or part of the unitary matrix Q from
/// an LQ factorization determined by CGELQF.
void unglq(f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cunglq_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}
void unglq(f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zunglq_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

/// Generates all or part of the orthogonal matrix Q from
/// a QL factorization determined by SGEQLF.
void orgql(f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sorgql_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}
void orgql(f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dorgql_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

/// Generates all or part of the unitary matrix Q from
/// a QL factorization determined by CGEQLF.
void ungql(f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cungql_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}
void ungql(f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zungql_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

/// Generates all or part of the orthogonal matrix Q from
/// a QR factorization determined by SGEQRF.
void orgqr(f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}
void orgqr(f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dorgqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

/// Generates all or part of the unitary matrix Q from
/// a QR factorization determined by CGEQRF.
void ungqr(f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cungqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}
void ungqr(f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zungqr_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

/// Generates all or part of the orthogonal matrix Q from
/// an RQ factorization determined by SGERQF.
void orgrq(f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sorgrq_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}
void orgrq(f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dorgrq_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

/// Generates all or part of the unitary matrix Q from
/// an RQ factorization determined by CGERQF.
void ungrq(f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cungrq_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}
void ungrq(f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zungrq_(&m, &n, &k, a, &lda, tau, work, &lwork, &info);
}

/// Generates the orthogonal transformation matrix from
/// a reduction to tridiagonal form determined by SSYTRD.
void orgtr(char uplo, f_int n, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    sorgtr_(&uplo, &n, a, &lda, tau, work, &lwork, &info, 1);
}
void orgtr(char uplo, f_int n, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dorgtr_(&uplo, &n, a, &lda, tau, work, &lwork, &info, 1);
}

/// Generates the unitary transformation matrix from
/// a reduction to tridiagonal form determined by CHETRD.
void ungtr(char uplo, f_int n, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    cungtr_(&uplo, &n, a, &lda, tau, work, &lwork, &info, 1);
}
void ungtr(char uplo, f_int n, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zungtr_(&uplo, &n, a, &lda, tau, work, &lwork, &info, 1);
}

/// Multiplies a general matrix by one of the orthogonal
/// transformation  matrices from a reduction to bidiagonal form
/// determined by SGEBRD.
void ormbr(char vect, char side, char trans, f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, f_int lwork, ref f_int info) {
    sormbr_(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1, 1);
}
void ormbr(char vect, char side, char trans, f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, f_int lwork, ref f_int info) {
    dormbr_(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1, 1);
}

/// Multiplies a general matrix by one of the unitary
/// transformation matrices from a reduction to bidiagonal form
/// determined by CGEBRD.
void unmbr(char vect, char side, char trans, f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, f_int lwork, ref f_int info) {
    cunmbr_(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1, 1);
}
void unmbr(char vect, char side, char trans, f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, f_int lwork, ref f_int info) {
    zunmbr_(&vect, &side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1, 1);
}

/// Multiplies a general matrix by the orthogonal transformation
/// matrix from a reduction to Hessenberg form determined by SGEHRD.
void ormhr(char side, char trans, f_int m, f_int n, f_int ilo, f_int ihi, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, f_int lwork, ref f_int info) {
    sormhr_(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void ormhr(char side, char trans, f_int m, f_int n, f_int ilo, f_int ihi, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, f_int lwork, ref f_int info) {
    dormhr_(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the unitary transformation
/// matrix from a reduction to Hessenberg form determined by CGEHRD.
void unmhr(char side, char trans, f_int m, f_int n, f_int ilo, f_int ihi, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, f_int lwork, ref f_int info) {
    cunmhr_(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void unmhr(char side, char trans, f_int m, f_int n, f_int ilo, f_int ihi, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, f_int lwork, ref f_int info) {
    zunmhr_(&side, &trans, &m, &n, &ilo, &ihi, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the orthogonal matrix
/// from an LQ factorization determined by SGELQF.
void ormlq(char side, char trans, f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, f_int lwork, ref f_int info) {
    sormlq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void ormlq(char side, char trans, f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, f_int lwork, ref f_int info) {
    dormlq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the unitary matrix
/// from an LQ factorization determined by CGELQF.
void unmlq(char side, char trans, f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, f_int lwork, ref f_int info) {
    cunmlq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void unmlq(char side, char trans, f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, f_int lwork, ref f_int info) {
    zunmlq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the orthogonal matrix
/// from a QL factorization determined by SGEQLF.
void ormql(char side, char trans, f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, f_int lwork, ref f_int info) {
    sormql_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void ormql(char side, char trans, f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, f_int lwork, ref f_int info) {
    dormql_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the unitary matrix
/// from a QL factorization determined by CGEQLF.
void unmql(char side, char trans, f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, f_int lwork, ref f_int info) {
    cunmql_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void unmql(char side, char trans, f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, f_int lwork, ref f_int info) {
    zunmql_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the orthogonal matrix
/// from a QR factorization determined by SGEQRF.
void ormqr(char side, char trans, f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, f_int lwork, ref f_int info) {
    sormqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void ormqr(char side, char trans, f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, f_int lwork, ref f_int info) {
    dormqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the unitary matrix
/// from a QR factorization determined by CGEQRF.
void unmqr(char side, char trans, f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, f_int lwork, ref f_int info) {
    cunmqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void unmqr(char side, char trans, f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, f_int lwork, ref f_int info) {
    zunmqr_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiples a general matrix by the orthogonal matrix
/// from an RZ factorization determined by STZRZF.
void ormr3(char side, char trans, f_int m, f_int n, f_int k, f_int l, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, ref f_int info) {
    sormr3_(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &info, 1, 1);
}
void ormr3(char side, char trans, f_int m, f_int n, f_int k, f_int l, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, ref f_int info) {
    dormr3_(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &info, 1, 1);
}

/// Multiples a general matrix by the unitary matrix
/// from an RZ factorization determined by CTZRZF.
void unmr3(char side, char trans, f_int m, f_int n, f_int k, f_int l, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, ref f_int info) {
    cunmr3_(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &info, 1, 1);
}
void unmr3(char side, char trans, f_int m, f_int n, f_int k, f_int l, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, ref f_int info) {
    zunmr3_(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &info, 1, 1);
}

/// Multiplies a general matrix by the orthogonal matrix
/// from an RQ factorization determined by SGERQF.
void ormrq(char side, char trans, f_int m, f_int n, f_int k, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, f_int lwork, ref f_int info) {
    sormrq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void ormrq(char side, char trans, f_int m, f_int n, f_int k, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, f_int lwork, ref f_int info) {
    dormrq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the unitary matrix
/// from an RQ factorization determined by CGERQF.
void unmrq(char side, char trans, f_int m, f_int n, f_int k, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, f_int lwork, ref f_int info) {
    cunmrq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void unmrq(char side, char trans, f_int m, f_int n, f_int k, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, f_int lwork, ref f_int info) {
    zunmrq_(&side, &trans, &m, &n, &k, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiples a general matrix by the orthogonal matrix
/// from an RZ factorization determined by STZRZF.
void ormrz(char side, char trans, f_int m, f_int n, f_int k, f_int l, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, f_int lwork, ref f_int info) {
    sormrz_(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void ormrz(char side, char trans, f_int m, f_int n, f_int k, f_int l, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, f_int lwork, ref f_int info) {
    dormrz_(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiples a general matrix by the unitary matrix
/// from an RZ factorization determined by CTZRZF.
void unmrz(char side, char trans, f_int m, f_int n, f_int k, f_int l, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, f_int lwork, ref f_int info) {
    cunmrz_(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}
void unmrz(char side, char trans, f_int m, f_int n, f_int k, f_int l, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, f_int lwork, ref f_int info) {
    zunmrz_(&side, &trans, &m, &n, &k, &l, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1);
}

/// Multiplies a general matrix by the orthogonal
/// transformation matrix from a reduction to tridiagonal form
/// determined by SSYTRD.
void ormtr(char side, char uplo, char trans, f_int m, f_int n, f_float *a, f_int lda, f_float *tau, f_float *c, f_int ldc, f_float *work, f_int lwork, ref f_int info) {
    sormtr_(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1, 1);
}
void ormtr(char side, char uplo, char trans, f_int m, f_int n, f_double *a, f_int lda, f_double *tau, f_double *c, f_int ldc, f_double *work, f_int lwork, ref f_int info) {
    dormtr_(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1, 1);
}

/// Multiplies a general matrix by the unitary
/// transformation matrix from a reduction to tridiagonal form
/// determined by CHETRD.
void unmtr(char side, char uplo, char trans, f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, f_int lwork, ref f_int info) {
    cunmtr_(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1, 1);
}
void unmtr(char side, char uplo, char trans, f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, f_int lwork, ref f_int info) {
    zunmtr_(&side, &uplo, &trans, &m, &n, a, &lda, tau, c, &ldc, work, &lwork, &info, 1, 1, 1);
}

/// Estimates the reciprocal of the condition number of a
/// symmetric positive definite band matrix, using the
/// Cholesky factorization computed by SPBTRF.
void pbcon(char uplo, f_int n, f_int kd, f_float *ab, f_int ldab, f_float *anorm, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    spbcon_(&uplo, &n, &kd, ab, &ldab, anorm, &rcond, work, iwork, &info, 1);
}
void pbcon(char uplo, f_int n, f_int kd, f_double *ab, f_int ldab, f_double *anorm, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dpbcon_(&uplo, &n, &kd, ab, &ldab, anorm, &rcond, work, iwork, &info, 1);
}
void pbcon(char uplo, f_int n, f_int kd, f_cfloat *ab, f_int ldab, f_float *anorm, f_float rcond, f_cfloat *work, f_float *rwork, ref f_int info) {
    cpbcon_(&uplo, &n, &kd, ab, &ldab, anorm, &rcond, work, rwork, &info, 1);
}
void pbcon(char uplo, f_int n, f_int kd, f_cdouble *ab, f_int ldab, f_double *anorm, f_double rcond, f_cdouble *work, f_double *rwork, ref f_int info) {
    zpbcon_(&uplo, &n, &kd, ab, &ldab, anorm, &rcond, work, rwork, &info, 1);
}

/// Computes row and column scalings to equilibrate a symmetric
/// positive definite band matrix and reduce its condition number.
void pbequ(char uplo, f_int n, f_int kd, f_float *ab, f_int ldab, f_float *s, f_float *scond, f_float *amax, ref f_int info) {
    spbequ_(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info, 1);
}
void pbequ(char uplo, f_int n, f_int kd, f_double *ab, f_int ldab, f_double *s, f_double *scond, f_double *amax, ref f_int info) {
    dpbequ_(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info, 1);
}
void pbequ(char uplo, f_int n, f_int kd, f_cfloat *ab, f_int ldab, f_float *s, f_float *scond, f_float *amax, ref f_int info) {
    cpbequ_(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info, 1);
}
void pbequ(char uplo, f_int n, f_int kd, f_cdouble *ab, f_int ldab, f_double *s, f_double *scond, f_double *amax, ref f_int info) {
    zpbequ_(&uplo, &n, &kd, ab, &ldab, s, scond, amax, &info, 1);
}

/// Improves the computed solution to a symmetric positive
/// definite banded system of linear equations AX=B, and provides
/// forward and backward error bounds for the solution.
void pbrfs(char uplo, f_int n, f_int kd, f_int nrhs, f_float *ab, f_int ldab, f_float *afb, f_int ldafb, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    spbrfs_(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void pbrfs(char uplo, f_int n, f_int kd, f_int nrhs, f_double *ab, f_int ldab, f_double *afb, f_int ldafb, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dpbrfs_(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void pbrfs(char uplo, f_int n, f_int kd, f_int nrhs, f_cfloat *ab, f_int ldab, f_cfloat *afb, f_int ldafb, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cpbrfs_(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void pbrfs(char uplo, f_int n, f_int kd, f_int nrhs, f_cdouble *ab, f_int ldab, f_cdouble *afb, f_int ldafb, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zpbrfs_(&uplo, &n, &kd, &nrhs, ab, &ldab, afb, &ldafb, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Computes a split Cholesky factorization of a real symmetric positive
/// definite band matrix.
void pbstf(char uplo, f_int n, f_int kd, f_float *ab, f_int ldab, ref f_int info) {
    spbstf_(&uplo, &n, &kd, ab, &ldab, &info, 1);
}
void pbstf(char uplo, f_int n, f_int kd, f_double *ab, f_int ldab, ref f_int info) {
    dpbstf_(&uplo, &n, &kd, ab, &ldab, &info, 1);
}
void pbstf(char uplo, f_int n, f_int kd, f_cfloat *ab, f_int ldab, ref f_int info) {
    cpbstf_(&uplo, &n, &kd, ab, &ldab, &info, 1);
}
void pbstf(char uplo, f_int n, f_int kd, f_cdouble *ab, f_int ldab, ref f_int info) {
    zpbstf_(&uplo, &n, &kd, ab, &ldab, &info, 1);
}

/// Computes the Cholesky factorization of a symmetric
/// positive definite band matrix.
void pbtrf(char uplo, f_int n, f_int kd, f_float *ab, f_int ldab, ref f_int info) {
    spbtrf_(&uplo, &n, &kd, ab, &ldab, &info, 1);
}
void pbtrf(char uplo, f_int n, f_int kd, f_double *ab, f_int ldab, ref f_int info) {
    dpbtrf_(&uplo, &n, &kd, ab, &ldab, &info, 1);
}
void pbtrf(char uplo, f_int n, f_int kd, f_cfloat *ab, f_int ldab, ref f_int info) {
    cpbtrf_(&uplo, &n, &kd, ab, &ldab, &info, 1);
}
void pbtrf(char uplo, f_int n, f_int kd, f_cdouble *ab, f_int ldab, ref f_int info) {
    zpbtrf_(&uplo, &n, &kd, ab, &ldab, &info, 1);
}

/// Solves a symmetric positive definite banded system
/// of linear equations AX=B, using the Cholesky factorization
/// computed by SPBTRF.
void pbtrs(char uplo, f_int n, f_int kd, f_int nrhs, f_float *ab, f_int ldab, f_float *b, f_int ldb, ref f_int info) {
    spbtrs_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1);
}
void pbtrs(char uplo, f_int n, f_int kd, f_int nrhs, f_double *ab, f_int ldab, f_double *b, f_int ldb, ref f_int info) {
    dpbtrs_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1);
}
void pbtrs(char uplo, f_int n, f_int kd, f_int nrhs, f_cfloat *ab, f_int ldab, f_cfloat *b, f_int ldb, ref f_int info) {
    cpbtrs_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1);
}
void pbtrs(char uplo, f_int n, f_int kd, f_int nrhs, f_cdouble *ab, f_int ldab, f_cdouble *b, f_int ldb, ref f_int info) {
    zpbtrs_(&uplo, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1);
}

/// Estimates the reciprocal of the condition number of a
/// symmetric positive definite matrix, using the
/// Cholesky factorization computed by SPOTRF.
void pocon(char uplo, f_int n, f_float *a, f_int lda, f_float *anorm, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    spocon_(&uplo, &n, a, &lda, anorm, &rcond, work, iwork, &info, 1);
}
void pocon(char uplo, f_int n, f_double *a, f_int lda, f_double *anorm, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dpocon_(&uplo, &n, a, &lda, anorm, &rcond, work, iwork, &info, 1);
}
void pocon(char uplo, f_int n, f_cfloat *a, f_int lda, f_float *anorm, f_float rcond, f_cfloat *work, f_float *rwork, ref f_int info) {
    cpocon_(&uplo, &n, a, &lda, anorm, &rcond, work, rwork, &info, 1);
}
void pocon(char uplo, f_int n, f_cdouble *a, f_int lda, f_double *anorm, f_double rcond, f_cdouble *work, f_double *rwork, ref f_int info) {
    zpocon_(&uplo, &n, a, &lda, anorm, &rcond, work, rwork, &info, 1);
}

/// Computes row and column scalings to equilibrate a symmetric
/// positive definite matrix and reduce its condition number.
void poequ(f_int n, f_float *a, f_int lda, f_float *s, f_float *scond, f_float *amax, ref f_int info) {
    spoequ_(&n, a, &lda, s, scond, amax, &info);
}
void poequ(f_int n, f_double *a, f_int lda, f_double *s, f_double *scond, f_double *amax, ref f_int info) {
    dpoequ_(&n, a, &lda, s, scond, amax, &info);
}
void poequ(f_int n, f_cfloat *a, f_int lda, f_float *s, f_float *scond, f_float *amax, ref f_int info) {
    cpoequ_(&n, a, &lda, s, scond, amax, &info);
}
void poequ(f_int n, f_cdouble *a, f_int lda, f_double *s, f_double *scond, f_double *amax, ref f_int info) {
    zpoequ_(&n, a, &lda, s, scond, amax, &info);
}

/// Improves the computed solution to a symmetric positive
/// definite system of linear equations AX=B, and provides forward
/// and backward error bounds for the solution.
void porfs(char uplo, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *af, f_int ldaf, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    sporfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void porfs(char uplo, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *af, f_int ldaf, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dporfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void porfs(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *af, f_int ldaf, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cporfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void porfs(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *af, f_int ldaf, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zporfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Computes the Cholesky factorization of a symmetric
/// positive definite matrix.
void potrf(char uplo, f_int n, f_float *a, f_int lda, ref f_int info) {
    spotrf_(&uplo, &n, a, &lda, &info, 1);
}
void potrf(char uplo, f_int n, f_double *a, f_int lda, ref f_int info) {
    dpotrf_(&uplo, &n, a, &lda, &info, 1);
}
void potrf(char uplo, f_int n, f_cfloat *a, f_int lda, ref f_int info) {
    cpotrf_(&uplo, &n, a, &lda, &info, 1);
}
void potrf(char uplo, f_int n, f_cdouble *a, f_int lda, ref f_int info) {
    zpotrf_(&uplo, &n, a, &lda, &info, 1);
}

/// Computes the inverse of a symmetric positive definite
/// matrix, using the Cholesky factorization computed by SPOTRF.
void potri(char uplo, f_int n, f_float *a, f_int lda, ref f_int info) {
    spotri_(&uplo, &n, a, &lda, &info, 1);
}
void potri(char uplo, f_int n, f_double *a, f_int lda, ref f_int info) {
    dpotri_(&uplo, &n, a, &lda, &info, 1);
}
void potri(char uplo, f_int n, f_cfloat *a, f_int lda, ref f_int info) {
    cpotri_(&uplo, &n, a, &lda, &info, 1);
}
void potri(char uplo, f_int n, f_cdouble *a, f_int lda, ref f_int info) {
    zpotri_(&uplo, &n, a, &lda, &info, 1);
}

/// Solves a symmetric positive definite system of linear
/// equations AX=B, using the Cholesky factorization computed by
/// SPOTRF.
void potrs(char uplo, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, ref f_int info) {
    spotrs_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
}
void potrs(char uplo, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, ref f_int info) {
    dpotrs_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
}
void potrs(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, ref f_int info) {
    cpotrs_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
}
void potrs(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, ref f_int info) {
    zpotrs_(&uplo, &n, &nrhs, a, &lda, b, &ldb, &info, 1);
}

/// Estimates the reciprocal of the condition number of a
/// symmetric positive definite matrix in packed storage,
/// using the Cholesky factorization computed by SPPTRF.
void ppcon(char uplo, f_int n, f_float *ap, f_float *anorm, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    sppcon_(&uplo, &n, ap, anorm, &rcond, work, iwork, &info, 1);
}
void ppcon(char uplo, f_int n, f_double *ap, f_double *anorm, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dppcon_(&uplo, &n, ap, anorm, &rcond, work, iwork, &info, 1);
}
void ppcon(char uplo, f_int n, f_cfloat *ap, f_float *anorm, f_float rcond, f_cfloat *work, f_float *rwork, ref f_int info) {
    cppcon_(&uplo, &n, ap, anorm, &rcond, work, rwork, &info, 1);
}
void ppcon(char uplo, f_int n, f_cdouble *ap, f_double *anorm, f_double rcond, f_cdouble *work, f_double *rwork, ref f_int info) {
    zppcon_(&uplo, &n, ap, anorm, &rcond, work, rwork, &info, 1);
}

/// Computes row and column scalings to equilibrate a symmetric
/// positive definite matrix in packed storage and reduce its condition
/// number.
void ppequ(char uplo, f_int n, f_float *ap, f_float *s, f_float *scond, f_float *amax, ref f_int info) {
    sppequ_(&uplo, &n, ap, s, scond, amax, &info, 1);
}
void ppequ(char uplo, f_int n, f_double *ap, f_double *s, f_double *scond, f_double *amax, ref f_int info) {
    dppequ_(&uplo, &n, ap, s, scond, amax, &info, 1);
}
void ppequ(char uplo, f_int n, f_cfloat *ap, f_float *s, f_float *scond, f_float *amax, ref f_int info) {
    cppequ_(&uplo, &n, ap, s, scond, amax, &info, 1);
}
void ppequ(char uplo, f_int n, f_cdouble *ap, f_double *s, f_double *scond, f_double *amax, ref f_int info) {
    zppequ_(&uplo, &n, ap, s, scond, amax, &info, 1);
}

/// Improves the computed solution to a symmetric positive
/// definite system of linear equations AX=B, where A is held in
/// packed storage, and provides forward and backward error bounds
/// for the solution.
void pprfs(char uplo, f_int n, f_int nrhs, f_float *ap, f_float *afp, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    spprfs_(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void pprfs(char uplo, f_int n, f_int nrhs, f_double *ap, f_double *afp, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dpprfs_(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void pprfs(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *afp, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cpprfs_(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void pprfs(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *afp, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zpprfs_(&uplo, &n, &nrhs, ap, afp, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Computes the Cholesky factorization of a symmetric
/// positive definite matrix in packed storage.
void pptrf(char uplo, f_int n, f_float *ap, ref f_int info) {
    spptrf_(&uplo, &n, ap, &info, 1);
}
void pptrf(char uplo, f_int n, f_double *ap, ref f_int info) {
    dpptrf_(&uplo, &n, ap, &info, 1);
}
void pptrf(char uplo, f_int n, f_cfloat *ap, ref f_int info) {
    cpptrf_(&uplo, &n, ap, &info, 1);
}
void pptrf(char uplo, f_int n, f_cdouble *ap, ref f_int info) {
    zpptrf_(&uplo, &n, ap, &info, 1);
}

/// Computes the inverse of a symmetric positive definite
/// matrix in packed storage, using the Cholesky factorization computed
/// by SPPTRF.
void pptri(char uplo, f_int n, f_float *ap, ref f_int info) {
    spptri_(&uplo, &n, ap, &info, 1);
}
void pptri(char uplo, f_int n, f_double *ap, ref f_int info) {
    dpptri_(&uplo, &n, ap, &info, 1);
}
void pptri(char uplo, f_int n, f_cfloat *ap, ref f_int info) {
    cpptri_(&uplo, &n, ap, &info, 1);
}
void pptri(char uplo, f_int n, f_cdouble *ap, ref f_int info) {
    zpptri_(&uplo, &n, ap, &info, 1);
}

/// Solves a symmetric positive definite system of linear
/// equations AX=B, where A is held in packed storage, using the
/// Cholesky factorization computed by SPPTRF.
void pptrs(char uplo, f_int n, f_int nrhs, f_float *ap, f_float *b, f_int ldb, ref f_int info) {
    spptrs_(&uplo, &n, &nrhs, ap, b, &ldb, &info, 1);
}
void pptrs(char uplo, f_int n, f_int nrhs, f_double *ap, f_double *b, f_int ldb, ref f_int info) {
    dpptrs_(&uplo, &n, &nrhs, ap, b, &ldb, &info, 1);
}
void pptrs(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *b, f_int ldb, ref f_int info) {
    cpptrs_(&uplo, &n, &nrhs, ap, b, &ldb, &info, 1);
}
void pptrs(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *b, f_int ldb, ref f_int info) {
    zpptrs_(&uplo, &n, &nrhs, ap, b, &ldb, &info, 1);
}

/// Computes the reciprocal of the condition number of a
/// symmetric positive definite tridiagonal matrix,
/// using the LDL**H factorization computed by SPTTRF.
void ptcon(f_int n, f_float *d, f_float *e, f_float *anorm, f_float rcond, f_float *work, ref f_int info) {
    sptcon_(&n, d, e, anorm, &rcond, work, &info);
}
void ptcon(f_int n, f_double *d, f_double *e, f_double *anorm, f_double rcond, f_double *work, ref f_int info) {
    dptcon_(&n, d, e, anorm, &rcond, work, &info);
}
void ptcon(f_int n, f_float *d, f_cfloat *e, f_float *anorm, f_float rcond, f_float *rwork, ref f_int info) {
    cptcon_(&n, d, e, anorm, &rcond, rwork, &info);
}
void ptcon(f_int n, f_double *d, f_cdouble *e, f_double *anorm, f_double rcond, f_double *rwork, ref f_int info) {
    zptcon_(&n, d, e, anorm, &rcond, rwork, &info);
}

/// Computes all eigenvalues and eigenvectors of a real symmetric
/// positive definite tridiagonal matrix, by computing the SVD of
/// its bidiagonal Cholesky factor.
void pteqr(char compz, f_int n, f_float *d, f_float *e, f_float *z, f_int ldz, f_float *work, ref f_int info) {
    spteqr_(&compz, &n, d, e, z, &ldz, work, &info, 1);
}
void pteqr(char compz, f_int n, f_double *d, f_double *e, f_double *z, f_int ldz, f_double *work, ref f_int info) {
    dpteqr_(&compz, &n, d, e, z, &ldz, work, &info, 1);
}
void pteqr(char compz, f_int n, f_float *d, f_float *e, f_cfloat *z, f_int ldz, f_float *work, ref f_int info) {
    cpteqr_(&compz, &n, d, e, z, &ldz, work, &info, 1);
}
void pteqr(char compz, f_int n, f_double *d, f_double *e, f_cdouble *z, f_int ldz, f_double *work, ref f_int info) {
    zpteqr_(&compz, &n, d, e, z, &ldz, work, &info, 1);
}

/// Improves the computed solution to a symmetric positive
/// definite tridiagonal system of linear equations AX=B, and provides
/// forward and backward error bounds for the solution.
void ptrfs(f_int n, f_int nrhs, f_float *d, f_float *e, f_float *df, f_float *ef, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, ref f_int info) {
    sptrfs_(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, &info);
}
void ptrfs(f_int n, f_int nrhs, f_double *d, f_double *e, f_double *df, f_double *ef, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, ref f_int info) {
    dptrfs_(&n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, &info);
}
void ptrfs(char uplo, f_int n, f_int nrhs, f_float *d, f_cfloat *e, f_float *df, f_cfloat *ef, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cptrfs_(&uplo, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void ptrfs(char uplo, f_int n, f_int nrhs, f_double *d, f_cdouble *e, f_double *df, f_cdouble *ef, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zptrfs_(&uplo, &n, &nrhs, d, e, df, ef, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Computes the LDL**H factorization of a symmetric
/// positive definite tridiagonal matrix.
void pttrf(f_int n, f_float *d, f_float *e, ref f_int info) {
    spttrf_(&n, d, e, &info);
}
void pttrf(f_int n, f_double *d, f_double *e, ref f_int info) {
    dpttrf_(&n, d, e, &info);
}
void pttrf(f_int n, f_float *d, f_cfloat *e, ref f_int info) {
    cpttrf_(&n, d, e, &info);
}
void pttrf(f_int n, f_double *d, f_cdouble *e, ref f_int info) {
    zpttrf_(&n, d, e, &info);
}

/// Solves a symmetric positive definite tridiagonal
/// system of linear equations, using the LDL**H factorization
/// computed by SPTTRF.
void pttrs(f_int n, f_int nrhs, f_float *d, f_float *e, f_float *b, f_int ldb, ref f_int info) {
    spttrs_(&n, &nrhs, d, e, b, &ldb, &info);
}
void pttrs(f_int n, f_int nrhs, f_double *d, f_double *e, f_double *b, f_int ldb, ref f_int info) {
    dpttrs_(&n, &nrhs, d, e, b, &ldb, &info);
}
void pttrs(char uplo, f_int n, f_int nrhs, f_float *d, f_cfloat *e, f_cfloat *b, f_int ldb, ref f_int info) {
    cpttrs_(&uplo, &n, &nrhs, d, e, b, &ldb, &info, 1);
}
void pttrs(char uplo, f_int n, f_int nrhs, f_double *d, f_cdouble *e, f_cdouble *b, f_int ldb, ref f_int info) {
    zpttrs_(&uplo, &n, &nrhs, d, e, b, &ldb, &info, 1);
}

/// Reduces a real symmetric-definite banded generalized eigenproblem
/// A x = lambda B x to standard form, where B has been factorized by
/// SPBSTF (Crawford's algorithm).
void sbgst(char vect, char uplo, f_int n, f_int ka, f_int kb, f_float *ab, f_int ldab, f_float *bb, f_int ldbb, f_float *x, f_int ldx, f_float *work, ref f_int info) {
    ssbgst_(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, &info, 1, 1);
}
void sbgst(char vect, char uplo, f_int n, f_int ka, f_int kb, f_double *ab, f_int ldab, f_double *bb, f_int ldbb, f_double *x, f_int ldx, f_double *work, ref f_int info) {
    dsbgst_(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, &info, 1, 1);
}

/// Reduces a complex Hermitian-definite banded generalized eigenproblem
/// A x = lambda B x to standard form, where B has been factorized by
/// CPBSTF (Crawford's algorithm).
void hbgst(char vect, char uplo, f_int n, f_int ka, f_int kb, f_cfloat *ab, f_int ldab, f_cfloat *bb, f_int ldbb, f_cfloat *x, f_int ldx, f_cfloat *work, f_float *rwork, ref f_int info) {
    chbgst_(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, rwork, &info, 1, 1);
}
void hbgst(char vect, char uplo, f_int n, f_int ka, f_int kb, f_cdouble *ab, f_int ldab, f_cdouble *bb, f_int ldbb, f_cdouble *x, f_int ldx, f_cdouble *work, f_double *rwork, ref f_int info) {
    zhbgst_(&vect, &uplo, &n, &ka, &kb, ab, &ldab, bb, &ldbb, x, &ldx, work, rwork, &info, 1, 1);
}

/// Reduces a symmetric band matrix to real symmetric
/// tridiagonal form by an orthogonal similarity transformation.
void sbtrd(char vect, char uplo, f_int n, f_int kd, f_float *ab, f_int ldab, f_float *d, f_float *e, f_float *q, f_int ldq, f_float *work, ref f_int info) {
    ssbtrd_(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info, 1, 1);
}
void sbtrd(char vect, char uplo, f_int n, f_int kd, f_double *ab, f_int ldab, f_double *d, f_double *e, f_double *q, f_int ldq, f_double *work, ref f_int info) {
    dsbtrd_(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info, 1, 1);
}

/// Reduces a Hermitian band matrix to real symmetric
/// tridiagonal form by a unitary similarity transformation.
void hbtrd(char vect, char uplo, f_int n, f_int kd, f_cfloat *ab, f_int ldab, f_float *d, f_float *e, f_cfloat *q, f_int ldq, f_cfloat *work, ref f_int info) {
    chbtrd_(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info, 1, 1);
}
void hbtrd(char vect, char uplo, f_int n, f_int kd, f_cdouble *ab, f_int ldab, f_double *d, f_double *e, f_cdouble *q, f_int ldq, f_cdouble *work, ref f_int info) {
    zhbtrd_(&vect, &uplo, &n, &kd, ab, &ldab, d, e, q, &ldq, work, &info, 1, 1);
}

/// Estimates the reciprocal of the condition number of a
/// real symmetric indefinite
/// matrix in packed storage, using the factorization computed
/// by SSPTRF.
void spcon(char uplo, f_int n, f_float *ap, f_int *ipiv, f_float *anorm, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    sspcon_(&uplo, &n, ap, ipiv, anorm, &rcond, work, iwork, &info, 1);
}
void spcon(char uplo, f_int n, f_double *ap, f_int *ipiv, f_double *anorm, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dspcon_(&uplo, &n, ap, ipiv, anorm, &rcond, work, iwork, &info, 1);
}
void spcon(char uplo, f_int n, f_cfloat *ap, f_int *ipiv, f_float *anorm, f_float rcond, f_cfloat *work, ref f_int info) {
    cspcon_(&uplo, &n, ap, ipiv, anorm, &rcond, work, &info, 1);
}
void spcon(char uplo, f_int n, f_cdouble *ap, f_int *ipiv, f_double *anorm, f_double rcond, f_cdouble *work, ref f_int info) {
    zspcon_(&uplo, &n, ap, ipiv, anorm, &rcond, work, &info, 1);
}

/// Estimates the reciprocal of the condition number of a
/// complex Hermitian indefinite
/// matrix in packed storage, using the factorization computed
/// by CHPTRF.
void hpcon(char uplo, f_int n, f_cfloat *ap, f_int *ipiv, f_float *anorm, f_float rcond, f_cfloat *work, ref f_int info) {
    chpcon_(&uplo, &n, ap, ipiv, anorm, &rcond, work, &info, 1);
}
void hpcon(char uplo, f_int n, f_cdouble *ap, f_int *ipiv, f_double *anorm, f_double rcond, f_cdouble *work, ref f_int info) {
    zhpcon_(&uplo, &n, ap, ipiv, anorm, &rcond, work, &info, 1);
}

/// Reduces a symmetric-definite generalized eigenproblem
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard
/// form,  where A and B are held in packed storage, and B has been
/// factorized by SPPTRF.
void spgst(f_int itype, char uplo, f_int n, f_float *ap, f_float *bp, ref f_int info) {
    sspgst_(&itype, &uplo, &n, ap, bp, &info, 1);
}
void spgst(f_int itype, char uplo, f_int n, f_double *ap, f_double *bp, ref f_int info) {
    dspgst_(&itype, &uplo, &n, ap, bp, &info, 1);
}

/// Reduces a Hermitian-definite generalized eigenproblem
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard
/// form,  where A and B are held in packed storage, and B has been
/// factorized by CPPTRF.
void hpgst(f_int itype, char uplo, f_int n, f_cfloat *ap, f_cfloat *bp, ref f_int info) {
    chpgst_(&itype, &uplo, &n, ap, bp, &info, 1);
}
void hpgst(f_int itype, char uplo, f_int n, f_cdouble *ap, f_cdouble *bp, ref f_int info) {
    zhpgst_(&itype, &uplo, &n, ap, bp, &info, 1);
}

/// Improves the computed solution to a real
/// symmetric indefinite system of linear equations
/// AX=B, where A is held in packed storage, and provides forward
/// and backward error bounds for the solution.
void sprfs(char uplo, f_int n, f_int nrhs, f_float *ap, f_float *afp, f_int *ipiv, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    ssprfs_(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void sprfs(char uplo, f_int n, f_int nrhs, f_double *ap, f_double *afp, f_int *ipiv, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dsprfs_(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void sprfs(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *afp, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    csprfs_(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void sprfs(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *afp, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zsprfs_(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Improves the computed solution to a complex
/// Hermitian indefinite system of linear equations
/// AX=B, where A is held in packed storage, and provides forward
/// and backward error bounds for the solution.
void hprfs(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *afp, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    chprfs_(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void hprfs(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *afp, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zhprfs_(&uplo, &n, &nrhs, ap, afp, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Reduces a symmetric matrix in packed storage to real
/// symmetric tridiagonal form by an orthogonal similarity
/// transformation.
void sptrd(char uplo, f_int n, f_float *ap, f_float *d, f_float *e, f_float *tau, ref f_int info) {
    ssptrd_(&uplo, &n, ap, d, e, tau, &info, 1);
}
void sptrd(char uplo, f_int n, f_double *ap, f_double *d, f_double *e, f_double *tau, ref f_int info) {
    dsptrd_(&uplo, &n, ap, d, e, tau, &info, 1);
}

/// Reduces a Hermitian matrix in packed storage to real
/// symmetric tridiagonal form by a unitary similarity
/// transformation.
void hptrd(char uplo, f_int n, f_cfloat *ap, f_float *d, f_float *e, f_cfloat *tau, ref f_int info) {
    chptrd_(&uplo, &n, ap, d, e, tau, &info, 1);
}
void hptrd(char uplo, f_int n, f_cdouble *ap, f_double *d, f_double *e, f_cdouble *tau, ref f_int info) {
    zhptrd_(&uplo, &n, ap, d, e, tau, &info, 1);
}

/// Computes the factorization of a real
/// symmetric-indefinite matrix in packed storage,
/// using the diagonal pivoting method.
void sptrf(char uplo, f_int n, f_float *ap, f_int *ipiv, ref f_int info) {
    ssptrf_(&uplo, &n, ap, ipiv, &info, 1);
}
void sptrf(char uplo, f_int n, f_double *ap, f_int *ipiv, ref f_int info) {
    dsptrf_(&uplo, &n, ap, ipiv, &info, 1);
}
void sptrf(char uplo, f_int n, f_cfloat *ap, f_int *ipiv, ref f_int info) {
    csptrf_(&uplo, &n, ap, ipiv, &info, 1);
}
void sptrf(char uplo, f_int n, f_cdouble *ap, f_int *ipiv, ref f_int info) {
    zsptrf_(&uplo, &n, ap, ipiv, &info, 1);
}

/// Computes the factorization of a complex
/// Hermitian-indefinite matrix in packed storage,
/// using the diagonal pivoting method.
void hptrf(char uplo, f_int n, f_cfloat *ap, f_int *ipiv, ref f_int info) {
    chptrf_(&uplo, &n, ap, ipiv, &info, 1);
}
void hptrf(char uplo, f_int n, f_cdouble *ap, f_int *ipiv, ref f_int info) {
    zhptrf_(&uplo, &n, ap, ipiv, &info, 1);
}

/// Computes the inverse of a real symmetric
/// indefinite matrix in packed storage, using the factorization
/// computed by SSPTRF.
void sptri(char uplo, f_int n, f_float *ap, f_int *ipiv, f_float *work, ref f_int info) {
    ssptri_(&uplo, &n, ap, ipiv, work, &info, 1);
}
void sptri(char uplo, f_int n, f_double *ap, f_int *ipiv, f_double *work, ref f_int info) {
    dsptri_(&uplo, &n, ap, ipiv, work, &info, 1);
}
void sptri(char uplo, f_int n, f_cfloat *ap, f_int *ipiv, f_cfloat *work, ref f_int info) {
    csptri_(&uplo, &n, ap, ipiv, work, &info, 1);
}
void sptri(char uplo, f_int n, f_cdouble *ap, f_int *ipiv, f_cdouble *work, ref f_int info) {
    zsptri_(&uplo, &n, ap, ipiv, work, &info, 1);
}

/// Computes the inverse of a complex
/// Hermitian indefinite matrix in packed storage, using the factorization
/// computed by CHPTRF.
void hptri(char uplo, f_int n, f_cfloat *ap, f_int *ipiv, f_cfloat *work, ref f_int info) {
    chptri_(&uplo, &n, ap, ipiv, work, &info, 1);
}
void hptri(char uplo, f_int n, f_cdouble *ap, f_int *ipiv, f_cdouble *work, ref f_int info) {
    zhptri_(&uplo, &n, ap, ipiv, work, &info, 1);
}

/// Solves a real symmetric
/// indefinite system of linear equations AX=B, where A is held
/// in packed storage, using the factorization computed
/// by SSPTRF.
void sptrs(char uplo, f_int n, f_int nrhs, f_float *ap, f_int *ipiv, f_float *b, f_int ldb, ref f_int info) {
    ssptrs_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}
void sptrs(char uplo, f_int n, f_int nrhs, f_double *ap, f_int *ipiv, f_double *b, f_int ldb, ref f_int info) {
    dsptrs_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}
void sptrs(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    csptrs_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}
void sptrs(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zsptrs_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}

/// Solves a complex Hermitian
/// indefinite system of linear equations AX=B, where A is held
/// in packed storage, using the factorization computed
/// by CHPTRF.
void hptrs(char uplo, f_int n, f_int nrhs, f_cfloat *ap, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    chptrs_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}
void hptrs(char uplo, f_int n, f_int nrhs, f_cdouble *ap, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zhptrs_(&uplo, &n, &nrhs, ap, ipiv, b, &ldb, &info, 1);
}

/// Computes selected eigenvalues of a real symmetric tridiagonal
/// matrix by bisection.
void stebz(char range, char order, f_int n, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_float *d, f_float *e, f_int m, f_int nsplit, f_float *w, f_int iblock, f_int isplit, f_float *work, f_int *iwork, ref f_int info) {
    sstebz_(&range, &order, &n, vl, vu, &il, &iu, abstol, d, e, &m, &nsplit, w, &iblock, &isplit, work, iwork, &info, 1, 1);
}
void stebz(char range, char order, f_int n, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_double *d, f_double *e, f_int m, f_int nsplit, f_double *w, f_int iblock, f_int isplit, f_double *work, f_int *iwork, ref f_int info) {
    dstebz_(&range, &order, &n, vl, vu, &il, &iu, abstol, d, e, &m, &nsplit, w, &iblock, &isplit, work, iwork, &info, 1, 1);
}

/// Computes all eigenvalues and, optionally, eigenvectors of a
/// symmetric tridiagonal matrix using the divide and conquer algorithm.
void stedc(char compz, f_int n, f_float *d, f_float *e, f_float *z, f_int ldz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    sstedc_(&compz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info, 1);
}
void stedc(char compz, f_int n, f_double *d, f_double *e, f_double *z, f_int ldz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dstedc_(&compz, &n, d, e, z, &ldz, work, &lwork, iwork, &liwork, &info, 1);
}
void stedc(char compz, f_int n, f_float *d, f_float *e, f_cfloat *z, f_int ldz, f_cfloat *work, f_int lwork, f_float *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    cstedc_(&compz, &n, d, e, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1);
}
void stedc(char compz, f_int n, f_double *d, f_double *e, f_cdouble *z, f_int ldz, f_cdouble *work, f_int lwork, f_double *rwork, f_int lrwork, f_int *iwork, f_int liwork, ref f_int info) {
    zstedc_(&compz, &n, d, e, z, &ldz, work, &lwork, rwork, &lrwork, iwork, &liwork, &info, 1);
}

/// Computes selected eigenvalues and, optionally, eigenvectors of a
/// symmetric tridiagonal matrix.  The eigenvalues are computed by the
/// dqds algorithm, while eigenvectors are computed from various "good"
/// LDL^T representations (also known as Relatively Robust Representations.)
void stegr(char jobz, char range, f_int n, f_float *d, f_float *e, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_float *z, f_int ldz, f_int isuppz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    sstegr_(&jobz, &range, &n, d, e, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, &isuppz, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void stegr(char jobz, char range, f_int n, f_double *d, f_double *e, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_double *z, f_int ldz, f_int isuppz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dstegr_(&jobz, &range, &n, d, e, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, &isuppz, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void stegr(char jobz, char range, f_int n, f_float *d, f_float *e, f_float *vl, f_float *vu, f_int il, f_int iu, f_float *abstol, f_int m, f_float *w, f_cfloat *z, f_int ldz, f_int isuppz, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    cstegr_(&jobz, &range, &n, d, e, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, &isuppz, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void stegr(char jobz, char range, f_int n, f_double *d, f_double *e, f_double *vl, f_double *vu, f_int il, f_int iu, f_double *abstol, f_int m, f_double *w, f_cdouble *z, f_int ldz, f_int isuppz, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    zstegr_(&jobz, &range, &n, d, e, vl, vu, &il, &iu, abstol, &m, w, z, &ldz, &isuppz, work, &lwork, iwork, &liwork, &info, 1, 1);
}

/// Computes selected eigenvectors of a real symmetric tridiagonal
/// matrix by inverse iteration.
void stein(f_int n, f_float *d, f_float *e, f_int m, f_float *w, f_int iblock, f_int isplit, f_float *z, f_int ldz, f_float *work, f_int *iwork, f_int ifail, ref f_int info) {
    sstein_(&n, d, e, &m, w, &iblock, &isplit, z, &ldz, work, iwork, &ifail, &info);
}
void stein(f_int n, f_double *d, f_double *e, f_int m, f_double *w, f_int iblock, f_int isplit, f_double *z, f_int ldz, f_double *work, f_int *iwork, f_int ifail, ref f_int info) {
    dstein_(&n, d, e, &m, w, &iblock, &isplit, z, &ldz, work, iwork, &ifail, &info);
}
void stein(f_int n, f_float *d, f_float *e, f_int m, f_float *w, f_int iblock, f_int isplit, f_cfloat *z, f_int ldz, f_float *work, f_int *iwork, f_int ifail, ref f_int info) {
    cstein_(&n, d, e, &m, w, &iblock, &isplit, z, &ldz, work, iwork, &ifail, &info);
}
void stein(f_int n, f_double *d, f_double *e, f_int m, f_double *w, f_int iblock, f_int isplit, f_cdouble *z, f_int ldz, f_double *work, f_int *iwork, f_int ifail, ref f_int info) {
    zstein_(&n, d, e, &m, w, &iblock, &isplit, z, &ldz, work, iwork, &ifail, &info);
}

/// Computes all eigenvalues and eigenvectors of a real symmetric
/// tridiagonal matrix, using the implicit QL or QR algorithm.
void steqr(char compz, f_int n, f_float *d, f_float *e, f_float *z, f_int ldz, f_float *work, ref f_int info) {
    ssteqr_(&compz, &n, d, e, z, &ldz, work, &info, 1);
}
void steqr(char compz, f_int n, f_double *d, f_double *e, f_double *z, f_int ldz, f_double *work, ref f_int info) {
    dsteqr_(&compz, &n, d, e, z, &ldz, work, &info, 1);
}
void steqr(char compz, f_int n, f_float *d, f_float *e, f_cfloat *z, f_int ldz, f_float *work, ref f_int info) {
    csteqr_(&compz, &n, d, e, z, &ldz, work, &info, 1);
}
void steqr(char compz, f_int n, f_double *d, f_double *e, f_cdouble *z, f_int ldz, f_double *work, ref f_int info) {
    zsteqr_(&compz, &n, d, e, z, &ldz, work, &info, 1);
}

/// Computes all eigenvalues of a real symmetric tridiagonal matrix,
/// using a root-free variant of the QL or QR algorithm.
void sterf(f_int n, f_float *d, f_float *e, ref f_int info) {
    ssterf_(&n, d, e, &info);
}
void sterf(f_int n, f_double *d, f_double *e, ref f_int info) {
    dsterf_(&n, d, e, &info);
}

/// Estimates the reciprocal of the condition number of a
/// real symmetric indefinite matrix,
/// using the factorization computed by SSYTRF.
void sycon(char uplo, f_int n, f_float *a, f_int lda, f_int *ipiv, f_float *anorm, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    ssycon_(&uplo, &n, a, &lda, ipiv, anorm, &rcond, work, iwork, &info, 1);
}
void sycon(char uplo, f_int n, f_double *a, f_int lda, f_int *ipiv, f_double *anorm, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dsycon_(&uplo, &n, a, &lda, ipiv, anorm, &rcond, work, iwork, &info, 1);
}
void sycon(char uplo, f_int n, f_cfloat *a, f_int lda, f_int *ipiv, f_float *anorm, f_float rcond, f_cfloat *work, ref f_int info) {
    csycon_(&uplo, &n, a, &lda, ipiv, anorm, &rcond, work, &info, 1);
}
void sycon(char uplo, f_int n, f_cdouble *a, f_int lda, f_int *ipiv, f_double *anorm, f_double rcond, f_cdouble *work, ref f_int info) {
    zsycon_(&uplo, &n, a, &lda, ipiv, anorm, &rcond, work, &info, 1);
}

/// Estimates the reciprocal of the condition number of a
/// complex Hermitian indefinite matrix,
/// using the factorization computed by CHETRF.
void hecon(char uplo, f_int n, f_cfloat *a, f_int lda, f_int *ipiv, f_float *anorm, f_float rcond, f_cfloat *work, ref f_int info) {
    checon_(&uplo, &n, a, &lda, ipiv, anorm, &rcond, work, &info, 1);
}
void hecon(char uplo, f_int n, f_cdouble *a, f_int lda, f_int *ipiv, f_double *anorm, f_double rcond, f_cdouble *work, ref f_int info) {
    zhecon_(&uplo, &n, a, &lda, ipiv, anorm, &rcond, work, &info, 1);
}

/// Reduces a symmetric-definite generalized eigenproblem
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard
/// form, where B has been factorized by SPOTRF.
void sygst(f_int itype, char uplo, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, ref f_int info) {
    ssygst_(&itype, &uplo, &n, a, &lda, b, &ldb, &info, 1);
}
void sygst(f_int itype, char uplo, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, ref f_int info) {
    dsygst_(&itype, &uplo, &n, a, &lda, b, &ldb, &info, 1);
}

/// Reduces a Hermitian-definite generalized eigenproblem
/// Ax= lambda Bx,  ABx= lambda x,  or BAx= lambda x, to standard
/// form, where B has been factorized by CPOTRF.
void hegst(f_int itype, char uplo, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, ref f_int info) {
    chegst_(&itype, &uplo, &n, a, &lda, b, &ldb, &info, 1);
}
void hegst(f_int itype, char uplo, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, ref f_int info) {
    zhegst_(&itype, &uplo, &n, a, &lda, b, &ldb, &info, 1);
}

/// Improves the computed solution to a real
/// symmetric indefinite system of linear equations
/// AX=B, and provides forward and backward error bounds for the
/// solution.
void syrfs(char uplo, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *af, f_int ldaf, f_int *ipiv, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    ssyrfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void syrfs(char uplo, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *af, f_int ldaf, f_int *ipiv, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dsyrfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1);
}
void syrfs(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *af, f_int ldaf, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    csyrfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void syrfs(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *af, f_int ldaf, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zsyrfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Improves the computed solution to a complex
/// Hermitian indefinite system of linear equations
/// AX=B, and provides forward and backward error bounds for the
/// solution.
void herfs(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *af, f_int ldaf, f_int *ipiv, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    cherfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}
void herfs(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *af, f_int ldaf, f_int *ipiv, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    zherfs_(&uplo, &n, &nrhs, a, &lda, af, &ldaf, ipiv, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1);
}

/// Reduces a symmetric matrix to real symmetric tridiagonal
/// form by an orthogonal similarity transformation.
void sytrd(char uplo, f_int n, f_float *a, f_int lda, f_float *d, f_float *e, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    ssytrd_(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info, 1);
}
void sytrd(char uplo, f_int n, f_double *a, f_int lda, f_double *d, f_double *e, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dsytrd_(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info, 1);
}

/// Reduces a Hermitian matrix to real symmetric tridiagonal
/// form by an orthogonal/unitary similarity transformation.
void hetrd(char uplo, f_int n, f_cfloat *a, f_int lda, f_float *d, f_float *e, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    chetrd_(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info, 1);
}
void hetrd(char uplo, f_int n, f_cdouble *a, f_int lda, f_double *d, f_double *e, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    zhetrd_(&uplo, &n, a, &lda, d, e, tau, work, &lwork, &info, 1);
}

/// Computes the factorization of a real symmetric-indefinite matrix,
/// using the diagonal pivoting method.
void sytrf(char uplo, f_int n, f_float *a, f_int lda, f_int *ipiv, f_float *work, f_int lwork, ref f_int info) {
    ssytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info, 1);
}
void sytrf(char uplo, f_int n, f_double *a, f_int lda, f_int *ipiv, f_double *work, f_int lwork, ref f_int info) {
    dsytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info, 1);
}
void sytrf(char uplo, f_int n, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *work, f_int lwork, ref f_int info) {
    csytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info, 1);
}
void sytrf(char uplo, f_int n, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *work, f_int lwork, ref f_int info) {
    zsytrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info, 1);
}

/// Computes the factorization of a complex Hermitian-indefinite matrix,
/// using the diagonal pivoting method.
void hetrf(char uplo, f_int n, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *work, f_int lwork, ref f_int info) {
    chetrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info, 1);
}
void hetrf(char uplo, f_int n, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *work, f_int lwork, ref f_int info) {
    zhetrf_(&uplo, &n, a, &lda, ipiv, work, &lwork, &info, 1);
}

/// Computes the inverse of a real symmetric indefinite matrix,
/// using the factorization computed by SSYTRF.
void sytri(char uplo, f_int n, f_float *a, f_int lda, f_int *ipiv, f_float *work, ref f_int info) {
    ssytri_(&uplo, &n, a, &lda, ipiv, work, &info, 1);
}
void sytri(char uplo, f_int n, f_double *a, f_int lda, f_int *ipiv, f_double *work, ref f_int info) {
    dsytri_(&uplo, &n, a, &lda, ipiv, work, &info, 1);
}
void sytri(char uplo, f_int n, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *work, ref f_int info) {
    csytri_(&uplo, &n, a, &lda, ipiv, work, &info, 1);
}
void sytri(char uplo, f_int n, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *work, ref f_int info) {
    zsytri_(&uplo, &n, a, &lda, ipiv, work, &info, 1);
}

/// Computes the inverse of a complex Hermitian indefinite matrix,
/// using the factorization computed by CHETRF.
void hetri(char uplo, f_int n, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *work, ref f_int info) {
    chetri_(&uplo, &n, a, &lda, ipiv, work, &info, 1);
}
void hetri(char uplo, f_int n, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *work, ref f_int info) {
    zhetri_(&uplo, &n, a, &lda, ipiv, work, &info, 1);
}

/// Solves a real symmetric indefinite system of linear equations AX=B,
/// using the factorization computed by SSPTRF.
void sytrs(char uplo, f_int n, f_int nrhs, f_float *a, f_int lda, f_int *ipiv, f_float *b, f_int ldb, ref f_int info) {
    ssytrs_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}
void sytrs(char uplo, f_int n, f_int nrhs, f_double *a, f_int lda, f_int *ipiv, f_double *b, f_int ldb, ref f_int info) {
    dsytrs_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}
void sytrs(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    csytrs_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}
void sytrs(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zsytrs_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}

/// Solves a complex Hermitian indefinite system of linear equations AX=B,
/// using the factorization computed by CHPTRF.
void hetrs(char uplo, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_int *ipiv, f_cfloat *b, f_int ldb, ref f_int info) {
    chetrs_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}
void hetrs(char uplo, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_int *ipiv, f_cdouble *b, f_int ldb, ref f_int info) {
    zhetrs_(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, &info, 1);
}

/// Estimates the reciprocal of the condition number of a triangular
/// band matrix, in either the 1-norm or the infinity-norm.
void tbcon(char norm, char uplo, char diag, f_int n, f_int kd, f_float *ab, f_int ldab, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    stbcon_(&norm, &uplo, &diag, &n, &kd, ab, &ldab, &rcond, work, iwork, &info, 1, 1, 1);
}
void tbcon(char norm, char uplo, char diag, f_int n, f_int kd, f_double *ab, f_int ldab, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dtbcon_(&norm, &uplo, &diag, &n, &kd, ab, &ldab, &rcond, work, iwork, &info, 1, 1, 1);
}
void tbcon(char norm, char uplo, char diag, f_int n, f_int kd, f_cfloat *ab, f_int ldab, f_float rcond, f_cfloat *work, f_float *rwork, ref f_int info) {
    ctbcon_(&norm, &uplo, &diag, &n, &kd, ab, &ldab, &rcond, work, rwork, &info, 1, 1, 1);
}
void tbcon(char norm, char uplo, char diag, f_int n, f_int kd, f_cdouble *ab, f_int ldab, f_double rcond, f_cdouble *work, f_double *rwork, ref f_int info) {
    ztbcon_(&norm, &uplo, &diag, &n, &kd, ab, &ldab, &rcond, work, rwork, &info, 1, 1, 1);
}

/// Provides forward and backward error bounds for the solution
/// of a triangular banded system of linear equations AX=B,
/// A**T X=B or A**H X=B.
void tbrfs(char uplo, char trans, char diag, f_int n, f_int kd, f_int nrhs, f_float *ab, f_int ldab, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    stbrfs_(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void tbrfs(char uplo, char trans, char diag, f_int n, f_int kd, f_int nrhs, f_double *ab, f_int ldab, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dtbrfs_(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void tbrfs(char uplo, char trans, char diag, f_int n, f_int kd, f_int nrhs, f_cfloat *ab, f_int ldab, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    ctbrfs_(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1, 1, 1);
}
void tbrfs(char uplo, char trans, char diag, f_int n, f_int kd, f_int nrhs, f_cdouble *ab, f_int ldab, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    ztbrfs_(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1, 1, 1);
}

/// Solves a triangular banded system of linear equations AX=B,
/// A**T X=B or A**H X=B.
void tbtrs(char uplo, char trans, char diag, f_int n, f_int kd, f_int nrhs, f_float *ab, f_int ldab, f_float *b, f_int ldb, ref f_int info) {
    stbtrs_(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1, 1, 1);
}
void tbtrs(char uplo, char trans, char diag, f_int n, f_int kd, f_int nrhs, f_double *ab, f_int ldab, f_double *b, f_int ldb, ref f_int info) {
    dtbtrs_(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1, 1, 1);
}
void tbtrs(char uplo, char trans, char diag, f_int n, f_int kd, f_int nrhs, f_cfloat *ab, f_int ldab, f_cfloat *b, f_int ldb, ref f_int info) {
    ctbtrs_(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1, 1, 1);
}
void tbtrs(char uplo, char trans, char diag, f_int n, f_int kd, f_int nrhs, f_cdouble *ab, f_int ldab, f_cdouble *b, f_int ldb, ref f_int info) {
    ztbtrs_(&uplo, &trans, &diag, &n, &kd, &nrhs, ab, &ldab, b, &ldb, &info, 1, 1, 1);
}

/// Computes some or all of the right and/or left generalized eigenvectors
/// of a pair of upper triangular matrices.
void tgevc(char side, char howmny, f_int select, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_int mm, f_int m, f_float *work, ref f_int info) {
    stgevc_(&side, &howmny, &select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info, 1, 1);
}
void tgevc(char side, char howmny, f_int select, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_int mm, f_int m, f_double *work, ref f_int info) {
    dtgevc_(&side, &howmny, &select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info, 1, 1);
}
void tgevc(char side, char howmny, f_int select, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_int mm, f_int m, f_cfloat *work, f_float *rwork, ref f_int info) {
    ctgevc_(&side, &howmny, &select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &info, 1, 1);
}
void tgevc(char side, char howmny, f_int select, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_int mm, f_int m, f_cdouble *work, f_double *rwork, ref f_int info) {
    ztgevc_(&side, &howmny, &select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &info, 1, 1);
}

/// Reorders the generalized real Schur decomposition of a real
/// matrix pair (A,B) using an orthogonal equivalence transformation
/// so that the diagonal block of (A,B) with row index IFST is moved
/// to row ILST.
void tgexc(f_int wantq, f_int wantz, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *q, f_int ldq, f_float *z, f_int ldz, f_int ifst, f_int ilst, f_float *work, f_int lwork, ref f_int info) {
    stgexc_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &ifst, &ilst, work, &lwork, &info);
}
void tgexc(f_int wantq, f_int wantz, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *q, f_int ldq, f_double *z, f_int ldz, f_int ifst, f_int ilst, f_double *work, f_int lwork, ref f_int info) {
    dtgexc_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &ifst, &ilst, work, &lwork, &info);
}
void tgexc(f_int wantq, f_int wantz, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *q, f_int ldq, f_cfloat *z, f_int ldz, f_int ifst, f_int ilst, ref f_int info) {
    ctgexc_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &ifst, &ilst, &info);
}
void tgexc(f_int wantq, f_int wantz, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *q, f_int ldq, f_cdouble *z, f_int ldz, f_int ifst, f_int ilst, ref f_int info) {
    ztgexc_(&wantq, &wantz, &n, a, &lda, b, &ldb, q, &ldq, z, &ldz, &ifst, &ilst, &info);
}

/// Reorders the generalized real Schur decomposition of a real
/// matrix pair (A, B) so that a selected cluster of eigenvalues
/// appears in the leading diagonal blocks of the upper quasi-triangular
/// matrix A and the upper triangular B.
void tgsen(f_int ijob, f_int wantq, f_int wantz, f_int select, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *alphar, f_float *alphai, f_float *betav, f_float *q, f_int ldq, f_float *z, f_int ldz, f_int m, f_float *pl, f_float *pr, f_float *dif, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    stgsen_(&ijob, &wantq, &wantz, &select, &n, a, &lda, b, &ldb, alphar, alphai, betav, q, &ldq, z, &ldz, &m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
}
void tgsen(f_int ijob, f_int wantq, f_int wantz, f_int select, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *alphar, f_double *alphai, f_double *betav, f_double *q, f_int ldq, f_double *z, f_int ldz, f_int m, f_double *pl, f_double *pr, f_double *dif, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dtgsen_(&ijob, &wantq, &wantz, &select, &n, a, &lda, b, &ldb, alphar, alphai, betav, q, &ldq, z, &ldz, &m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
}
void tgsen(f_int ijob, f_int wantq, f_int wantz, f_int select, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *alphav, f_cfloat *betav, f_cfloat *q, f_int ldq, f_cfloat *z, f_int ldz, f_int m, f_float *pl, f_float *pr, f_float *dif, f_cfloat *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    ctgsen_(&ijob, &wantq, &wantz, &select, &n, a, &lda, b, &ldb, alphav, betav, q, &ldq, z, &ldz, &m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
}
void tgsen(f_int ijob, f_int wantq, f_int wantz, f_int select, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *alphav, f_cdouble *betav, f_cdouble *q, f_int ldq, f_cdouble *z, f_int ldz, f_int m, f_double *pl, f_double *pr, f_double *dif, f_cdouble *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    ztgsen_(&ijob, &wantq, &wantz, &select, &n, a, &lda, b, &ldb, alphav, betav, q, &ldq, z, &ldz, &m, pl, pr, dif, work, &lwork, iwork, &liwork, &info);
}

/// Computes the generalized singular value decomposition of two real
/// upper triangular (or trapezoidal) matrices as output by SGGSVP.
void tgsja(char jobu, char jobv, char jobq, f_int m, f_int p, f_int n, f_int k, f_int l, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *tola, f_float *tolb, f_float *alphav, f_float *betav, f_float *u, f_int ldu, f_float *v, f_int ldv, f_float *q, f_int ldq, f_float *work, f_int ncycle, ref f_int info) {
    stgsja_(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, tola, tolb, alphav, betav, u, &ldu, v, &ldv, q, &ldq, work, &ncycle, &info, 1, 1, 1);
}
void tgsja(char jobu, char jobv, char jobq, f_int m, f_int p, f_int n, f_int k, f_int l, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *tola, f_double *tolb, f_double *alphav, f_double *betav, f_double *u, f_int ldu, f_double *v, f_int ldv, f_double *q, f_int ldq, f_double *work, f_int ncycle, ref f_int info) {
    dtgsja_(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, tola, tolb, alphav, betav, u, &ldu, v, &ldv, q, &ldq, work, &ncycle, &info, 1, 1, 1);
}
void tgsja(char jobu, char jobv, char jobq, f_int m, f_int p, f_int n, f_int k, f_int l, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_float *tola, f_float *tolb, f_float *alphav, f_float *betav, f_cfloat *u, f_int ldu, f_cfloat *v, f_int ldv, f_cfloat *q, f_int ldq, f_cfloat *work, f_int ncycle, ref f_int info) {
    ctgsja_(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, tola, tolb, alphav, betav, u, &ldu, v, &ldv, q, &ldq, work, &ncycle, &info, 1, 1, 1);
}
void tgsja(char jobu, char jobv, char jobq, f_int m, f_int p, f_int n, f_int k, f_int l, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_double *tola, f_double *tolb, f_double *alphav, f_double *betav, f_cdouble *u, f_int ldu, f_cdouble *v, f_int ldv, f_cdouble *q, f_int ldq, f_cdouble *work, f_int ncycle, ref f_int info) {
    ztgsja_(&jobu, &jobv, &jobq, &m, &p, &n, &k, &l, a, &lda, b, &ldb, tola, tolb, alphav, betav, u, &ldu, v, &ldv, q, &ldq, work, &ncycle, &info, 1, 1, 1);
}

/// Estimates reciprocal condition numbers for specified
/// eigenvalues and/or eigenvectors of a matrix pair (A, B) in
/// generalized real Schur canonical form, as returned by SGGES.
void tgsna(char job, char howmny, f_int select, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_float *s, f_float *dif, f_int mm, f_int m, f_float *work, f_int lwork, f_int *iwork, ref f_int info) {
    stgsna_(&job, &howmny, &select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, &m, work, &lwork, iwork, &info, 1, 1);
}
void tgsna(char job, char howmny, f_int select, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_double *s, f_double *dif, f_int mm, f_int m, f_double *work, f_int lwork, f_int *iwork, ref f_int info) {
    dtgsna_(&job, &howmny, &select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, &m, work, &lwork, iwork, &info, 1, 1);
}
void tgsna(char job, char howmny, f_int select, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_float *s, f_float *dif, f_int mm, f_int m, f_cfloat *work, f_int lwork, f_int *iwork, ref f_int info) {
    ctgsna_(&job, &howmny, &select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, &m, work, &lwork, iwork, &info, 1, 1);
}
void tgsna(char job, char howmny, f_int select, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_double *s, f_double *dif, f_int mm, f_int m, f_cdouble *work, f_int lwork, f_int *iwork, ref f_int info) {
    ztgsna_(&job, &howmny, &select, &n, a, &lda, b, &ldb, vl, &ldvl, vr, &ldvr, s, dif, &mm, &m, work, &lwork, iwork, &info, 1, 1);
}

/// Solves the generalized Sylvester equation.
void tgsyl(char trans, f_int ijob, f_int m, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *c, f_int ldc, f_float *d, f_int ldd, f_float *e, f_int lde, f_float *f, f_int ldf, f_float *scale, f_float *dif, f_float *work, f_int lwork, f_int *iwork, ref f_int info) {
    stgsyl_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work, &lwork, iwork, &info, 1);
}
void tgsyl(char trans, f_int ijob, f_int m, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *c, f_int ldc, f_double *d, f_int ldd, f_double *e, f_int lde, f_double *f, f_int ldf, f_double *scale, f_double *dif, f_double *work, f_int lwork, f_int *iwork, ref f_int info) {
    dtgsyl_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work, &lwork, iwork, &info, 1);
}
void tgsyl(char trans, f_int ijob, f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *c, f_int ldc, f_cfloat *d, f_int ldd, f_cfloat *e, f_int lde, f_cfloat *f, f_int ldf, f_float *scale, f_float *dif, f_cfloat *work, f_int lwork, f_int *iwork, ref f_int info) {
    ctgsyl_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work, &lwork, iwork, &info, 1);
}
void tgsyl(char trans, f_int ijob, f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *c, f_int ldc, f_cdouble *d, f_int ldd, f_cdouble *e, f_int lde, f_cdouble *f, f_int ldf, f_double *scale, f_double *dif, f_cdouble *work, f_int lwork, f_int *iwork, ref f_int info) {
    ztgsyl_(&trans, &ijob, &m, &n, a, &lda, b, &ldb, c, &ldc, d, &ldd, e, &lde, f, &ldf, scale, dif, work, &lwork, iwork, &info, 1);
}

/// Estimates the reciprocal of the condition number of a triangular
/// matrix in packed storage, in either the 1-norm or the infinity-norm.
void tpcon(char norm, char uplo, char diag, f_int n, f_float *ap, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    stpcon_(&norm, &uplo, &diag, &n, ap, &rcond, work, iwork, &info, 1, 1, 1);
}
void tpcon(char norm, char uplo, char diag, f_int n, f_double *ap, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dtpcon_(&norm, &uplo, &diag, &n, ap, &rcond, work, iwork, &info, 1, 1, 1);
}
void tpcon(char norm, char uplo, char diag, f_int n, f_cfloat *ap, f_float rcond, f_cfloat *work, f_float *rwork, ref f_int info) {
    ctpcon_(&norm, &uplo, &diag, &n, ap, &rcond, work, rwork, &info, 1, 1, 1);
}
void tpcon(char norm, char uplo, char diag, f_int n, f_cdouble *ap, f_double rcond, f_cdouble *work, f_double *rwork, ref f_int info) {
    ztpcon_(&norm, &uplo, &diag, &n, ap, &rcond, work, rwork, &info, 1, 1, 1);
}

/// Provides forward and backward error bounds for the solution
/// of a triangular system of linear equations AX=B, A**T X=B or
/// A**H X=B, where A is held in packed storage.
void tprfs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_float *ap, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    stprfs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void tprfs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_double *ap, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dtprfs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void tprfs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    ctprfs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1, 1, 1);
}
void tprfs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    ztprfs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1, 1, 1);
}

///  Computes the inverse of a triangular matrix in packed storage.
void tptri(char uplo, char diag, f_int n, f_float *ap, ref f_int info) {
    stptri_(&uplo, &diag, &n, ap, &info, 1, 1);
}
void tptri(char uplo, char diag, f_int n, f_double *ap, ref f_int info) {
    dtptri_(&uplo, &diag, &n, ap, &info, 1, 1);
}
void tptri(char uplo, char diag, f_int n, f_cfloat *ap, ref f_int info) {
    ctptri_(&uplo, &diag, &n, ap, &info, 1, 1);
}
void tptri(char uplo, char diag, f_int n, f_cdouble *ap, ref f_int info) {
    ztptri_(&uplo, &diag, &n, ap, &info, 1, 1);
}

/// Solves a triangular system of linear equations AX=B,
/// A**T X=B or A**H X=B, where A is held in packed storage.
void tptrs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_float *ap, f_float *b, f_int ldb, ref f_int info) {
    stptrs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info, 1, 1, 1);
}
void tptrs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_double *ap, f_double *b, f_int ldb, ref f_int info) {
    dtptrs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info, 1, 1, 1);
}
void tptrs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_cfloat *ap, f_cfloat *b, f_int ldb, ref f_int info) {
    ctptrs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info, 1, 1, 1);
}
void tptrs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_cdouble *ap, f_cdouble *b, f_int ldb, ref f_int info) {
    ztptrs_(&uplo, &trans, &diag, &n, &nrhs, ap, b, &ldb, &info, 1, 1, 1);
}

/// Estimates the reciprocal of the condition number of a triangular
/// matrix, in either the 1-norm or the infinity-norm.
void trcon(char norm, char uplo, char diag, f_int n, f_float *a, f_int lda, f_float rcond, f_float *work, f_int *iwork, ref f_int info) {
    strcon_(&norm, &uplo, &diag, &n, a, &lda, &rcond, work, iwork, &info, 1, 1, 1);
}
void trcon(char norm, char uplo, char diag, f_int n, f_double *a, f_int lda, f_double rcond, f_double *work, f_int *iwork, ref f_int info) {
    dtrcon_(&norm, &uplo, &diag, &n, a, &lda, &rcond, work, iwork, &info, 1, 1, 1);
}
void trcon(char norm, char uplo, char diag, f_int n, f_cfloat *a, f_int lda, f_float rcond, f_cfloat *work, f_float *rwork, ref f_int info) {
    ctrcon_(&norm, &uplo, &diag, &n, a, &lda, &rcond, work, rwork, &info, 1, 1, 1);
}
void trcon(char norm, char uplo, char diag, f_int n, f_cdouble *a, f_int lda, f_double rcond, f_cdouble *work, f_double *rwork, ref f_int info) {
    ztrcon_(&norm, &uplo, &diag, &n, a, &lda, &rcond, work, rwork, &info, 1, 1, 1);
}

/// Computes some or all of the right and/or left eigenvectors of
/// an upper quasi-triangular matrix.
void trevc(char side, char howmny, f_int select, f_int n, f_float *t, f_int ldt, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_int mm, f_int m, f_float *work, ref f_int info) {
    strevc_(&side, &howmny, &select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info, 1, 1);
}
void trevc(char side, char howmny, f_int select, f_int n, f_double *t, f_int ldt, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_int mm, f_int m, f_double *work, ref f_int info) {
    dtrevc_(&side, &howmny, &select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, &info, 1, 1);
}
void trevc(char side, char howmny, f_int select, f_int n, f_cfloat *t, f_int ldt, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_int mm, f_int m, f_cfloat *work, f_float *rwork, ref f_int info) {
    ctrevc_(&side, &howmny, &select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &info, 1, 1);
}
void trevc(char side, char howmny, f_int select, f_int n, f_cdouble *t, f_int ldt, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_int mm, f_int m, f_cdouble *work, f_double *rwork, ref f_int info) {
    ztrevc_(&side, &howmny, &select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, &mm, &m, work, rwork, &info, 1, 1);
}

/// Reorders the Schur factorization of a matrix by an orthogonal
/// similarity transformation.
void trexc(char compq, f_int n, f_float *t, f_int ldt, f_float *q, f_int ldq, f_int ifst, f_int ilst, f_float *work, ref f_int info) {
    strexc_(&compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, work, &info, 1);
}
void trexc(char compq, f_int n, f_double *t, f_int ldt, f_double *q, f_int ldq, f_int ifst, f_int ilst, f_double *work, ref f_int info) {
    dtrexc_(&compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, work, &info, 1);
}
void trexc(char compq, f_int n, f_cfloat *t, f_int ldt, f_cfloat *q, f_int ldq, f_int ifst, f_int ilst, ref f_int info) {
    ctrexc_(&compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, &info, 1);
}
void trexc(char compq, f_int n, f_cdouble *t, f_int ldt, f_cdouble *q, f_int ldq, f_int ifst, f_int ilst, ref f_int info) {
    ztrexc_(&compq, &n, t, &ldt, q, &ldq, &ifst, &ilst, &info, 1);
}

/// Provides forward and backward error bounds for the solution
/// of a triangular system of linear equations A X=B, A**T X=B or
/// A**H X=B.
void trrfs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *x, f_int ldx, f_float *ferr, f_float *berr, f_float *work, f_int *iwork, ref f_int info) {
    strrfs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void trrfs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *x, f_int ldx, f_double *ferr, f_double *berr, f_double *work, f_int *iwork, ref f_int info) {
    dtrrfs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, iwork, &info, 1, 1, 1);
}
void trrfs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *x, f_int ldx, f_float *ferr, f_float *berr, f_cfloat *work, f_float *rwork, ref f_int info) {
    ctrrfs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1, 1, 1);
}
void trrfs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *x, f_int ldx, f_double *ferr, f_double *berr, f_cdouble *work, f_double *rwork, ref f_int info) {
    ztrrfs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, x, &ldx, ferr, berr, work, rwork, &info, 1, 1, 1);
}

/// Reorders the Schur factorization of a matrix in order to find
/// an orthonormal basis of a right invariant subspace corresponding
/// to selected eigenvalues, and returns reciprocal condition numbers
/// (sensitivities) of the average of the cluster of eigenvalues
/// and of the invariant subspace.
void trsen(char job, char compq, f_int select, f_int n, f_float *t, f_int ldt, f_float *q, f_int ldq, f_float *wr, f_float *wi, f_int m, f_float *s, f_float *sep, f_float *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    strsen_(&job, &compq, &select, &n, t, &ldt, q, &ldq, wr, wi, &m, s, sep, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void trsen(char job, char compq, f_int select, f_int n, f_double *t, f_int ldt, f_double *q, f_int ldq, f_double *wr, f_double *wi, f_int m, f_double *s, f_double *sep, f_double *work, f_int lwork, f_int *iwork, f_int liwork, ref f_int info) {
    dtrsen_(&job, &compq, &select, &n, t, &ldt, q, &ldq, wr, wi, &m, s, sep, work, &lwork, iwork, &liwork, &info, 1, 1);
}
void trsen(char job, char compq, f_int select, f_int n, f_cfloat *t, f_int ldt, f_cfloat *q, f_int ldq, f_cfloat *w, f_int m, f_float *s, f_float *sep, f_cfloat *work, f_int lwork, ref f_int info) {
    ctrsen_(&job, &compq, &select, &n, t, &ldt, q, &ldq, w, &m, s, sep, work, &lwork, &info, 1, 1);
}
void trsen(char job, char compq, f_int select, f_int n, f_cdouble *t, f_int ldt, f_cdouble *q, f_int ldq, f_cdouble *w, f_int m, f_double *s, f_double *sep, f_cdouble *work, f_int lwork, ref f_int info) {
    ztrsen_(&job, &compq, &select, &n, t, &ldt, q, &ldq, w, &m, s, sep, work, &lwork, &info, 1, 1);
}

/// Estimates the reciprocal condition numbers (sensitivities)
/// of selected eigenvalues and eigenvectors of an upper
/// quasi-triangular matrix.
void trsna(char job, char howmny, f_int select, f_int n, f_float *t, f_int ldt, f_float *vl, f_int ldvl, f_float *vr, f_int ldvr, f_float *s, f_float *sep, f_int mm, f_int m, f_float *work, f_int ldwork, f_int *iwork, ref f_int info) {
    strsna_(&job, &howmny, &select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, &m, work, &ldwork, iwork, &info, 1, 1);
}
void trsna(char job, char howmny, f_int select, f_int n, f_double *t, f_int ldt, f_double *vl, f_int ldvl, f_double *vr, f_int ldvr, f_double *s, f_double *sep, f_int mm, f_int m, f_double *work, f_int ldwork, f_int *iwork, ref f_int info) {
    dtrsna_(&job, &howmny, &select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, &m, work, &ldwork, iwork, &info, 1, 1);
}
void trsna(char job, char howmny, f_int select, f_int n, f_cfloat *t, f_int ldt, f_cfloat *vl, f_int ldvl, f_cfloat *vr, f_int ldvr, f_float *s, f_float *sep, f_int mm, f_int m, f_cfloat *work, f_int ldwork, f_float *rwork, ref f_int info) {
    ctrsna_(&job, &howmny, &select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, &m, work, &ldwork, rwork, &info, 1, 1);
}
void trsna(char job, char howmny, f_int select, f_int n, f_cdouble *t, f_int ldt, f_cdouble *vl, f_int ldvl, f_cdouble *vr, f_int ldvr, f_double *s, f_double *sep, f_int mm, f_int m, f_cdouble *work, f_int ldwork, f_double *rwork, ref f_int info) {
    ztrsna_(&job, &howmny, &select, &n, t, &ldt, vl, &ldvl, vr, &ldvr, s, sep, &mm, &m, work, &ldwork, rwork, &info, 1, 1);
}

/// Solves the Sylvester matrix equation A X +/- X B=C where A
/// and B are upper quasi-triangular, and may be transposed.
void trsyl(char trana, char tranb, f_int isgn, f_int m, f_int n, f_float *a, f_int lda, f_float *b, f_int ldb, f_float *c, f_int ldc, f_float *scale, ref f_int info) {
    strsyl_(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info, 1, 1);
}
void trsyl(char trana, char tranb, f_int isgn, f_int m, f_int n, f_double *a, f_int lda, f_double *b, f_int ldb, f_double *c, f_int ldc, f_double *scale, ref f_int info) {
    dtrsyl_(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info, 1, 1);
}
void trsyl(char trana, char tranb, f_int isgn, f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, f_cfloat *c, f_int ldc, f_float *scale, ref f_int info) {
    ctrsyl_(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info, 1, 1);
}
void trsyl(char trana, char tranb, f_int isgn, f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, f_cdouble *c, f_int ldc, f_double *scale, ref f_int info) {
    ztrsyl_(&trana, &tranb, &isgn, &m, &n, a, &lda, b, &ldb, c, &ldc, scale, &info, 1, 1);
}

/// Computes the inverse of a triangular matrix.
void trtri(char uplo, char diag, f_int n, f_float *a, f_int lda, ref f_int info) {
    strtri_(&uplo, &diag, &n, a, &lda, &info, 1, 1);
}
void trtri(char uplo, char diag, f_int n, f_double *a, f_int lda, ref f_int info) {
    dtrtri_(&uplo, &diag, &n, a, &lda, &info, 1, 1);
}
void trtri(char uplo, char diag, f_int n, f_cfloat *a, f_int lda, ref f_int info) {
    ctrtri_(&uplo, &diag, &n, a, &lda, &info, 1, 1);
}
void trtri(char uplo, char diag, f_int n, f_cdouble *a, f_int lda, ref f_int info) {
    ztrtri_(&uplo, &diag, &n, a, &lda, &info, 1, 1);
}

/// Solves a triangular system of linear equations AX=B,
/// A**T X=B or A**H X=B.
void trtrs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_float *a, f_int lda, f_float *b, f_int ldb, ref f_int info) {
    strtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info, 1, 1, 1);
}
void trtrs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_double *a, f_int lda, f_double *b, f_int ldb, ref f_int info) {
    dtrtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info, 1, 1, 1);
}
void trtrs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_cfloat *a, f_int lda, f_cfloat *b, f_int ldb, ref f_int info) {
    ctrtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info, 1, 1, 1);
}
void trtrs(char uplo, char trans, char diag, f_int n, f_int nrhs, f_cdouble *a, f_int lda, f_cdouble *b, f_int ldb, ref f_int info) {
    ztrtrs_(&uplo, &trans, &diag, &n, &nrhs, a, &lda, b, &ldb, &info, 1, 1, 1);
}

/// Computes an RQ factorization of an upper trapezoidal matrix.
void tzrqf(f_int m, f_int n, f_float *a, f_int lda, f_float *tau, ref f_int info) {
    stzrqf_(&m, &n, a, &lda, tau, &info);
}
void tzrqf(f_int m, f_int n, f_double *a, f_int lda, f_double *tau, ref f_int info) {
    dtzrqf_(&m, &n, a, &lda, tau, &info);
}
void tzrqf(f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *tau, ref f_int info) {
    ctzrqf_(&m, &n, a, &lda, tau, &info);
}
void tzrqf(f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *tau, ref f_int info) {
    ztzrqf_(&m, &n, a, &lda, tau, &info);
}

/// Computes an RZ factorization of an upper trapezoidal matrix
/// (blocked version of STZRQF).
void tzrzf(f_int m, f_int n, f_float *a, f_int lda, f_float *tau, f_float *work, f_int lwork, ref f_int info) {
    stzrzf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void tzrzf(f_int m, f_int n, f_double *a, f_int lda, f_double *tau, f_double *work, f_int lwork, ref f_int info) {
    dtzrzf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void tzrzf(f_int m, f_int n, f_cfloat *a, f_int lda, f_cfloat *tau, f_cfloat *work, f_int lwork, ref f_int info) {
    ctzrzf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}
void tzrzf(f_int m, f_int n, f_cdouble *a, f_int lda, f_cdouble *tau, f_cdouble *work, f_int lwork, ref f_int info) {
    ztzrzf_(&m, &n, a, &lda, tau, work, &lwork, &info);
}


/// Multiplies a general matrix by the unitary
/// transformation matrix from a reduction to tridiagonal form
/// determined by CHPTRD.
void upmtr(char side, char uplo, char trans, f_int m, f_int n, f_cfloat *ap, f_cfloat *tau, f_cfloat *c, f_int ldc, f_cfloat *work, ref f_int info) {
    cupmtr_(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info, 1, 1, 1);
}
void upmtr(char side, char uplo, char trans, f_int m, f_int n, f_cdouble *ap, f_cdouble *tau, f_cdouble *c, f_int ldc, f_cdouble *work, ref f_int info) {
    zupmtr_(&side, &uplo, &trans, &m, &n, ap, tau, c, &ldc, work, &info, 1, 1, 1);
}


//------------------------------------
//     ----- MISC routines -----
//------------------------------------

f_int ilaenv(f_int ispec, char *name, char *opts, f_int n1, f_int n2, f_int n3, f_int n4, f_int len_name, f_int len_opts) {
    return ilaenv_(&ispec, name, opts, &n1, &n2, &n3, &n4, len_name, len_opts);
}
void ilaenvset(f_int ispec, char *name, char *opts, f_int n1, f_int n2, f_int n3, f_int n4, f_int nvalue, ref f_int info, f_int len_name, f_int len_opts) {
    // hmm this doesn't seem to exist in the lib in -g debug builds for some reason
    //ilaenvset_(&ispec, name, opts, &n1, &n2, &n3, &n4, &nvalue, &info, len_name, len_opts);
}

///
f_float slamch(char[]cmach) {
    return slamch_(cmach.ptr, toInt(cmach.length));
}
f_double dlamch(char[]cmach) {
    return dlamch_(cmach.ptr, toInt(cmach.length));
}

///
lapack_float_ret_t second() {
    return second_();
}
f_double secnd() {
    return dsecnd_();
}


