// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/minpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0
*/
module scid.ports.minpack.qrfac;


private import std.algorithm: min, max;
private import std.math: sqrt;

private import scid.ports.minpack.enorm;




/** This subroutine uses Householder transformations with column
    pivoting (optional) to compute a QR factorization of the
    m by n matrix A. That is, qrfac determines an orthogonal
    matrix Q, a permutation matrix P, and an upper trapezoidal
    matrix R with diagonal elements of nonincreasing magnitude,
    such that AP = QR. The Householder transformation for
    column k, k = 1,2,...,min(m,n), is of the form
    ---
                          t
          i - (1/u(k)) u u
    ---
    where u has zeros in the first k-1 positions. The form of
    this transformation and the method of pivoting first
    appeared in the corresponding LINPACK subroutine.

    Params:

      m = a positive integer input variable set to the number
        of rows of a.

      n = a positive integer input variable set to the number
        of columns of a.

      a = an m by n array. on input a contains the matrix for
        which the qr factorization is to be computed. on output
        the strict upper trapezoidal part of a contains the strict
        upper trapezoidal part of r, and the lower trapezoidal
        part of a contains a factored form of q (the non-trivial
        elements of the u vectors described above).

      lda = a positive integer input variable not less than m
        which specifies the leading dimension of the array a.

      pivot = a logical input variable. if pivot is set true,
        then column pivoting is enforced. if pivot is set false,
        then no column pivoting is done.

      ipvt = an integer output array of length lipvt. ipvt
        defines the permutation matrix p such that a*p = q*r.
        column j of p is column ipvt(j) of the identity matrix.
        if pivot is false, ipvt is not referenced.

      lipvt = a positive integer input variable. if pivot is false,
        then lipvt may be as small as 1. if pivot is true, then
        lipvt must be at least n.

      rdiag = an output array of length n which contains the
        diagonal elements of r.

      acnorm = an output array of length n which contains the
        norms of the corresponding columns of the input matrix a.
        if this information is not needed, then acnorm can coincide
        with rdiag.

      wa = a work array of length n. if pivot is false, then wa
        can coincide with rdiag.
*/
void qrfac(Real)(size_t m, size_t n, Real* a, size_t lda,
    bool pivot, size_t* ipvt, size_t lipvt, Real* rdiag, Real* acnorm,
    Real* wa)
{
    size_t i, j, ij, ikmax, jp1, k, kmax, minmn;
    Real ajnorm, sum, temp;
    
    enum : Real
    {
        one = 1.0,
        p05 = 0.05,
        zero = 0.0,

        // epsmch is the machine precision
        epsmch = Real.epsilon
    }


    // Compute the initial column norms and initialize several arrays.
    for (j=0; j<n; j++)
    {
        acnorm[j] = enorm(m, a+j*m);
        rdiag[j] = acnorm[j];
        wa[j] = rdiag[j];
        if (pivot) ipvt[j] = j;
    }


    // Reduce A to R with Householder transformations.
    minmn = min(m, n);
    for (j=0; j<minmn; j++)
    {
        if (pivot)
        {
            // Bring the column of largest norm into the pivot position.
            kmax = j;
            for (k=j; k<n; k++)
                if (rdiag[k] > rdiag[kmax]) kmax = k;

            if (kmax != j)
            {
                for (i=0; i<m; i++)
                {
                    ij = i + j*lda;
                    ikmax = i + kmax*lda;
                    temp = a[ij];
                    a[ij] = a[ikmax];
                    a[ikmax] = temp;
                }

                rdiag[kmax] = rdiag[j];
                wa[kmax] = wa[j];
                k = ipvt[j];
                ipvt[j] = ipvt[kmax];
                ipvt[kmax] = k;
            }
        }


        // Compute the Householder transformation to reduce the
        // j-th column of A to a multiple of the j-th unit vector.
        ajnorm = enorm(m-j, a+j+j*lda);
        if (ajnorm != zero)
        {
            if (a[j+j*lda] < zero) ajnorm = -ajnorm;
            for (i=j; i<m; i++)
                a[i+j*lda] /= ajnorm;
            a[j+j*lda] += one;


            // Apply the transformation to the remaining columns
            // and update the norms.
            jp1 = j + 1;
            if (n >= jp1)
            {
                for (k=jp1; k<n; k++)
                {
                    sum = zero;
                    for (i=j; i<m; i++)
                        sum += a[i+j*lda]*a[i+k*lda];
                    temp = sum/a[j+j*lda];
                    for (i=j; i<m; i++)
                        a[i+k*lda] -= temp*a[i+j*lda];
                    
                    if (pivot  &&  rdiag[k] != zero)
                    {
                        temp = a[j+k*lda]/rdiag[k];
                        rdiag[k] *= sqrt(max(zero, one-temp*temp));
                        if (p05*((rdiag[k]/wa[k])^^2) <= epsmch)
                        {
                            rdiag[k] = enorm(m-j-1, a+jp1+k*lda);
                            wa[k] = rdiag[k];
                        }
                    }
                }
            }
        }
        rdiag[j] = -ajnorm;
    }
}


unittest { alias qrfac!(real) rqrfac; }
