// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/minpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0
*/
module scid.ports.minpack.dogleg;


import std.algorithm: max, min;
import std.math;

import scid.ports.minpack.enorm;



/** Given an m by n matrix A, an n by n nonsingular diagonal matrix D,
    an m-vector b, and a positive number delta, the problem is to
    determine the convex combination x of the Gauss-Newton and scaled
    gradient directions that minimizes (Ax-b) in the least squares
    sense, subject to the restriction that the Euclidean norm of Dx
    be at most delta.

    This function completes the solution of the problem if it is
    provided with the necessary information from the QR factorization
    of A. That is, if A=QR, where Q has orthogonal columns and R is
    an upper triangular matrix, then dogleg expects the full upper
    triangle of R and the first n components of (Q^T)b.

    Params:
        n = a positive integer variable set to the order of R.
        r = the upper triangular matrix R.
        lr = (not documented)
        diag = a vector of length n which must contain the diagonal
            elements of the matrix D.
        qtb = a vector of length n which must contain the first n
            elements of the vector (Q^T)b.
        delta = a positive variable which specifies an upper bound
            on the Euclidean norm of Dx.
        x = an output vector of length n which contains the desired
            convex combination of the Gauss-Newton direction and
            the scaled gradient direction.
        wa1 = work array of length n.
        wa2 = work array of length n.
*/
void dogleg(Real)(size_t n, Real* r, size_t lr, Real* diag, Real* qtb,
    Real delta, Real* x, Real* wa1, Real* wa2)
{
    size_t i, j, jj, jp1, k, l;
    Real alpha, bnorm, gnorm, qnorm, sgnorm, sum, temp;
    enum : Real
    {
        one = 1.0,
        zero = 0.0,

        // epsmch is the machine precision.
        epsmch = Real.epsilon
    }


    
    // First, calculate the Gauss-Newton direction.
    jj = (n*(n+1))/2;
    for (k=0; k<n; k++)
    {
        j = n - k - 1;
        jp1 = j + 1;
        jj -= k + 1;
        l = jj + 1;

        sum = zero;
        if (n > jp1)
        {
            for (i=jp1; i<n; i++)
            {
                sum += r[l]*x[i];
                l++;
            }
        }
        
        temp = r[jj];
        if (temp == zero)
        {
            l = j;
            for (i=0; i<=j; i++)
            {
                temp = max(temp, abs(r[l]));
                l += n - i - 1;
            }
            temp *= epsmch;
            if (temp == zero) temp = epsmch;
        }

        x[j] = (qtb[j] - sum)/temp;
    }


    // Test whether the Gauss-Newton direction is acceptable.
    for (j=0; j<n; j++)
    {
        wa1[j] = zero;
        wa2[j] = diag[j]*x[j];
    }
    qnorm = enorm(n, wa2);
    if (qnorm <= delta) return;


    // The Gauss-Newton direction is not acceptable.
    // Next, calculate the scaled gradient direction.
    l = 0;
    for (j=0; j<n; j++)
    {
        temp = qtb[j];
        for (i=j; i<n; i++)
        {
            wa1[i] = wa1[i] + r[l]*temp;
            l++;
        }

        wa1[j] = wa1[j]/diag[j];
    }


    // Calculate the norm of the scaled gradient and test for
    // the special case in which the scaled gradient is zero.
    gnorm = enorm(n, wa1);
    sgnorm = zero;
    alpha = delta/qnorm;
    if (gnorm != zero)
    {
        // Calculate the point along the scaled gradient
        // at which the quadratic is minimized.
        for (j=0; j<n; j++)
            wa1[j] = (wa1[j]/gnorm)/diag[j];

        l = 0;
        for (j=0; j<n; j++)
        {
            sum = zero;
            for (i=j; i<n; i++)
            {
                sum += r[l]*wa1[i];
                l++;
            }
            wa2[j] = sum;
        }

        temp = enorm(n, wa2);
        sgnorm = (gnorm/temp)/temp;


        // Test whether the scaled gradient direction is acceptable.
        alpha = zero;
        if (sgnorm < delta)
        {
            // The scaled gradient direction is not acceptable.
            // Finally, calculate the point along the dogleg
            // at which the quadratic is minimized.
            bnorm = enorm(n, qtb);
            temp = (bnorm/gnorm)*(bnorm/qnorm)*(sgnorm/delta);
            temp = temp - (delta/qnorm)*((sgnorm/delta)^^2)
                + sqrt((temp-delta/qnorm)^^2
                    + (one-((delta/qnorm)^^2))*(one-((sgnorm/delta)^^2)));
            alpha = ((delta/qnorm)*(one - ((sgnorm/delta)^^2)))/temp;
        }
    }


    // Form appropriate convex combination of the Gauss-Newton
    // direction and the scaled gradient direction.
    temp = (one - alpha)*min(sgnorm, delta);
    for (j=0; j<n; j++)
        x[j] = temp*wa1[j] + alpha*x[j];

    return;
}


unittest { alias dogleg!(real) realDog; }
