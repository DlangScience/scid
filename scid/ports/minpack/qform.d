/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0 
*/
module scid.ports.minpack.qform;


private import std.algorithm: min;




/** This subroutine proceeds from the computed QR factorization of
    an m by n matrix a to accumulate the m by m orthogonal matrix
    Q from its factored form.

    Params:
        m = a positive integer input variable set to the number
            of rows of a and the order of q.
        n = a positive integer input variable set to the number
            of columns of a.
        q = an array of length m^2. On input the full lower trapezoid in
            the first min(m,n) columns of Q contains the factored form.
            on output Q has been accumulated into a square matrix.
        ldq = a positive integer input variable not less than m
            which specifies the leading dimension of the array q.
        wa = a work array of length m.
*/
void qform(Real)(size_t m, size_t n, Real* q, size_t ldq, Real* wa)
{
    size_t i, j, k, l;

    enum : Real
    {
        one = 1.0,
        zero = 0.0
    }


    // Zero out upper triangle of Q in the first min(m,n) columns.
    size_t minmn = min(m, n);
    if (minmn >= 2)
    {
        for (j=1; j<minmn; j++)
        {
            for (i=0; i<j; i++)
                q[i+j*ldq] = zero;
        }
    }


    // Initialize remaining columns to those of the identity matrix.
    if (m > n)
    {
        for (j=n; j<m; j++)
        {
            for (i=0; i<m; i++)
                q[i+j*ldq] = zero;
            q[j+j*ldq] = one;
        }
    }


    // Accumulate Q from its factored form.
    Real sum, temp;
    for (l=0; l<minmn; l++)
    {
        k = minmn - l - 1;
        for (i=k; i<m; i++)
        {
            wa[i] = q[i+k*ldq];
            q[i+k*ldq] = zero;
        }
        q[k+k*ldq] = one;

        if (wa[k] != zero)
        {
            for (j=k; j<m; j++)
            {
                sum = zero;
                for (i=k; i<m; i++)
                    sum += q[i+j*ldq]*wa[i];
                temp = sum/wa[k];
                for (i=k; i<m; i++)
                    q[i+j*ldq] = q[i+j*ldq] - temp*wa[i];
            }
        }
    }
}


unittest { alias qform!(real) rqform; }
