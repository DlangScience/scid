// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/minpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0
*/
module scid.ports.minpack.r1mpyq;


private import std.math: abs, sqrt;




/** Given an m by n matrix A, this subroutine computes AQ where
    Q is the product of 2*(n-1) transformations
    ---
    gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1)
    ---
    and gv(i), gw(i) are Givens rotations in the (i,n) plane which
    eliminate elements in the i-th and n-th planes, respectively.
    Q itself is not given, rather the information to recover the
    gv, gw rotations is supplied.


    Params:
        m = a positive integer input variable set to the number
            of rows of A.
        
        n = a positive integer input variable set to the number
            of columns of A.
        
        a = an array of length m^2. On input a must contain the matrix
            A to be postmultiplied by the orthogonal matrix Q
            described above. On output AQ has replaced A.
            the first min(m,n) columns of Q contains the factored form.
            on output Q has been accumulated into a square matrix.
        
        lda = a positive integer input variable not less than m
            which specifies the leading dimension of the array a.
        
        v = an input array of length n. v(i) must contain the
            information necessary to recover the Givens rotation gv(i)
            described above.

        v = an input array of length n. v(i) must contain the
            information necessary to recover the Givens rotation gw(i)
            described above.
*/
void r1mpyq(Real)(size_t m, size_t n, Real* a, size_t lda, Real* v, Real* w)
{
    size_t i, j, ij, in_, nmj;
    Real cos, sin, temp;
    
    enum  Real one = 1.0;

    
    // Apply the first set of Givens rotations to A.
    if (n > 1)
    {
        size_t nm1 = n - 1;
        for (nmj=2; nmj<=n; nmj++)
        {
            j = n - nmj;
            if (abs(v[j]) > one)
            {
                cos = one/v[j];
                sin = sqrt(one-cos*cos);
            }
            else
            {
                sin = v[j];
                cos = sqrt(one-sin*sin);
            }

            for (i=0; i<m; i++)
            {
                ij = i + j*lda;
                in_ = i + nm1*lda;
                temp = cos*a[ij] - sin*a[in_];
                a[in_] = sin*a[ij] + cos*a[in_];
                a[ij] = temp;
            }
        }


        // Apply the second set of Givens rotations to A.
        for (j=0; j<nm1; j++)
        {
            if (abs(w[j]) > one)
            {
                cos = one/w[j];
                sin = sqrt(one-cos*cos);
            }
            else
            {
                sin = w[j];
                cos = sqrt(one-sin*sin);
            }

            for (i=0; i<m; i++)
            {
                ij = i + j*lda;
                in_ = i + nm1*lda;
                temp = cos*a[ij] + sin*a[in_];
                a[in_] = -sin*a[ij] + cos*a[in_];
                a[ij] = temp;
            }
        }
    } 
}


unittest { alias r1mpyq!(real) rr1mpyq; }
