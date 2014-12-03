// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/minpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0
*/
module scid.ports.minpack.fdjac1;


private import std.algorithm: max;
private import std.math: abs, sqrt;



/** This subroutine computes a forward-difference approximation
    to the n by n Jacobian matrix associated with a specified
    problem of n functions in n variables. If the Jacobian has
    a banded form, then function evaluations are saved by only
    approximating the nonzero terms.


    Params:
        fcn = 
            the name of the user-supplied function or delegate which
            calculates the functions. fcn should be written as follows.
            ---
            void f(size_t n, real* x, real* fvec, int iflag)
            {
                // calculate the functions at x and
                // return this vector in fvec.
            }
            ---
            The value of iflag should not be changed by fcn unless
            the user wants to terminate execution of fdjac1.
            In this case set iflag to a negative integer.

        n = 
            a positive integer input variable set to the number
            of functions and variables.

        x = 
            an input array of length n.

        fvec = 
            an input array of length n which must contain the
            functions evaluated at x.

        fjac =
            an output n by n array which contains the
            approximation to the Jacobian matrix evaluated at x.

        ldfjac =
            a positive integer input variable not less than n
            which specifies the leading dimension of the array fjac.

        iflag =
            an integer variable which can be used to terminate
            the execution of fdjac1. See description of fcn.

        ml = 
            a nonnegative integer input variable which specifies
            the number of subdiagonals within the band of the
            Jacobian matrix. If the Jacobian is not banded, set
            ml to at least n - 1.

        epsfcn =
            an input variable used in determining a suitable
            step length for the forward-difference approximation. This
            approximation assumes that the relative errors in the
            functions are of the order of epsfcn. If epsfcn is less
            than the machine precision, it is assumed that the relative
            errors in the functions are of the order of the machine
            precision.

        mu =
            a nonnegative integer input variable which specifies
            the number of superdiagonals within the band of the
            Jacobian matrix. If the jacobian is not banded, set
            mu to at least n - 1.

        wa1 =
            work array of length n.

        wa2 =
            work array of length n. If ml + mu + 1 is at
            least n, then the Jacobian is considered dense, and wa2 is
            not referenced.
*/
void fdjac1(Real, Func)(Func fcn, size_t n, Real* x, Real* fvec,
    Real* fjac, size_t ldfjac, int iflag, size_t ml, Real epsfcn,
    size_t mu, Real* wa1, Real* wa2)
{
    size_t i, j, k;
    Real temp, h;


    enum Real zero = 0.0;

    // epsmch is the machine precision.
    enum Real epsmch = Real.epsilon;


    immutable Real eps = sqrt(max(epsfcn, epsmch));
    immutable size_t msum = ml + mu + 1;

    if (msum >= n)
    {
        // Computation of dense approximate Jacobian.
        for (j=0; j<n; j++)
        {
            temp = x[j];
            h = eps*abs(temp);
            if (h == zero) h = eps;
            x[j] = temp + h;

            fcn(n, x, wa1, iflag);
            if (iflag < 0)  return;

            x[j] = temp;
            for (i=0; i<n; i++)
                fjac[i+j*ldfjac] = (wa1[i] - fvec[i])/h;
        }
        return;
    }


    // Computation of banded approximate Jacobian.
    size_t ij;
    for (k=0; k<msum; k++)
    {
        for (j=k; j<n; j+= msum)
        {
            wa2[j] = x[j];
            h = eps*abs(wa2[j]);
            if (h == zero) h = eps;
            x[j] = wa2[j] + h;
        }

        fcn(n, x, wa1, iflag);
        if (iflag < 0)  return;

        for (j=k; j<n; j+=msum)
        {
            x[j] = wa2[j];
            h = eps*abs(wa2[j]);
            if (h == zero) h = eps;
            
            for (i=0; i<n; i++)
            {
                ij = i+j*ldfjac;
                fjac[ij] = zero;
                if (i >= j-mu  &&  i <= j+ml)
                    fjac[ij] = (wa1[i] - fvec[i])/h;
            }
        }
    }
}


unittest
{
    void f(size_t n, double* a, double* fvec, ref int iflag)
    {
        auto x = a[0], y = a[1];

        fvec[0] = 2.0*x*y;
        fvec[1] = x-y;
    }

    bool close(double x, double y)
    {
        return (abs(x-y) < 1e-6);
    }

    double[] x = [ 1.0, 2.0 ];
    double[] fx = new double[2];
    double[] j = new double[4];
    double[] w1 = new double[2];
    double[] w2 = new double[2];
    int iflag = 1;

    f(2, x.ptr, fx.ptr, iflag);
    fdjac1(&f, 2, x.ptr, fx.ptr, j.ptr, 2, iflag, 2, 1e-8, 2, w1.ptr, w2.ptr);

    assert (close(j[0], 4.0) && close(j[2], 2.0)
        &&  close(j[1], 1.0) && close(j[3], -1.0));
}
