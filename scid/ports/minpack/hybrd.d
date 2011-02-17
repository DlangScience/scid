/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost Licence 1.0 
*/
module scid.ports.minpack.hybrd;


private import std.algorithm: min, max;
private import std.math: abs;

private import scid.ports.minpack.dogleg;
private import scid.ports.minpack.enorm;
private import scid.ports.minpack.fdjac1;
private import scid.ports.minpack.qform;
private import scid.ports.minpack.qrfac;
private import scid.ports.minpack.r1mpyq;
private import scid.ports.minpack.r1updt;




/** The purpose of hybrd is to find a zero of a system of
    n nonlinear functions in n variables by a modification
    of the Powell Hybrid Method. The user must provide a
    function which calculates the functions. The Jacobian is
    then calculated by a forward-difference approximation.

    Params:
        fcn = the user-supplied function or delegate which
            calculates the functions. The function referenced
            by fcn must be declared
            in an external statement in the user calling
            program, and should be written as follows.
            ---
            void foo (size_t n,real* x, real* fvec, int iflag)
            {
                // calculate the functions at x and
                // return this vector in fvec.
            }
            ---
            The value of iflag should not be changed by the function unless
            the user wants to terminate execution of hybrd.
            In this case set iflag to a negative integer.

        n = a positive integer input variable set to the number
            of functions and variables.

        x = an array of length n. On input x must contain
            an initial estimate of the solution vector. On output x
            contains the final estimate of the solution vector.

        fvec = an output array of length n which contains
            the functions evaluated at the output x.

        xtol = a nonnegative input variable. Termination
            occurs when the relative error between two consecutive
            iterates is at most xtol.

        maxfev = a positive integer input variable. Termination
            occurs when the number of calls to fcn is at least maxfev
            by the end of an iteration.

        ml = a nonnegative integer input variable which specifies
            the number of subdiagonals within the band of the
            jacobian matrix. If the Jacobian is not banded, set
            ml to at least n - 1.

        mu = a nonnegative integer input variable which specifies
            the number of superdiagonals within the band of the
            jacobian matrix. If the Jacobian is not banded, set
            mu to at least n - 1.

        epsfcn = an input variable used in determining a suitable
            step length for the forward-difference approximation. This
            approximation assumes that the relative errors in the
            functions are of the order of epsfcn. If epsfcn is less
            than the machine precision, it is assumed that the relative
            errors in the functions are of the order of the machine
            precision.

        diag = an array of length n. If mode = 1 (see
            below), diag is internally set. If mode = 2, diag
            must contain positive entries that serve as
            multiplicative scale factors for the variables.

        mode = an integer input variable. If mode = 1, the
            variables will be scaled internally. If mode = 2,
            the scaling is specified by the input diag. Other
            values of mode are equivalent to mode = 1.

        factor = a positive input variable used in determining the
            initial step bound. This bound is set to the product of
            factor and the Euclidean norm of diag*x if nonzero, or else
            to factor itself. In most cases factor should lie in the
            interval (.1,100.). 100. is a generally recommended value.

        nprint = an integer input variable that enables controlled
            printing of iterates if it is positive. In this case,
            fcn is called with iflag = 0 at the beginning of the first
            iteration and every nprint iterations thereafter and
            immediately prior to return, with x and fvec available
            for printing. If nprint is not positive, no special calls
            of fcn with iflag = 0 are made.

        info = an integer output variable. If the user has
            terminated execution, info is set to the (negative)
            value of iflag. See description of fcn. Otherwise,
            info is set as follows.
            ---
            0:  improper input parameters.

            1:  relative error between two consecutive iterates
                is at most xtol.

            2:  number of calls to fcn has reached or exceeded
                maxfev.

            3:  xtol is too small. no further improvement in
                the approximate solution x is possible.

            4:  iteration is not making good progress, as
                measured by the improvement from the last
                five Jacobian evaluations.

            5:  iteration is not making good progress, as
                measured by the improvement from the last
                ten iterations.
            ---

        nfev = an integer output variable set to the number of
            calls to fcn.

        fjac = an output n by n array which contains the
            orthogonal matrix Q produced by the QR factorization
            of the final approximate Jacobian.

        ldfjac = a positive integer input variable not less than n
            which specifies the leading dimension of the array fjac.

        r = an output array of length lr which contains the
            upper triangular matrix produced by the QR factorization
            of the final approximate Jacobian, stored rowwise.

        lr = a positive integer input variable not less than
            (n*(n+1))/2.

        qtf = an output array of length n which contains
            the vector (Q transpose)*fvec.

        wa1 = work array of length n.
        wa2 = work array of length n.
        wa3 = work array of length n.
        wa4 = work array of length n.
*/
void hybrd(Real, Func)(Func fcn, size_t n, Real* x, Real* fvec,
    Real xtol, uint maxfev, size_t ml, size_t mu, Real epsfcn,
    Real* diag, int mode, Real factor, int nprint, out int info,
    out uint nfev, Real* fjac, size_t ldfjac, Real* r, size_t lr,
    Real* qtf, Real* wa1, Real* wa2, Real* wa3, Real* wa4)
{
    size_t i, j, l;
    ptrdiff_t jm1;
    long msum;
    int iflag;
    size_t[1] iwa;
    bool jeval, sing;
    Real actred, delta, fnorm, fnorm1, pnorm, prered, ratio, sum, temp, xnorm;

    enum : Real
    {
        one   = 1.0,
        p1    = 1e-1,
        p5    = 5e-1,
        p001  = 1e-3,
        p0001 = 1e-4,
        zero  = 0.0,

        // epsmch is the machine precision.
        epsmch  = Real.epsilon
    }

    info = 0;
    iflag = 0;
    nfev = 0;


    // Check the input parameters for errors.
    if (xtol < zero  ||  factor <= zero  ||  ldfjac < n  ||  lr < (n*(n+1)/2))
        goto end;

    if (mode == 2)
    {
        for (j=0; j<n; j++)
            if (diag[j] <= zero) goto end;
    }


    // Evaluate the function at the starting point
    // and calculate its norm.
    iflag = 1;
    fcn(n, x, fvec, iflag);
    nfev = 1;
    if (iflag < 0) goto end;
    fnorm = enorm(n, fvec);


    // Determine the number of calls to fcn needed to compute
    // the Jacobian matrix.
    msum = min(ml+mu+1, n);


    // Initialize iteration counter and monitors.
    int iter = 1;
    int nsuc = 0;
    int ncfail = 0;
    int nslow1 = 0;
    int nslow2 = 0;


    // Beginning of the outer loop.
    outer: while(true)
    {
        jeval = true;

        // Calculate the Jacobian matrix.
        fdjac1!(Real, Func)
            (fcn, n, x, fvec, fjac, ldfjac, iflag, ml, epsfcn, mu, wa1, wa2);
        nfev += msum;
        if (iflag < 0) break outer;


        // Compute the QR factorization of the Jacobian.
        qrfac(n, n, fjac, ldfjac, false, iwa.ptr, 1, wa1, wa2, wa3);


        // On the first iteration and if mode is 1, scale according
        // to the norms of the columns of the initial Jacobian.
        if (iter == 1)
        {
            if (mode != 2)
            {
                for (j=0; j<n; j++)
                {
                    diag[j] = wa2[j];
                    if (wa2[j] == zero) diag[j] = one;
                }
            }


            // On the first iteration, calculate the norm of the scaled x
            // and initialize the step bound delta.
            for (j=0; j<n; j++)
                wa3[j] = diag[j]*x[j];
            xnorm = enorm(n, wa3);
            delta = factor*xnorm;
            if (delta == zero) delta = factor;
        }


        // Form (Q transpose)*fvec and store in qtf.
        for (i=0; i<n; i++)
            qtf[i] = fvec[i];
        for (j=0; j<n; j++)
        {
            if (fjac[j+j*ldfjac] != zero)
            {
                sum = zero;
                for (i=j; i<n; i++)
                    sum += fjac[i+j*ldfjac]*qtf[i];
                temp = -sum/fjac[j+j*ldfjac];
                for (i=j; i<n; i++)
                    qtf[i] = qtf[i] + fjac[i+j*ldfjac]*temp;
            }
        }


        // Copy the triangular factor of the QR factorization into r.
        sing = false;
        for (j=0; j<n; j++)
        {
            l = j;
            jm1 = j - 1;
            if (jm1 >= 0)
            {
                for (i=0; i<=jm1; i++)
                {
                    r[l] = fjac[i+j*ldfjac];
                    l += n - i - 1;
                }
            }

            r[l] = wa1[j];
            if (wa1[j] == zero) sing = true;
        }


        // Accumulate the orthogonal factor in fjac.
        qform!(Real)(n, n, fjac, ldfjac, wa1);


        // Rescale if necessary.
        if (mode != 2)
        {
            for (j=0; j<n; j++)
                diag[j] = max(diag[j], wa2[j]);
        }
        

        // Beginning of inner loop.
        inner: while(true)
        {
            // If requested, call fcn to enable printing of iterates.
            if (nprint > 0)
            {
                iflag = 0;
                if ((iter-1 % nprint) == 0)  fcn(n, x, fvec, iflag);
                if (iflag < 0)  break outer;
            }


            // Determine the direction p.
            dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3);


            // Store the direction p and x + p. Calculate the norm of p.
            for (j=0; j<n; j++)
            {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j]*wa1[j];
            }
            pnorm = enorm(n, wa3);


            // On the first iteration, adjust the initial step bound.
            if (iter == 1)  delta = min(delta, pnorm);


            // Evaluate the function at x + p and calculate its norm.
            iflag = 1;
            fcn(n, wa2, wa4, iflag);
            nfev++;
            if (iflag < 0)  break outer;
            fnorm1 = enorm(n, wa4);


            // Compute the scaled actual reduction.
            actred = -one;
            if (fnorm1 < fnorm)  actred = one - (fnorm1/fnorm)^^2;


            // Compute the scaled predicted reduction.
            l = 0;
            for (i=0; i<n; i++)
            {
                sum = zero;
                for (j=i; j<n; j++)
                {
                    sum += r[l]*wa1[j];
                    l++;
                }
                wa3[i] = qtf[i] + sum;
            }

            temp = enorm(n, wa3);
            prered = zero;
            if (temp < fnorm)  prered = one - (temp/fnorm)^^2;


            // Compute the ratio of the actual to the predicted
            // reduction.
            ratio = zero;
            if (prered > zero)  ratio = actred/prered;


            // Update the step bound.
            if (ratio < p1)
            {
                nsuc = 0;
                ncfail++;
                delta *= p5;
            }
            else
            {
                ncfail = 0;
                nsuc++;
                if (ratio >= p5  ||  nsuc > 1)
                    delta = max(delta, pnorm/p5);
                if (abs(ratio-one) <= p1)  delta = pnorm/p5;
            }


            // Test for successful iteration.
            if (ratio >= p0001)
            {
                // Successful iteration. Update x, fvec, and their norms.
                for (j=0; j<n; j++)
                {
                    x[j] = wa2[j];
                    wa2[j] = diag[j]*x[j];
                    fvec[j] = wa4[j];
                }

                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }


            // Determine the progress of the iteration.
            nslow1++;
            if (actred >= p001) nslow1 = 0;
            if (jeval) nslow2++;
            if (actred >= p1) nslow2 = 0;


            // Test for convergence.
            if (delta <= xtol*xnorm  ||  fnorm == zero)
            {
                info = 1;
                break outer;
            }


            // Tests for termination and stringent tolerances.
            if (nfev >= maxfev)  info = 2;
            if (p1*max(p1*delta, pnorm) <= epsmch*xnorm)  info = 3;
            if (nslow2 == 5)  info = 4;
            if (nslow1 == 10)  info = 5;
            if (info != 0)  break outer;


            // Criterion for recalculating Jacobian approximation
            // by forward differences.
            if (ncfail == 2)  break inner;
            
            // Calculate the rank one modification to the Jacobian
            // and update qtf if necessary.
            for (j=0; j<n; j++)
            {
                sum = zero;
                for (i=0; i<n; i++)
                    sum += fjac[i+j*ldfjac]*wa4[i];
                wa2[j] = (sum - wa3[j])/pnorm;
                wa1[j] = diag[j]*((diag[j]*wa1[j])/pnorm);
                if (ratio >= p0001)  qtf[j] = sum;
            }


            // Compute the QR factorization of the updated Jacobian.
            r1updt(n, n, r, lr, wa1, wa2, wa3, sing);
            r1mpyq(n, n, fjac, ldfjac, wa2, wa3);
            r1mpyq(1, n, qtf, 1, wa2, wa3);


            // End of the inner loop.
            jeval = false;
        }


        // End of the outer loop.
    }

end:
    // Termination, either normal or user imposed.
    if (iflag < 0) info = iflag;

    iflag = 0;
    if (nprint > 0) fcn(n, x, fvec, iflag);
}


version (unittest) { import std.stdio; }

unittest
{
    real a = 1.0;
    real b = 10.0;

    void rosenbrock(size_t n, real* x, real* fx, ref int iflag)
    {
        assert (n == 2);
        fx[0] = a * (1 - x[0]);
        fx[1] = b * (x[1] - x[0]*x[0]);
        //writefln("f(%s, %s) = (%s, %s)", x[0], x[1], fx[0], fx[1]);
    }

    alias hybrd!(real, typeof(&rosenbrock)) rhybrd;

    enum size_t n = 2;
    enum size_t lr = (n*(n+1))/2;

    real[n] x = [ -10.0, -5.0 ];
    real[n] fvec, diag, qtf, wa1, wa2, wa3, wa4;
    real[lr] r;
    real[n*n] fjac;

    real xtol = 0.0;
    uint maxfev = 1000;
    real epsfcn = 1e-10;

    uint mode = 1;
    real factor = 100.0;
    int nprint = 0;
    int info;
    uint nfev;

    rhybrd(&rosenbrock, n, x.ptr, fvec.ptr, xtol, maxfev, n, n, epsfcn,
        diag.ptr, mode, factor, nprint, info, nfev, fjac.ptr, n, r.ptr,
        lr, qtf.ptr, wa1.ptr, wa2.ptr, wa3.ptr, wa4.ptr);
    //writefln("Status: %s  Result: (%s , %s)  nfev: %s", info, x[0], x[1], nfev);
    assert (info == 1);
    assert (abs(1.0-x[0]) <= xtol  &&  abs(1.0-x[1]) <= xtol);

}

