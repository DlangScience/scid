// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/napack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.napack.quasi;


import scid.core.fortran;

import scid.ports.napack.addchg;
import scid.ports.napack.stopit;




//
//      ________________________________________________________
//     |                                                        |
//     |SOLVE  F SUB I (X) = 0, I = 1 TO N, USING A QUASI-NEWTON|
//     |              SCHEME (BROYDEN'S GOOD METHOD)            |
//     |                                                        |
//     |    INPUT:                                              |
//     |                                                        |
//     |         X     --STARTING GUESS                         |
//     |                                                        |
//     |         H     --STARTING GUESS FOR INVERSE JACOBIAN    |
//     |                                                        |
//     |         LH    --LEADING (ROW) DIMENSION OF ARRAY H     |
//     |                                                        |
//     |         N     --NUMBER OF EQUATIONS                    |
//     |                                                        |
//     |         NDIGIT--DESIRED NUMBER CORRECT DIGITS          |
//     |                                                        |
//     |         LIMIT --MAXIMUM NUMBER OF ITERATIONS           |
//     |                                                        |
//     |         FUNC  --NAME OF FUNCTION EVALUATION SUBROUTINE |
//     |                 (EXTERNAL IN MAIN PROG.) FUNC(F,X) PUTS|
//     |                 IN ARRAY F THE FUNCTION VALUE AT X     |
//     |                                                        |
//     |         W     -WORK ARRAY (LENGTH AT LEAST 3N)         |
//     |                                                        |
//     |    OUTPUT:                                             |
//     |                                                        |
//     |         X     --SOLUTION                               |
//     |                                                        |
//     |         DIF   --INPUT FOR SUBROUTINE WHATIS            |
//     |                                                        |
//     |         SIZE  --INPUT FOR SUBROUTINE WHATIS            |
//     |                                                        |
//     |    PACKAGE SUBROUTINES: ADDCHG,STOPIT                  |
//     |________________________________________________________|
///
void quasi(Real, Func)(Real[] x_, Real[] h_, int lh, int n, out Real dif,
    out Real size, int ndigit, int limit, Func func, Real[] w_)
{
    int i, j;
    Real s, t, u;
    auto x = dimension(x_.ptr, n);
    auto h = dimension(h_.ptr, lh, n);
    auto w = dimension(w_.ptr, n, 3);
    for (i=1; i<=n; i++)
    {
        w[i,1] = 0.0;
    }
    goto l30;
//     --------------------
//     |*** ADD S TO X ***|
//     --------------------
l20:addchg(dif, size, x_, w_, n);
    stopit(dif, size, ndigit, limit);
    if (dif <= 0.0) return;
//     -----------------
//     |*** FORM  U ***|
//     -----------------
l30:for (i=1; i<=n; i++)
    {
        w[i,2] = w[i,1];
        w[i,1] = 0.0;
    }
    func(w_[2*lh .. 3*lh], x_);
    for (j=1; j<=n; j++)
    {
        t = w[j,3];
        for (i=1; i<=n; i++)
        {
            w[i,1] = w[i,1] - t*h[i,j];
        }
    }
//     -----------------------------------
//     |*** FORM  S DOT S AND S DOT U ***|
//     -----------------------------------
    t = 0.0;
    u = 0.0;
    for (i=1; i<=n; i++)
    {
        s = w[i,2];
        t = t + s*w[i,1];
        u = u + s*s;
    }
    t = u - t;
    if (t == 0.0) goto l20;
//     ------------------
//     |*** UPDATE H ***|
//     ------------------
    for (j=1; j<=n; j++)
    {
        s = 0.0;
        for (i=1; i<=n; i++)
        {
            s += h[i,j]*w[i,2];
        }
        s /= t;
        for (i=1; i<=n; i++)
        {
            h[i,j] = h[i,j] + s*w[i,1];
        }
    }
    u /= t;
    for (i=1; i<=n; i++)
    {
        w[i,1] = u*w[i,1];
    }
    goto l20;
}


unittest
{
    alias quasi!(float, void delegate(float[], float[])) fquasi;
    alias quasi!(double, void delegate(double[], double[])) dquasi;
    alias quasi!(double, void function(double[], double[])) dfquasi;
    alias quasi!(real, void delegate(real[], real[])) rquasi;
}


version(unittest)
{
    import std.math;
    import scid.core.testing;
}

unittest
{
    void rosenbrock(double[] fval, double[] x)
    {
        fval[0] = (1 - x[0]);
        fval[1] = 10.0 * (x[1] - x[0]^^2);
    }

    double[] guess = [-10.0, -5.0];
    double[] invJacAtGuess = [-1.0, 20.0, 0.0, 0.1];
    double[] work = new double[6];

    double dif, size;
    quasi(guess, invJacAtGuess, 2, 2, dif, size, 6, 500, &rosenbrock, work);
    check (guess == [1.0, 1.0]);
}
