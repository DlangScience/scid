/** Functions related to the solving of nonlinear equations, i.e. finding
    roots of nonlinear functions.

    See_also:
    $(LINK2 http://www.digitalmars.com/d/2.0/phobos/std_numeric.html#findRoot,std.numeric.findRoot()),
    for solving nonlinear equations in one variable.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009-2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.nonlinear;


import std.traits;
import std.typecons;

import scid.core.memory;
import scid.core.traits;
import scid.ports.minpack.hybrd;
import scid.ports.napack.quasi;
import scid.calculus;
import scid.exception;
import scid.linalg;
import scid.util;

version (unittest)
{
    import std.math;
    import scid.core.testing;
}




/** Searches for a root of N functions of N variables using a variant
    of the Powell Hybrid Method (the HYBRD routine from MINPACK).

    Params:
        f = The set of equations, given as a function, delegate or
            functor that takes an array of length N as input and
            returns an array of length N. If the function has
            an additional array input parameter this will be
            assumed to be a buffer for the output value.
        guess = A starting point for the algorithm. The closer this
            guess is to the true root, the greater chance that the
            algorithm converges.
        epsRel = Success criterion: The algorithm stops when
            the relative error between two consecutive iterations
            is at most epsRel.
        maxFuncEvals = (optional) The maximum number of function evaluations.
            If maxFuncEvals<1, it is set to 200*(N+1).
        buffer = (optional) A buffer of length at least N, for the return value.
            

    Example:
    The Rosenbrock function is a commonly used test problem for
    optimisation algorithms. It has a global minimum at (1,1) that is hard
    to locate numerically because it lies in a long, narrow valley.
    Instead of using an optimisation algorithm, let us try to locate
    the minimum by finding the root of the Rosenbrock function's
    gradient.
    ---
    // The Rosenbrock function is defined as
    //     f(x,y) = (1-x)^2 + 100 (y-x^2)^2.
    // Thus, its gradient is:
    real[] dRosenbrock(real[] v, real[] buf)
    {
        auto x = v[0], y = v[1];
        buf[0] = -2*(1-x) - 400*x*(y-x*x);
        buf[1] = 200 * (y-x*x);

        return buf;
    }

    real[] guess = [ 2.0, 2.0 ];
    auto root = findRoot(&dRosenbrock, guess, 0.0L);

    writeln(root);  // Prints "1 1". Yay!
    ---
*/
Real[] findRoot (Real, Func)
    (Func f, Real[] guess, Real epsRel, int maxFuncEvals = 0,
    Real[] buffer=null)
in
{
    assert (guess.length > 0, "findRoot: empty guess vector given");
}
body
{
    mixin (newFrame);

    static assert (isFloatingPoint!Real,
        "findRoot: Not a floating-point type: "~T.stringof);
    static assert (isVectorField!(Func, Real),
        "findRoot: Invalid function type ("~Func.stringof~"), or "
       ~"function type doesn't match parameter type ("~Real.stringof~")");

    // Wrap the user-supplied function.
    void fcn(size_t m, Real* x, Real* fvec, ref int iflag)
    {
        static if (isBufferVectorField!(Func, Real))
        {
            // When the function takes a buffer, we check that it
            // actually uses it.
            auto fTest = f(x[0 .. m], fvec[0 .. m]);
            assert (fTest.ptr == fvec);
        }
        else
        {
            auto fTest = f(x[0 .. m]);
        }
        assert (fTest.length == m, "findRoot: The number of "
            ~"equations must be equal to the number of variables");
    };

    immutable int n = guess.length;
    immutable int wslen = (n*(3*n + 15))/2;
    if (maxFuncEvals < 1) maxFuncEvals = 200*(n+1);

    // Copy the guessed vector into the buffer.
    buffer.length = n;
    buffer[] = guess[];

    // There are a lot of parameters to the hybrd function, and we set them
    // in the "correct" order and use the "correct" names.
    Real* x = buffer.ptr;
    Real* fvec = cast(Real*) TempAlloc.malloc(wslen);
    alias epsRel xtol;
    alias maxFuncEvals maxfev;
    size_t ml_mu = n-1;
    Real epsfcn = 0.0;
    Real* diag = fvec + n;  diag[0 .. n] = 1.0;
    int mode = 2;
    enum Real factor = 100.0;
    int nprint = 0;
    int info = 0;
    uint nfev;
    Real* fjac = diag + n;
    alias n ldfjac;
    Real* r = fjac + ldfjac*n;
    size_t lr = (n*(n+1))/2;
    Real* qtf = r + lr;
    Real* wa1 = qtf + n;
    Real* wa2 = wa1 + n;
    Real* wa3 = wa2 + n;
    Real* wa4 = wa3 + n;
    
    // Phew! Call hybrd() now.
    hybrd!(Real, typeof(&fcn))(&fcn, n, x, fvec, xtol, maxfev, ml_mu, ml_mu,
        epsfcn, diag, mode, factor, nprint, info, nfev, fjac, ldfjac, r, lr,
        qtf, wa1, wa2, wa3, wa4);

    switch (info)
    {
        case 1: // Success!
            return x[0 .. n];
        
        case 0:
            throw new NumericsException(NE.InvalidInput);
        case 2:
            throw new NumericsException(NE.Limit);
        case 3:
            throw new NumericsException(NE.Accuracy);
        case 4:
        case 5:
            throw new NumericsException(NE.Convergence);
        default:
            throw new NumericsException;
    }
}

unittest
{
    real[] dRosenbrock(real[] v, real[] fx=null)
    {
        assert (v.length == 2 && fx.length == 2);
        auto x = v[0], y = v[1];
        fx[0] = -2*(1-x) - 400*x*(y-x*x);
        fx[1] = 200 * (y - x*x);
        return fx;
    }

    real[] guess = [ 2.0, 2.0 ];
    auto root = findRoot(&dRosenbrock, guess, 0.0L);
    check (approxEqual(root, [1.0L, 1.0L].dup, 1e-6));
}




/** Divides the interval into the given number of equal-sized subintervals,
    checks whether any of the subintervals bracket a root, and returns
    the ones that do.  If an endpoint of a subinterval [a,b] $(I is) a root,
    i.e. f(a)=0, then the interval is returned as [a,a].

    A buffer of length at least nIntervals+1, for storing the brackets, may
    optionally be provided.  If not, one will be allocated.

    This function is primarily meant for use with
    $(LINK2 http://www.digitalmars.com/d/2.0/phobos/std_numeric.html#findRoot,std.numeric.findRoot()).
*/
T[2][] bracketRoots(T, Func)(Func f, T lower, T upper, uint nIntervals,
    T[2][] buffer = null)
{
    buffer.length = nIntervals+1;
    int numBrackets = 0;

    auto lo = lower;
    auto flo = f(lo);
    immutable step = (upper - lower)/nIntervals;

    foreach (i; 0 .. nIntervals)
    {
        immutable hi = (i < nIntervals - 1 ? lo + step : upper);
        immutable fhi = f(hi);

        if (flo == 0)
        {
            T[2] b;
            b[0] = lo;
            b[1] = lo;
            buffer[numBrackets] = b;
            ++numBrackets;
        }
        else
        {
            if (flo * fhi < 0)
            {
                T[2] b;
                b[0] = lo;
                b[1] = hi;
                buffer[numBrackets] = b;
                ++numBrackets;
            }
        }

        lo = hi;
        flo = fhi;
    }

    // Check for a root in the endpoint as well.
    if (flo == 0)
    {
        T[2] b;
        b[0] = lo;
        b[1] = lo;
        buffer[numBrackets] = b;
        ++numBrackets;
    }
    
    buffer.length = numBrackets;
    return buffer;
}


unittest
{
    real f(real x)
    {
        return (2+x) * (1+x) * x * (1-x) * (2-x);
    }

    bool has(real[2] intv, real x)
    {
        return intv[0] < x && intv[1] > x;
    }

    auto b = bracketRoots(&f, -2.0L, 2.0L, 15);
    check (b.length == 5);
    check (b[0][0] == -2 && b[0][1] == -2);
    check (has(b[1], -1));
    check (has(b[2], 0));
    check (has(b[3], 1));
    check (b[4][0] ==  2 && b[4][1] ==  2);
}




/** Given a function f and a starting point x0, this routine searches
    along the x-axis in the positive (scale>0) or negative (scale<0)
    direction until f(x) has an opposite sign from f(x0).
    (It first tries x0+scale, then, for each iteration, scale is
    multiplied by a constant factor.)

    If such a point is found, bracketFrom() returns a tuple containing
    x and f(x).  If not, an exception is thrown.
*/
Tuple!(T, "x", T, "fx") bracketFrom
    (T, Func)
    (Func f, T x0, T scale, int maxIterations=10)
in
{
    assert (scale != 0);
    assert (maxIterations > 0);
}
body
{
    immutable fx0 = f(x0);

    enum expandFactor = 1.6;
    real step = scale;
    foreach (i; 0 .. maxIterations)
    {
        immutable x = x0 + step;
        immutable fx = f(x);
        if (fx0 * fx < 0) return typeof(return)(x, fx);
        step *= expandFactor;
    }

    enforceNE(false, NE.Limit);
    assert(0);
}


unittest
{
    real f(real x) { return x^^2 - 100; }
    
    auto upper = bracketFrom(&f, 0.0L, 1.0L);
    check(f(upper.x) == upper.fx);
    check(upper.fx > 0);
}
