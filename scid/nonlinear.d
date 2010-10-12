/** Functions related to the solving of nonlinear equations, i.e. finding
    roots of nonlinear functions.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009-2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.nonlinear;


import std.range;
import std.traits;

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




/** An interval bracketing a root of some function.

    If a function f(x) is continuous on an interval (x1,x2),
    and f(x1) and f(x2) have opposite sign, we know the function
    must pass through zero somewhere in the interval.
    Such an interval is said to 'bracket' a root of the function.
*/
struct BracketingInterval(X, Y)
{
    /// The interval limits.
    X x1;
    X x2;   /// ditto

    /// The function value at x1 and x2, respectively
    Y y1;
    Y y2;   /// ditto
}




/** Find a root of the function f, given an interval that
    is known to bracket a root.

    This function is included for convenient use with
    BracketingIntervals returned by the bracketing functions
    in this module.  Under the hood it just forwards to the
    $(LINK2 http://www.digitalmars.com/d/2.0/phobos/std_numeric.html#findRoot,std.numeric.findRoot())
    function in Phobos.
*/
T findRoot(T, R)(R delegate(T) f, BracketingInterval!(T, R) bracket)
{
    return std.numeric.findRoot(f,
        bracket.x1, bracket.x2,
        bracket.y1, bracket.y2,
        (T a, T b) { return false; }    // Machine precision
        ).field[0];
}




/** Uses bracketRoots() to divide the interval [a,b] into subintervals
    and check which ones bracket roots.  Then, findRoot() is applied
    to each bracketing interval, and an array containing the roots
    is returned.

    A buffer of length at least nIntervals+1, for storing the roots,
    may optionally be provided.
*/
T[] findRoots(T, Func)(Func f, T a, T b, uint nIntervals,
    T[] buffer=null)
{
    mixin(scid.core.memory.newFrame);

    // Find bracketing subintervals.
    auto bracketBuffer =
        newStack!(BracketingInterval!(T, ReturnType!Func))(nIntervals+1);
    auto intervals =
        bracketRoots!(T,Func)(f, a, b, nIntervals, bracketBuffer);

    // Find all the bracketed roots.
    buffer.length = intervals.length;
    foreach (i, iv; intervals)
    {
        // If both endpoints are equal, this value *is* the root.
        if (iv.x1 == iv.x2)  buffer[i] = iv.x1;

        // If not, call findRoot() to locate the root.
        else buffer[i] = findRoot(f, iv);
    }

    return buffer;
}


unittest
{
    real f(real x)
    {
        return (2+x) * (1+x) * x * (1-x) * (2-x);
    }

    auto r = findRoots(&f, -2.0L, 2.0L, 15);
    check (r.length == 5);
    check (approxEqual(r, [-2.0, -1.0, 0.0, 1.0, 2.0], real.epsilon));
}




/** Divides the interval [a,b] into the given number of equal-sized
    subintervals,
    checks whether any of the subintervals bracket a root, and returns
    the ones that do, together with the function values at those points.
    If an endpoint of a subinterval [x1,x2] $(I is) a root,
    i.e. f(x1)=0, then the interval is returned as [x1,x1].

    A buffer of length at least nIntervals+1, for storing the brackets, may
    optionally be provided.  If not, one will be allocated.
*/
BracketingInterval!(T, ReturnType!Func)[] bracketRoots(T, Func)
    (Func f, T a, T b, uint nIntervals,
     BracketingInterval!(T, ReturnType!Func)[] buffer = null)
{
    static assert (is (typeof(buffer) == typeof(return)));
    alias ElementType!(typeof(return)) B;

    buffer.length = nIntervals+1;
    int numBrackets = 0;

    auto lo = a;
    auto flo = f(lo);
    immutable step = (b - a)/nIntervals;

    foreach (i; 0 .. nIntervals)
    {
        immutable hi = (i < nIntervals - 1 ? lo + step : b);
        immutable fhi = f(hi);

        if (flo == 0)
        {
            B br;
            br.x1 = lo;  br.y1 = flo;
            br.x2 = lo;  br.y2 = flo;
            buffer[numBrackets] = br;
            ++numBrackets;
        }
        else
        {
            if (flo * fhi < 0)
            {
                B br;
                br.x1 = lo;  br.y1 = flo;
                br.x2 = hi;  br.y2 = fhi;
                buffer[numBrackets] = br;
                ++numBrackets;
            }
        }

        lo = hi;
        flo = fhi;
    }

    // Check for a root in the endpoint as well.
    if (flo == 0)
    {
        B br;
        br.x1 = lo;  br.y1 = flo;
        br.x2 = lo;  br.y2 = flo;
        buffer[numBrackets] = br;
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

    bool has(B)(B intv, real x)
    {
        return intv.x1 < x && intv.x2 > x;
    }

    auto b = bracketRoots(&f, -2.0L, 2.0L, 15);
    check (b.length == 5);
    check (b[0].x1 == -2 && b[0].x2 == -2);
    check (has(b[1], -1));
    check (has(b[2], 0));
    check (has(b[3], 1));
    check (b[3].y1 == f(b[3].x1)  &&  b[3].y2 == f(b[3].x2));
    check (b[4].x1 ==  2 && b[4].x2 ==  2);
}




/** Given a function f and a starting point x1, this routine searches
    along the x-axis in the positive (scale>0) or negative (scale<0)
    direction until it reaches a point x2 where f(x2) has an opposite
    sign from f(x1).

    If such a point is found, a BracketingInterval containing the two
    points x1 and x2, as well as the function values in those points,
    is returned.  If not, an exception is thrown.

    On the first iteration, x2 = x1+scale.  Hence, scale should be a
    characteristic scale for the function (i.e. a scale over which the
    function changes significantly).  Thereafter, the interval is expanded
    geometrically by multiplying scale by a constant factor for each iteration.
*/
BracketingInterval!(T, ReturnType!Func) bracketFrom(T, Func)
    (Func f, T x1, T scale, int maxIterations=40)
in
{
    assert (scale != 0);
    assert (maxIterations > 0);
}
body
{
    immutable fx1 = f(x1);

    enum expandFactor = 1.6;
    real step = scale;
    foreach (i; 0 .. maxIterations)
    {
        immutable x2 = x1 + step;
        immutable fx2 = f(x2);
        if (fx1 * fx2 < 0) return typeof(return)(x1, x2, fx1, fx2);
        step *= expandFactor;
    }

    enforceNE(false, NE.Limit);
    assert(0);
}


unittest
{
    real f(real x) { return x^^2 - 100; }
    
    auto bracket = bracketFrom(&f, 0.0L, 1.0L);
    check(bracket.x1 == 0);
    check(f(bracket.x1) == bracket.y1);
    check(f(bracket.x2) == bracket.y2);
    check(bracket.x2 > 10);
    check(bracket.y2 > 0);
}
