/** Functions related to the solving of nonlinear equations, i.e. finding
    roots of nonlinear functions.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009-2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.nonlinear;


import std.algorithm;
import std.functional;
import std.math;
import std.range;
import std.traits;
import std.typecons;

import scid.core.memory;
import scid.core.traits;
import scid.ports.minpack.hybrd;
import scid.ports.napack.quasi;
import scid.calculus;
import scid.exception;
//import scid.types;
import scid.util;

version (unittest) { import scid.core.testing; }




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
    (scope Func f, Real[] guess, Real epsRel, int maxFuncEvals = 0,
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




/** Find a root of the function f.

    This function first calls $(LINK2 #bracketRoot,bracketRoot) to
    obtain an interval inside which there must be a root, and then calls
    $(LINK2 http://www.digitalmars.com/d/2.0/phobos/std_numeric.html#findRoot,std.numeric.findRoot())
    to pin down the location of the root.

    The parameters x0, scale, xMin, and xMax are just passed on to
    $(LINK2 #bracketRoot,bracketRoot), and they are described in detail
    in its documentation.  In brief, x0 should be an estimate of the
    root's location, while scale should be a characteristic scale for
    the function, i.e. a distance over which the function changes
    significantly.  [xMin,xMax] is the interval inside which the algorithm
    is allowed to search.

    You may specify the desired (minimum) number of digits of precision
    in the answer.  If this is left out, the algorithm will attempt
    to achieve full machine precision.
*/
T findRoot(F, T)(scope F f, T x0, T scale, T xMin, T xMax, int precision)
    if (isUnaryFunction!(F, T) && isFloatingPoint!T)
{
    return findRootImpl(f, x0, scale, xMin, xMax,
        (T a, T b) { return matchDigits(a, b, precision); });
}

/// ditto
T findRoot(F, T)(scope F f, T x0, T scale, int precision)
    if (isUnaryFunction!(F, T) && isFloatingPoint!T)
{
    return findRoot(f, x0, scale, -T.infinity, T.infinity, precision);
}

/// ditto
T findRoot(F, T)
    (scope F f, T x0, T scale, T xMin = -T.infinity, T xMax = T.infinity)
    if (isUnaryFunction!(F, T) && isFloatingPoint!T)
{
    return findRootImpl(f, x0, scale, xMin, xMax,
        (T a, T b) { return false; });
}


// Implementation of findRoot()
private T findRootImpl(F, T)
    (
        scope F f,
        T x0, T scale,
        T xMin, T xMax,
        scope bool delegate(T, T) tolerance
    )
    if (isUnaryFunction!(F, T) && isFloatingPoint!T)
{
    auto bracket = bracketRoot(f, x0, scale, xMin, xMax);
    if (bracket.y1 == 0) return bracket.x1;
    if (bracket.y2 == 0) return bracket.x2;

    // std.numeric.findRoot() only takes a delegate
    static if (is (F == delegate))
        auto dg = f;
    else static if (isFunctionPointer!F)
        auto dg = toDelegate(f);
    else static if (isFunctor!F)
        scope ReturnType!F delegate(ParameterTypeTuple!F) dg = &f.opCall;

    return std.numeric.findRoot(dg,
        bracket.x1, bracket.x2,
        bracket.y1, bracket.y2,
        tolerance
    )[0];
}


unittest
{
    real f(real x) { return log(x); }

    immutable inaccurateRoot =
        findRoot(&f, 0.5L, 1.0L, real.epsilon, real.infinity, 2);
    check(matchDigits(inaccurateRoot, 1.0, 2));
    check(!matchDigits(inaccurateRoot, 1.0, 10));

    immutable accurateRoot =
        findRoot(&f, 0.5L, 1.0L, real.epsilon, real.infinity);
    check(accurateRoot == 1.0);
}

unittest
{
    // Function
    static real f(real x) { return x; }
    check(findRoot(&f, 1.0L, 1.0L) == 0.0L);
}

unittest
{
    // Functor
    struct Functor { real opCall(real x) { return x^^3; } }
    Functor g;
    check(findRoot(g, 1.0L, 1.0L) == 0.0L);
}




/** Bracket a root of the function f.

    If a function f(x) is continuous on an interval [x1,x2],
    and f(x1) and f(x2) have opposite sign, we know the function
    must pass through zero somewhere in the interval.
    The points x1 and x2 are then said to 'bracket' the
    root.  This is usually the first step in locating the root
    of a function.

    If this function succeeds, it returns a RootBracket containing
    the points x1 and x2, together with the function values f(x1)
    and f(x2).  If it fails, an exception is thrown.
    Note that this library considers the points to be bracketing
    a root also if the root is located exactly at x1 and/or x2,
    i.e. if f(x1)=0 and/or f(x2)=0.

    Details:

    This function will start by evaluating f(x) in the points
    x0 and x0+scale and see if those
    points bracket a root of the given function.  If not, the interval
    is expanded geometrically (i.e. the distance between the points
    is multiplied by a constant factor), always in the direction where
    f(x) is smallest, until the points bracket a root.

    You may optionally specify a limiting interval [xMin, xMax], and the
    algorithm will never search outside it.  This is useful,
    for instance, for functions that are only defined for certain
    values of x.  If you do specify such an interval, the
    initial point x0 must lie inside it.

    It is usually worthwhile analysing the behaviour of the function
    in order to find appropriate values for x0 and scale.
    The closer x0 is to the actual root, the fewer steps (i.e. the
    fewer evaluations of f) this algorithm will require to succeed.
    If scale is too large, and the function has several roots,
    there is a chance that it will just step across both roots
    and not find any of them.  On the other hand, if it is too small,
    it may again cause the algorithm to take more steps than would
    otherwise be necessary.
*/
RootBracket!(T, ReturnType!F) bracketRoot(F, T)
    (
        scope F f,
        in T x0, in T scale,
        in T xMin = -T.infinity, in T xMax = T.infinity,
    )
    if (isFloatingPoint!T && isUnaryFunction!(F, T))
in
{
    assert (scale != 0, "scale must be nonzero");
    assert (xMin < xMax, "xMin must be smaller than xMax");
    assert (xMin <= x0 && x0 <= xMax, "x0 must be in the interval [xMin,xMax]");
}
body
{
    alias typeof(return) B;
    enum expandFactor = 1.6;


    // Function that searches upwards from xMin
    B upwards(real x, real fx, real dx)
    {
        immutable fxMin = f(xMin);
        for (;;)
        {
            if (fxMin * fx <= 0) return B(xMin, x, fxMin, fx);
            enforceNE(x != xMax, "Unable to bracket root");
            
            dx *= expandFactor;
            x = min(x + dx, xMax);
            fx = f(x);
        }
        assert(0);
    }


    // Function that searches downwards from xMax
    B downwards(real x, real fx, real dx)
    {
        immutable fxMax = f(xMax);
        for (;;)
        {
            if (fxMax * fx <= 0) return B(x, xMax, fx, fxMax);
            enforceNE(x != xMin, "Unable to bracket root");
            
            dx *= expandFactor;
            x = max(x - dx, xMin);
            fx = f(x);
        }
        assert(0);
    }


    // These are the initial points
    real x1 = x0;
    real x2 = x0 + scale;

    // If x0 is either endpoint of the allowed interval, or if x2
    // falls outside the interval, search only in one direction.
    if (x1 == xMin || x2 <= xMin)
    {
        immutable x = min(x1 + abs(scale), xMax);
        return upwards(x, f(x), abs(scale));
    }
    if (x1 == xMax || x2 >= xMax)
    {
        immutable x = max(x1 - abs(scale), xMin);
        return downwards(x, f(x), abs(scale));
    }

    
    // Both x1 and x2 fall inside [xMin, xMax], so we use bidirectional search
    if (x1 > x2) swap(x1, x2);
    real fx1 = f(x1);
    real fx2 = f(x2);

    for (;;)
    {
        // Check whether interval brackets a root
        if (fx1 * fx2 <= 0)  return B(x1, x2, fx1, fx2);

        // Expand interval in the direction where f(x) is closest to zero.
        if (fabs(fx1) < fabs(fx2))
        {
            x1 += expandFactor * (x1 - x2);
            if (x1 <= xMin) return upwards(x2, fx2, x2-x1);
            fx1 = f(x1);
        }
        else
        {
            x2 += expandFactor * (x2 - x1);
            if (x2 >= xMax) return downwards(x1, fx1, x2-x1);
            fx2 = f(x2);
        }
    }
    assert(0);
}


unittest
{
    real f(real x) { return 1 - x; }
    auto b = bracketRoot(&f, -100.0L, 1.0L);
    check (b.contains(1));
}

unittest
{
    real f(real x) { return log(x); }
    auto b = bracketRoot(&f, 2*real.epsilon, 0.1L, real.epsilon, real.infinity);
    check (b.contains(1));
}




/** A set of points that bracket a root of some function. */
struct RootBracket(X, Y)
{
    /// Two points that bracket a root.
    X x1;
    X x2;   /// ditto

    /// The function value at x1 and x2, respectively
    Y y1;
    Y y2;   /// ditto


    version(unittest) private bool contains(X point)
    {
        if (x1 <= x2) return (point >= x1 && point <= x2);
        else          return (point >= x2 && point <= x1);
    }
}




/** Uses bracketRoots() to divide the interval [a,b] into subintervals
    and check which ones bracket roots.  Then, findRoot() is applied
    to each bracketing interval, and an array containing the roots
    is returned.

    A buffer of length at least nIntervals+1, for storing the roots,
    may optionally be provided.
*/
T[] findRoots(T, Func)(scope Func f, T a, T b, uint nIntervals,
    T[] buffer=null)
{
    mixin(scid.core.memory.newFrame);

    // Find bracketing subintervals.
    auto bracketBuffer =
        newStack!(RootBracket!(T, ReturnType!Func))(nIntervals+1);
    auto intervals =
        bracketRoots!(T,Func)(f, a, b, nIntervals, bracketBuffer);

    // Find all the bracketed roots.
    buffer.length = intervals.length;
    foreach (i, iv; intervals)
    {
        // Check if a root is located at the "lower" endpoint.
        if (iv.y1 == 0) buffer[i] = iv.x1;

        // If not, call findRoot() to locate the root.
        // Note that if it is located at the "higher" end point it will
        // be caught in the next iteration.
        else if (iv.y2 != 0) buffer[i] =
            std.numeric.findRoot(f, iv.x1, iv.x2, iv.y1, iv.y2,
                (T a, T b) { return false; })[0];
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

    A buffer of length at least nIntervals+1, for storing the brackets, may
    optionally be provided.  If not, one will be allocated.
*/
RootBracket!(T, ReturnType!Func)[] bracketRoots(T, Func)
    (scope Func f, T a, T b, uint nIntervals,
     RootBracket!(T, ReturnType!Func)[] buffer = null)
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

        if (flo == 0 || flo * fhi < 0)
        {
            B br;
            br.x1 = lo;  br.y1 = flo;
            br.x2 = hi;  br.y2 = fhi;
            buffer[numBrackets] = br;
            ++numBrackets;
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

    auto b = bracketRoots(&f, -2.0L, 2.0L, 15);
    check(b.length == 5);
    foreach (i; b)
    {
        check(i.y1 == f(i.x1));
        check(i.y2 == f(i.x2));
    }
    check(b[0].x1 == -2);
    check(b[1].contains(-1));
    check(b[2].contains(0));
    check(b[3].contains(1));
    check(b[4].x1 ==  2);
}




/*  Given a function f and a starting point x1, this routine searches
    along the x-axis in the positive (scale>0) or negative (scale<0)
    direction until it reaches a point x2 where f(x2) has an opposite
    sign from f(x1).

    If such a point is found, a RootBracket containing the two
    points x1 and x2, as well as the function values in those points,
    is returned.  If not, an exception is thrown.

    If the function is exactly zero in one of the endpoints, a
    RootBracket starting and ending at that point is returned.

    On the first iteration, x2 = x1+scale.  Hence, scale should be a
    characteristic scale for the function (i.e. a scale over which the
    function changes significantly).  Thereafter, the interval is expanded
    geometrically by multiplying scale by a constant factor for each iteration.
*/
RootBracket!(T, ReturnType!F) bracketFrom(F, T)
    (scope F f, T x1, T scale, int maxIterations=40)
{
    return bracketFrom(f, x1, f(x1), scale, maxIterations);
}


//  ditto
RootBracket!(T, ReturnType!F) bracketFrom(F, T, R)
    (scope F f, T x1, R fx1, T scale, int maxIterations=40)
    if (is(ReturnType!F : R))
in
{
    assert (scale != 0, "scale must be nonzero");
    assert (maxIterations > 0, "maxIterations must be positive");
}
body
{
    enum expandFactor = 1.6;
    real step = scale;

    foreach (i; 0 .. maxIterations)
    {
        immutable x2 = x1 + step;
        immutable fx2 = f(x2);
        if (fx1 * fx2 <= 0) return typeof(return)(x1, x2, fx1, fx2);
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




enum HandleNaN : bool { yes = true, no = false }

/*  This function does essentially the same as bracketFrom(), except that
    it expands the interval in $(I both) directions until it brackets a root.

    Unlike bracketFrom(), this function can optionally handle NaNs, in
    the way that if the function returns a NaN at either endpoint the routine
    stops searching in that direction.  Instead it tries to locate the point
    where the function starts returning NaN, and uses that as one endpoint
    while continuing the search in the other direction using bracketFrom().
    This feature is disabled by default.  Pass HandleNaN.yes as the last
    argument to enable it.
*/
RootBracket!(T, ReturnType!Func) bracketOut(T, Func)
    (scope Func f, T x1, T scale, int maxIterations = 40,
     HandleNaN handleNaN = HandleNaN.no)
in
{
    assert (scale != 0, "scale must be nonzero");
    assert (maxIterations > 0, "maxIterations must be positive");
}
body
{
    enum expandFactor = 1.6;

    real x2 = x1 + scale;
    real fx1 = f(x1);
    real fx2 = f(x2);

    foreach (i; 0 .. maxIterations)
    {
        // Check whether interval brackets a root
        if (fx1 * fx2 <= 0)  return typeof(return)(x1, x2, fx1, fx2);

        // Check for NaN
        if (handleNaN)
        {
            if (isNaN(fx1))
            {
                auto n = findNaN(f, x2, x1, fx2, scale*1e-6f);
                return bracketFrom(f, n.xValid, x2-n.xValid, maxIterations-i);
            }
            else if (isNaN(fx2))
            {
                auto n = findNaN(f, x1, x2, fx1, scale*1e-6f);
                return bracketFrom(f, n.xValid, x1-n.xValid, maxIterations-i);
            }
        }

        // Expand interval in the direction where f(x) is closest to zero.
        if (fabs(fx1) < fabs(fx2))
        {
            x1 += expandFactor * (x1 - x2);
            fx1 = f(x1);
        }
        else
        {
            x2 += expandFactor * (x2 - x1);
            fx2 = f(x2);
        }
    }

    enforceNE(false, NE.Limit);
    assert(0);
}


unittest
{
    real f(real x) { return x^^3; }

    auto bracket = bracketOut(&f, 1.0L, 0.1L);
    check(f(bracket.x1) == bracket.y1);
    check(f(bracket.x2) == bracket.y2);
    check(bracket.y1 * bracket.y2 <= 0);
    check(bracket.contains(0));

    bracket = bracketOut(&f, -1.0L, -0.1L);
    check(f(bracket.x1) == bracket.y1);
    check(f(bracket.x2) == bracket.y2);
    check(bracket.y1 * bracket.y2 <= 0);
    check(bracket.contains(0));
}


unittest
{
    // Check the handleNaN feature
    real f(real x) { return x <= 0 ? real.nan : sin(x); }

    auto bracket = bracketOut(&f, 0.001L, 0.001L, 40, HandleNaN.yes);
    check(bracket.x1 > 0);
    check(f(bracket.x1) == bracket.y1);
    check(f(bracket.x2) == bracket.y2);
    check(bracket.y1 * bracket.y2 <= 0);
    check(bracket.contains(PI/2));
    check(isNaN(f(bracket.x1-0.000001)));
}




/** Use bisection to find the point where the given predicate goes from
    returning false to returning true.

    Params:
        f               =   The function.
        predicate       =   The predicate, which must take a point and
                            the function value at that point and return
                            a boolean.
        xTrue           =   A point where the predicate is true.
        xFalse          =   A point where the predicate is false.
        fTrue           =   (optional) The value of f at xTrue.
        fFalse          =   (optional) The value of f at xFalse.
        xTolerance      =   Success: When the absolute distance between
                            xTrue and xFalse is less than this number,
                            the function returns.
        maxIterations   =   Failure: When the algorithm has failed to
                            produce a result after maxIterations bisections,
                            an exception is thrown.

    Returns:
    A tuple containing values named xTrue, xFalse, fTrue, and fFalse, which
    satisfy
    ---
    f(xTrue) == fTrue
    f(xFalse) == fFalse
    predicate(xTrue, fTrue) == true
    predicate(xFalse, fFalse) == false
    abs(xTrue-xFalse) <= xTolerance
    ---

    Example:
    ---
    // Find a root by bisection
    auto r = bisect(
        (real x) { return x^^3; },
        (real x, real fx) { return fx < 0; },
        -1.0L, 1.5L, 1e-10L
        );

    // Let's check if we got the right answer.
    enum root = 0.0L;
    assert (abs(r.xTrue - root) <= 1e-10);
    assert (abs(r.xFalse - root) <= 1e-10);

    assert (r.fTrue < 0);
    assert (r.fFalse >= 0);
    assert (abs(r.xTrue - r.xFalse) <= 1e-10);
    assert (r.xNaN < 0);
    assert (abs(r.xValid - r.xNan) <= 1e-6);
    ---
*/
Tuple!(T, "xTrue", T, "xFalse", R, "fTrue", R, "fFalse")
bisect(F, T, R = ReturnType!F)
    (scope F f, bool delegate(T, R) predicate, T xTrue, T xFalse,
     T xTolerance, int maxIterations=40)
{
    return bisect(f, predicate, xTrue, xFalse, f(xTrue), f(xFalse),
        xTolerance, maxIterations);
}


/// ditto
Tuple!(T, "xTrue", T, "xFalse", R, "fTrue", R, "fFalse")
bisect(F, T, R = ReturnType!F)
    (scope F f, bool delegate(T, R) predicate, T xTrue, T xFalse,
     R fTrue, R fFalse, T xTolerance, int maxIterations=40)
    if (isFloatingPoint!T && isFloatingPoint!R)
in
{
    assert (predicate(xTrue, fTrue) == true, "Predicate is false at xTrue");
    assert (predicate(xFalse, fFalse) == false, "Predicate is true at xFalse");
    assert (xTolerance > 0, "xTolerance must be positive");
}
body
{
    foreach (i; 0 .. maxIterations)
    {
        if (fabs(xTrue-xFalse) <= xTolerance)
            return typeof(return)(xTrue, xFalse, fTrue, fFalse);

        immutable xMid = (xTrue + xFalse) / 2;
        immutable fMid = f(xMid);
        if (predicate(xMid, fMid))
        {
            xTrue = xMid;
            fTrue = fMid;
        }
        else
        {
            xFalse = xMid;
            fFalse = fMid;
        }
    }

    enforceNE(false, NE.Limit);
    assert(0);
}

 
unittest
{
    auto r = bisect(
        (real x) { return x^^3; },
        (real x, real fx) { return fx < 0; },
        -1.0L, 1.5L, 1e-10L
        );

    enum root = 0.0L;
    check(abs(r.xTrue - root) <= 1e-10);
    check(abs(r.xFalse - root) <= 1e-10);

    check(r.fTrue == r.xTrue^^3);
    check(r.fFalse == r.xFalse^^3);
    check(r.fTrue < 0);
    check(r.fFalse >= 0);
    check(abs(r.xTrue - r.xFalse) <= 1e-10);
}





/*  Use bisection to find the point where a function starts returning NaN.

    This function has been superseded (and is now implemented in terms of)
    the more general bisect().  It may be removed in the future.

    Params:
        f               =   The function.
        xValid          =   A point where the function is known to
                            return a valid result.
        xNaN            =   A point where the function is known to
                            return NaN.
        fValid          =   (optional) The value of f at xValid.
        xTolerance      =   Success: When the absolute distance between
                            xValid and xNaN is less than this number,
                            the function returns.
        maxIterations   =   Failure: When the algorithm has failed to
                            produce a result after maxIterations bisections,
                            an exception is thrown.

    Returns:
    A tuple containing the values of xValid, xNaN, and fValid, where
    xValid and xNaN now satisfy the success criterion
    ---
    abs(xValid-xNaN) <= xTolerance
    ---

    Example:
    ---
    auto r = findNaN(&log, 1.0, -1.0, 1e-6);
    assert (r.xValid >= 0);
    assert (r.xNaN < 0);
    assert (abs(r.xValid - r.xNan) <= 1e-6);
    ---
*/
Tuple!(T, "xValid", T, "xNaN", ReturnType!Func, "fValid") findNaN(Func, T)
    (scope Func f, T xValid, T xNaN, T xTolerance, int maxIterations=40)
in
{
    assert (isNaN(f(xNaN)), "f(xNaN) is not NaN");
}
body
{
    auto fValid = f(xValid);
    return findNaN(f, xValid, xNaN, fValid, xTolerance, maxIterations);
}


//  ditto
Tuple!(T, "xValid", T, "xNaN", R, "fValid") findNaN(Func, T, R)
    (scope Func f, T xValid, T xNaN, R fValid, T xTolerance, int maxIterations=40)
    if (isUnaryFunction!(Func, R, T) && isFloatingPoint!R)
in
{
    assert (!isNaN(fValid), "fValid is NaN");
    assert (isNaN(f(xNaN)), "f(xNaN) is not NaN");
    assert (xTolerance > 0, "xTolerance must be positive");
}
body
{
    auto b = bisect(
        f,
        (T x, R fx) { return isNaN(fx); },
        xNaN, xValid, R.nan, fValid,
        xTolerance, maxIterations);
    return typeof(return)(b.xFalse, b.xTrue, b.fFalse);
}

    
unittest
{
    auto r = findNaN(&log, 1.0, -1.0, 1e-6);
    check (!isNaN(r.fValid));
    check (log(r.xValid) == r.fValid);
    check (r.xValid >= 0);
    check (r.xNaN < 0);
    check (abs(r.xValid-r.xNaN) <= 1e-6);
}
