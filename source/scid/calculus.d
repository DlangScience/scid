/** Functions for numerical integration (quadrature) and differentiation.

    It is recommended to read the following introduction before using
    any of the functions in this module, as it discusses features and
    issues that are common to several of them.

    # Integration

    This module contains several different integration methods, suitable
    for a wide range of problems.  The choice of which method to use
    depends on several things, such as whether the integration
    interval is finite or infinite, whether the integrand oscillates, or
    whether it exhibits other difficult behaviour such as singularities or
    discontinuities in the integration interval.  Choosing the right method
    can both speed up the calculation and make it more accurate.
    Therefore, each method below has a short description of which types of
    integral it is suitable for.

    When you just want to "get things done" without any fuss,
    use the general-purpose $(D integrate()) function.  It tries to select
    the method which is most likely to succeed for the given interval.

    All the integration functions have a set of input parameters in common:
    The function to integrate, the limits of the integral, and some accuracy
    requirements.  All of these are described in the documentation for
    $(D integrate()).

    Most of the quadrature routines in SciD are ported from from
    $(LINK2 http://www.netlib.org/quadpack,QUADPACK)
    which is a well-known FORTRAN package for numerical evaluation of
    integrals using
    $(LINK2 http://en.wikipedia.org/wiki/Gaussian_quadrature,Gaussian quadrature).
    Gaussian quadrature is
    $(LINK2 http://mathworld.wolfram.com/NumericalIntegration.html,considered)
    to be the most accurate method available for numerical integration of
    smooth functions.  The QUADPACK integrators are named $(D integrateQ*()) below.

    In addition, SciD has support for
    $(LINK2 http://mathworld.wolfram.com/DoubleExponentialIntegration.html,double exponential integration).
    The DE algorithms are ported from Takuya Ooura's
    $(LINK2 http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html,DE-Quadrature package),
    and are available in this module through the $(D integrateDE*()) functions.

    # Differentiation

    The simplest, but fastest, form of numerical differentiation is using
    finite differences. The functions in this module use three different
    finite-difference formulas: Forward differences,
    ---
            f(x+h) - f(x)
    f'(x) = -------------
                  h
    ---
    backward differences,
    ---
            f(x) - f(x-h)
    f'(x) = -------------
                  h
    ---
    and central differences,
    ---
            f(x+h) - f(x-h)
    f'(x) = ---------------
                  2 h
    ---
    Of these three, the latter is the most accurate, but in the cases where
    $(D f(x)) is already known it requires one more function evaluation compared
    to the forward-/backward-difference formulas. This is of particular
    importance when one needs to evaluate several derivatives, as is the
    case with e.g. gradients and Jacobians.

    Some functions in this module, like $(D diff()), use another method by Ridders,
    which extrapolates the results of finite-difference formulas like the
    above to make the approximation more accurate – usually a lot more so.
    This requires several function evaluations and is therefore also a lot
    slower than simple finite differences.

    All the differentiation methods in this module take an optional
    $(D scale) parameter. This is a "characteristic scale" for the function,
    i.e. a scale over which the function changes significantly. It is
    worth experimenting a bit with this parameter, as it can have a
    drastic impact on accuracy. (For an extreme example, see the "Example"
    section for the $(D diff()) function below.) For practical reasons, $(D 1.0) is
    chosen as the default scale in this module. Another, probably more
    common, choice that sometimes works well is to set $(D scale = x), where
    $(D x) is the point at which the derivative is taken.


    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009-2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
    Macros:
        D = <b><code>$0</code></b>
        INFTY = &infin;
        SUP = <sup>$0</sup>
*/
module scid.calculus;


import std.algorithm;
import std.conv;
import std.math;
import std.traits;

import scid.core.fortran;
import scid.core.memory;
import scid.core.traits;
import scid.internal.calculus.integrate_qng;
import scid.ports.intde.intde1;
import scid.ports.quadpack.qage;
import scid.ports.quadpack.qagie;
import scid.ports.quadpack.qagpe;
import scid.ports.quadpack.qagse;
import scid.ports.quadpack.qawce;
import scid.ports.quadpack.qawfe;
import scid.ports.quadpack.qawoe;
import scid.ports.quadpack.qawse;
import scid.matrix;
import scid.types;
import scid.util;

version(unittest)
{
    import std.range;
    import scid.core.testing;
}




/** Integrate a function from $(D a) to $(D b).

    This module contains a lot of different integration routines, and
    $(D _integrate()) aims to be a simple interface to some of
    these. Currently, it only calls the $(D integrateQAGS()) and
    $(D integrateQAGI()) routines, but this may change in the future.
    In particular, it should take into account whether a function
    is oscillatory, if it has discontinuities and singularities in
    the integration interval, and so on.

    Params:
        f       = The function to _integrate.  This can be given as a
                    function pointer, a delegate, or a struct/class
                    which defines $(D opCall()).
        a       = The lower limit of integration.
        b       = The upper limit of integration.
        epsRel  = (optional) The requested relative accuracy.
        epsAbs  = (optional) The requested absolute accuracy.

    Note:
    This function is different from the other $(D integrateXXX()) routines
    in this module in that it allows $(D a) and $(D b) to take on the special
    values $(D Real.infinity) and $(D -Real.infinity), where $(D Real) is the
    floating-point type being used.

    Example:
    ---
    // Let's integrate the function
    //
    //              log(x)
    //     f(x) = ----------
    //                     2
    //            1 + 100 x
    //
    // from x=0 to x=infinity.

    real f(real x)  { return log(x)/(1 + 100*x*x); }
    auto result = integrate(&f, 0.0L, real.infinity, 0.001L);

    real exact = -PI * log(10.0)/20.0;
    assert (abs(result.value - exact) < result.error);
    ---
*/
Result!Real integrate (Func, Real)(scope Func f, Real a, Real b,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
{
    if (isFinite(a))
    {
        if (isFinite(b))
            return integrateQAGS(f, a, b, epsRel, epsAbs);

        else if (b == Real.infinity)
            return integrateQAGI(f, a, Infinite.upper, epsRel, epsAbs);

        else if (b == -Real.infinity)
            return -integrateQAGI(f, a, Infinite.lower, epsRel, epsAbs);
    }
    else if (a == -Real.infinity)
    {
        if (isFinite(b))
            return integrateQAGI(f, b, Infinite.lower, epsRel, epsAbs);

        else if (b == Real.infinity)
            return integrateQAGI(f, epsRel, epsAbs);

        else if (b == -Real.infinity)
            return Result!Real(0.0, 0.0);
    }
    else if (a == Real.infinity)
    {
        if (isFinite(b))
            return -integrateQAGI(f, b, Infinite.upper, epsRel, epsAbs);

        else if (b == -Real.infinity)
            return -integrateQAGI(f, epsRel, epsAbs);

        else if (b == Real.infinity)
            return Result!Real(0.0, 0.0);
    }

    // a or b is NaN.
    return Result!Real(Real.nan, Real.nan);
}


unittest
{
    // Tests for when both limits are finite.
    real logsqrt(real x) { return (x <= 0.0 ? 0.0 : log(x)/sqrt(x)); }
    auto logsqrtResult = integrate(&logsqrt, 0.0L, 1.0L, 0.0L, 1e-10L);
    assert (approxEqual(logsqrtResult.value, -4.0L, 0.0L, logsqrtResult.error));
    assert (logsqrtResult.error < abs(logsqrtResult.value*1e-10L));


    // Tests for when both limits are infinite.
    real gauss(real x)
    {
        return exp(-x*x);
    }

    auto gaussResult1 =
        integrate(&gauss, -real.infinity, real.infinity, 1e-10L);
    assert (approxEqual(gaussResult1.value, sqrt(PI), 0.0L, gaussResult1.error));
    assert (gaussResult1.error < abs(gaussResult1.value*1e-10));

    auto gaussResult2 =
        integrate(&gauss, real.infinity, -real.infinity, 1e-10L);
    assert (approxEqual(gaussResult2.value,-sqrt(PI), 0.0L, gaussResult2.error));

    auto gaussResult3 =
        integrate(&gauss, -real.infinity, -real.infinity, 1e-10L);
    assert (gaussResult3.value == 0.0);

    auto gaussResult4 =
        integrate(&gauss, real.infinity, real.infinity, 1e-10L);
    assert (gaussResult4.value == 0.0);


    // Tests for when one limit is infinite.
    real u(real x)  { return log( x)/(1 + 100*x*x); }
    real l(real x)  { return log(-x)/(1 + 100*x*x); }

    real ulExact = -PI * log(10.0)/20.0;
    auto uResult = integrate(&u, 0.0L, real.infinity, 0.001L);
    auto muResult = integrate(&u, real.infinity, 0.0L, 0.001L);
    auto lResult = integrate(&l, -real.infinity, 0.0L, 0.001L);
    auto mlResult = integrate(&l, 0.0L, -real.infinity, 0.001L);

    assert (approxEqual(uResult.value, ulExact, 0.0L, uResult.error));
    assert (approxEqual(lResult.value, ulExact, 0.0L, lResult.error));
    assert (uResult.value == -muResult.value);
    assert (lResult.value == -mlResult.value);
}




/** Calculate the integral of $(D f(x)) over the finite interval ($(D a),$(D b))
    using a simple
    non-adaptive automatic integrator, based on a sequence of rules with
    increasing degree of algebraic precision.

    This method should only be used for well-behaved integrands, or when
    speed is a lot more important than accuracy.

    Example:
    ---
    // Despite the statement above, integrateQNG() can handle some
    // difficulties, such as the endpoint singularity in the following
    // example:
    double f(double x) { return x^^2 * log(1/x); }
    auto i = integrateQNG(&f, 0.0, 1.0);
    ---
*/
Result!Real integrateQNG(Func, Real)(scope Func f, Real a, Real b,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
    in
    {
        assert (epsAbs > 0 || epsRel >= 50*Real.epsilon,
            "Requested accuracy is too small.");
    }
    body
{
    return scid.internal.calculus.integrate_qng.qng(f, a, b, epsRel, epsAbs);
}


unittest
{
    double f(double x) { return x^^2 * log(1/x); }
    auto result = integrateQNG(&f, 0.0, 1.0);
    assert (abs(result.value - 1.0/9) <= result.error);
}




// Integration rules for integrateQAG()
enum GaussKronrod { rule15 = 1, rule21, rule31, rule41, rule51, rule61 }




/** Calculate the integral of $(D f(x)) over the finite interval
    ($(D a),$(D b)) using a simple globally adaptive integrator.

    This method is suitable for functions without singularities
    or discontinuities which are too difficult for $(D integrateQNG()),
    and, in particular, for functions with oscillating behaviour of a
    non-specific type.

    It is possible to choose between 6 pairs of
    $(LINK2 http://en.wikipedia.org/wiki/Gauss–Kronrod_quadrature_formula,
    Gauss–Kronrod quadrature formulae)
    for the rule evaluation component.  The pairs of high
    degree of precision are suitable for handling integration
    difficulties due to a strongly oscillating integrand, whereas
    the lower-order rules are more efficient for well-behaved integrands.
    The $(D rule) parameter may take on the following values, corresponding to
    15-, 21-, 31-, 41-, 51-, and 61-point Gauss–Kronrod rules:
    ---
    enum GaussKronrod { rule15, rule21, rule31, rule41, rule51, rule61 }
    ---

    Example:
    ---
    // Strongly oscillating integrand, better use highest-order rule.
    double f(double x) { return cos(100 * sin(x)); }
    auto i = integrateQAG(&f, 0.0, cast(double) PI, GaussKronrod.rule61);
    ---
*/
Result!Real integrateQAG(Func, Real)
    (scope Func f, Real a, Real b, GaussKronrod rule = GaussKronrod.rule31,
     Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
    in
    {
        assert (epsAbs > 0 || epsRel >= 50*Real.epsilon,
            "Requested accuracy is too small.");
    }
    body
{
    Real result, abserr;
    int neval, ier, last;

    enum int limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord;

    qage(f, a, b, epsAbs, epsRel, rule, limit, result, abserr,
        neval, ier, alist.ptr, blist.ptr, rlist.ptr, elist.ptr,
        iord.ptr, last);
    checkQuadpackStatus(ier, result, abserr);

    return typeof(return)(result, abserr);
}


unittest
{
    double f(double x) { return cos(100 * sin(x)); }
    auto i = integrateQAG(&f, 0.0, cast(double) PI, GaussKronrod.rule61);
    assert (isAccurate(i.value, i.error, 0.062787400491492695655, 1e-6));
}




/** Calculate the integral of $(D f(x)) over the finite interval
    ($(D a),$(D b)) using a general-purpose integration algorithm.

    $(D _integrateQAGS()) is an integrator based on globally adaptive interval
    subdivision in connection with extrapolation by the Epsilon
    algorithm, which eliminates the effects of integrand singularities
    of several types.

    This is the most "intelligent" of the general-purpose finite-interval
    integration methods, and the one that best handles bad integrand
    behaviour.  It is fairly expensive in terms of processing time, though,
    so if that is an issue you may want to investigate the integrand
    in more detail and try one or a combination of the more specialised
    integration methods.

    Example:
    ---
    // This function has an internal singularity at x = 1/3.
    double f(double x) { return 1/sqrt(abs(x-1.0/3)); }
    auto i = integrateQAGS(&f, 0.0, 1.0);
    ---
*/
Result!Real integrateQAGS(Func, Real)(scope Func f, Real a, Real b,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
    in
    {
        assert (epsAbs > 0 || epsRel >= 50*Real.epsilon,
            "Requested accuracy is too small.");
    }
    body
{
    Real result, abserr;
    int neval, ier;
    int last;

    enum int limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord;

    qagse(f, a, b, epsAbs, epsRel, limit, result, abserr, neval, ier,
        alist.ptr, blist.ptr, rlist.ptr, elist.ptr, iord.ptr, last);
    checkQuadpackStatus(ier, result, abserr);

    return typeof(return)(result, abserr);
}


unittest
{
    double f(double x) { return 1/sqrt(abs(x-1.0/3)); }
    auto i = integrateQAGS(&f, 0.0, 1.0);
    double expected = 2 * (sqrt(2.0/3) + sqrt(1.0/3));
    assert (isAccurate(i, expected, 1e-8));
}




/** Calculate the integral of $(D f(x)) over the finite interval ($(D a),$(D b)),
    taking into account known points of special difficulty
    inside the interval.

    This routine uses the same integration method as $(D integrateQAGS()),
    but allows you to specify an array of points where the
    integrand has internal singularities, discontinuities or
    other types of bad behaviour.
*/
Result!Real integrateQAGP(Func, Real)(scope Func f, Real a, Real b,
    Real[] trouble, Real epsRel = cast(Real) 1e-6,
    Real epsAbs = cast(Real) 0)
    in
    {
        assert (epsAbs > 0 || epsRel >= 50*Real.epsilon,
            "Requested accuracy is too small.");
    }
    body
{
    if (trouble.length == 0)  return integrateQAGS(f, a, b, epsRel, epsAbs);

    mixin(newFrame);

    Real result, abserr;
    int neval, ier;
    int last;

    enum int limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord, level;

    // The array passed to qagp must have two unused slots.
    // Unfortunately, this will sometimes cause an allocation.
    int npts2 = toInt(trouble.length) + 2;
    trouble.length = npts2;
    auto pts = newStack!Real(npts2);
    auto ndin = newStack!int(npts2);

    qagpe(f, a, b, npts2, trouble.ptr, epsAbs, epsRel, limit,
        result, abserr, neval, ier,
        alist.ptr, blist.ptr, rlist.ptr, elist.ptr, pts.ptr,
        iord.ptr, level.ptr, ndin.ptr, last);
    checkQuadpackStatus(ier, result, abserr);

    return typeof(return)(result, abserr);
}


unittest
{
    real f(real x) { return abs(x - PI/4)^^(-0.4L); }
    auto i = integrateQAGP(&f, 0.0L, 1.0L, [PI/4], 1e-8L);
    auto expect = ((1-PI/4)^^0.6L + (PI/4)^^0.6L)/0.6L;
    assert (isAccurate(i, expect, 1e-8L));
}




// Integration ranges for integrateQAGI()
enum Infinite { upper = 1, lower = -1 }




/** Calculate the integral of $(D f(x)) over an infinite interval.

    The infinite range is mapped onto a finite interval and subsequently
    the same strategy as in $(D integrateQAGS()) is applied.

    To integrate $(D f(x)) over the interval (-$(INFTY),$(INFTY)), use the
    first form.  To integrate $(D f(x)) over the interval (-$(INFTY),$(D a)) or
    ($(D a),$(INFTY)) use the second form with $(D inf=Infinite.lower) or
    $(D inf=Infinite.upper), respectively.

    Example:
    ---
    // Slowly convergent integral over infinite interval,
    // integrand with endpoint singularity.
    double f(double x) { return (1 + 10*x)^^(-2) / sqrt(x); }
    auto i = integrateQAGI(&f, 0.0, Infinite.upper, 1e-8);
    ---
*/
Result!Real integrateQAGI(Func, Real)
    (scope Func f, Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
{
    return integrateQAGI(f, Real.init, cast(Infinite) 2, epsRel, epsAbs);
}


/// ditto
Result!Real integrateQAGI(Func, Real)(scope Func f, Real a, Infinite inf,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
    in
    {
        assert (epsAbs > 0 || epsRel >= 50*Real.epsilon,
            "Requested accuracy is too small.");
    }
    body
{
    Real result, abserr;
    int neval, ier, last;

    enum int limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord;

    qagie(f, a, inf, epsAbs, epsRel, limit, result, abserr, neval, ier,
        alist.ptr, blist.ptr, rlist.ptr, elist.ptr, iord.ptr, last);
    checkQuadpackStatus(ier, result, abserr);

    return typeof(return)(result, abserr);
}


unittest
{
    double f(double x) { return (1 + 10*x)^^(-2) / sqrt(x); }
    auto i = integrateQAGI(&f, 0.0, Infinite.upper);
    double expected = PI / (2 * sin(PI/2)) / sqrt(10.0);
    assert (isAccurate(i.value, i.error, expected, 1e-6));
}




/** Calculate the integral of an oscillatory function over the
    finite interval ($(D a),$(D b)).

    Use this to calculate the integral of $(D f(x)*cos(omega*x))
    or $(D f(x)*sin(omega*x))
    where $(D f(x)) is the (possibly singular) user-specified function
    and $(D omega) is a known constant.  The weight function is specified
    by setting $(D weight) to $(D Oscillation.cos) or $(D Oscillation.sin).

    The rule evaluation component is based on the modified
    $(LINK2 http://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature,
    Clenshaw–Curtis technique).
    An adaptive subdivision scheme is used in connection with
    an extrapolation procedure, which is a modification of that in
    $(D integrateQAGS()) and allows the algorithm to deal with
    singularities in $(D f(x)).

    See_Also:
    $(D integrateQAWF()), for similar integrals over an infinite interval.

    Example:
    ---
    // Integrate exp(20*(x-1))*sin(256*x) over the interval (0,1)
    real f(real x) { return exp(20*(x-1)); }
    auto i = integrateQAWO(&f, 0.0L, 1.0L, 256.0L, Oscillation.sin, 1e-15L);
    ---
*/
Result!Real integrateQAWO(Func, Real)(scope Func f, Real a, Real b, Real omega,
    Oscillation weight, Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
    in
    {
        assert (epsAbs > 0 || epsRel >= 50*Real.epsilon,
            "Requested accuracy is too small.");
    }
    body
{
    Real result, abserr;
    int neval, ier, last, momcom;

    enum limit = 500;
    enum icall = 1;
    enum maxp1 = 21;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord, nnlog;
    Real[maxp1*25] chebmo;

    qawoe!(Real, Func)(f, a, b, omega, weight, epsAbs, epsRel, limit,
        icall, maxp1, result, abserr, neval, ier, last,
        alist.ptr, blist.ptr, rlist.ptr, elist.ptr, iord.ptr, nnlog.ptr,
        momcom, chebmo.ptr);
    checkQuadpackStatus(ier, result, abserr);

    return typeof(return)(result, abserr);
}


enum Oscillation { cos = 1, sin = 2 }


unittest
{
    enum eps = 1e-15L;
    real f(real x) { return exp(20*(x-1)); }
    auto i = integrateQAWO(&f, 0.0L, 1.0L, 256.0L, Oscillation.sin, eps);
    real expect = (20*sin(256) - 256*cos(256) + 256*exp(-20.0L))/(400+65536);
    assert (isAccurate(i, expect, eps));
}




/** Calculates a Fourier transform integral.

    Use this to calculate the integral of $(D f(x)*cos(omega*x))
    or $(D f(x)*sin(omega*x))
    on the interval ($(D a),$(INFTY)).  The weight function is specified
    by setting $(D weight) to $(D Oscillation.cos) or $(D Oscillation.sin).

    The procedure of $(D integrateQAWO()) is applied on successive finite
    intervals, and convergence acceleration by means of the Epsilon
    algorithm is applied to the series of integral approximations.

    Note:
    This function is unique among the QUADPACK integrators in that
    you must specify an absolute accuracy, not a relative one.

    See_Also:
    $(D integrateQAWO()), for similar integrals over a finite interval.

    Example:
    ---
    // Calculate the integral of cos(x)*exp(-x/64)/sqrt(x) over
    // the interval (0, infinity).
    real f(real x) { return x <= 0 ? 0 : exp(-x/64) / sqrt(x); }
    auto i = integrateQAWF(&f, 0.0L, 1.0L, Oscillation.cos, 1e-15L);
    ---
*/
Result!Real integrateQAWF(Func, Real)(scope Func f, Real a, Real omega,
    Oscillation weight, Real epsAbs)
in
{
    assert (epsAbs > 0, "Requested accuracy is too small.");
}
body
{
    Real result, abserr;
    int neval, ier, lst;

    enum limit = 500;
    enum limlst = 50;
    enum maxp1 = 21;
    Real[limit] alist, blist, rlist, elist;
    Real[limlst] rslst, erlst;
    int[limit] iord, nnlog;
    int[limlst] ierlst;
    Real[maxp1*25] chebmo;

    qawfe!(Real, Func)(f, a, omega, weight, epsAbs, limlst, limit,
        maxp1, result, abserr, neval, ier, rslst.ptr, erlst.ptr, ierlst.ptr,
        lst, alist.ptr, blist.ptr, rlist.ptr, elist.ptr, iord.ptr, nnlog.ptr,
        chebmo.ptr);

    // QAWF has slightly different error codes, so we can't
    // use checkQuadpackStatus()
    string msg, extra;
    switch(ier)
    {
        case 0:
            return typeof(return)(result, abserr);
        case 1:
            msg = "The maximum number of cycles allowed has been reached";
            break;
        case 4:
            msg = "Convergence acceleration algorithm failed to reach the "
                ~"desired accuracy";
            break;
        case 6:
            msg = "Invalid input";
            extra = "If you get this error message, it is a bug in SciD.  "
                ~ "Please report.";
            break;
        case 7:
            msg = "Bad integrand behaviour within one or more of the cycles.";
            break;
        default:
            msg = "Invalid QUADPACK error code";
            extra = "If you get this error message, it is a bug in SciD.  "
                ~ "Please report.";
    }
    if (extra.length == 0)
    {
        extra = "There may be local integration difficulties "
            ~"(singularities, discontinuities, etc.).  If the position "
            ~"of such a difficulty can be determined, try to split "
            ~"up the interval at this point and call the appropriate "
            ~"integrators on the subintervals.";
    }
    throw new IntegrationException(msg, extra, result, abserr);
}


unittest
{
    enum eps = 1e-14L;
    real f(real x) { return x <= 0 ? 0 : exp(-x/64) / sqrt(x); }
    auto i = integrateQAWF(&f, 0.0L, 1.0L, Oscillation.cos, eps);
    real expect = sqrt(PI) * (1 + 1.0L/4096)^^(-0.25L) * cos(atan(64.0L)/2);
    assert (isAccurate(i, expect, 0.0L, eps));
}




// Weight functions for integrateQAWS()
enum Weight { unity = 1, logxa = 2, logbx = 3, logab = 4 }


/** Calculate an integral over the finite interval ($(D a),$(D b)),
    where the integrand has algebraic and/or logarithmic endpoint
    singularities of a known type.

    The integrand is taken to be on the form
    ---
                alpha        beta
    f(x) (x - a)      (b - x)     w(x)
    ---
    where $(D f(x)) is the given function and $(D w(x)) is
    specified by setting the $(D weight) parameter to one of the following:
    ---
    Weight.unity:  w(x) = 1
    Weight.logxa:  w(x) = log(x-a)
    Weight.logbx:  w(x) = log(b-x)
    Weight.logab:  w(x) = log(x-a) log(b-x)
    ---

    A globally adaptive subdivision strategy is applied,
    with modified
    $(LINK2 http://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature,
    Clenshaw–Curtis integration) on those subintervals
    which contain $(D a) or $(D b).

    Example:
    ---
    // Calculate the integral of 1/(sqrt(1-x^^2) * (x+1.5)).
    // Another way to write this integrand is
    //     (x-(-1))^^(-0.5) * (1-x)^^0.5 / (x+1.5),
    // so we set alpha = beta = -0.5.
    real f(real x) { return 1/(x + 1.5L); }
    auto i = integrateQAWS(&f, -1.0L, 1.0L, -0.5L, -0.5L, Weight.unity, 1e-15L);
    ---
*/
Result!Real integrateQAWS(Func, Real)(scope Func f, Real a, Real b,
    Real alpha, Real beta, Weight weight,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
    in
    {
        assert (epsAbs > 0 || epsRel >= 50*Real.epsilon,
            "Requested accuracy is too small.");
    }
    body
{
    Real result, abserr;
    int neval, ier, last;

    enum int limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord;

    int sign = 1;
    if (b <= a)
    {
        if (b == a) return typeof(return)(0, 0);
        swap(a, b);
        sign = -1;
    }

    qawse(f, a, b, alpha, beta, weight, epsAbs, epsRel, limit,
        result, abserr, neval, ier,
        alist.ptr, blist.ptr, rlist.ptr, elist.ptr, iord.ptr, last);
    checkQuadpackStatus(ier, result, abserr);

    return typeof(return)(sign*result, abserr);
}


unittest
{
    enum eps = 1e-15L;
    real f(real x) { return 1/(x + 1.5L); }
    auto i = integrateQAWS(&f, 1.0L, -1.0L, -0.5L, -0.5L, Weight.unity, eps);
    auto expected = -PI / sqrt(1.5L^^2 - 1);
    assert (isAccurate(i, expected, eps));
}




/** Calculate a Cauchy principal value integral.

    Use this to calculate the integral of $(D f(x)/(x-c)) over the finite
    interval ($(D a),$(D b)), where $(D f(x)) is smooth on the entire
    interval and $(D c) is not one of the endpoints.

    The strategy is globally adaptive.  Modified
    $(LINK2 http://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature,
    Clenshaw–Curtis integration) is used on those intervals containing the point
    $(D x = c).

    Example:
    ---
    // Integrate cos(x-1)/(x-1) over the interval (0,3)
    real f(real x) { return cos(x-1); }
    auto i = integrateQAWC(&f, 0.0L, 3.0L, 1.0L, 1e-15L);
    ---
*/
Result!Real integrateQAWC(Func, Real)(scope Func f, Real a, Real b, Real c,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
    in
    {
        assert (epsAbs > 0 || epsRel >= 50*Real.epsilon,
            "Requested accuracy is too small.");
        assert (c != a && c != b, "Singularity at interval endpoint.");
    }
    body
{
    Real result, abserr;
    int neval, ier, last;

    enum limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord;

    qawce!(Real, Func)(f, a, b, c, epsAbs, epsRel, limit, result, abserr,
        neval, ier, alist.ptr, blist.ptr, rlist.ptr, elist.ptr, iord.ptr,
        last);
    checkQuadpackStatus(ier, result, abserr);

    return typeof(return)(result, abserr);
}


unittest
{
    // This is just a simple test of basic functionality.
    // More thorough unittests are in the scid.ports.quadpack.qawc module.
    enum eps = 1e-15L;
    enum expect = 0.085576905873896861036L;
    static real f(real x) { return cos(x-1); }
    auto i = integrateQAWC(&f, 0.0L, 3.0L, 1.0L, eps);
    assert (abs(i.value-expect) < eps*expect && i.error <= eps);
}




/** Calculate the integral of $(D f(x)) over the finite interval
    ($(D a),$(D b)) using double exponential integration.

    Example:
    ---
    double f(double x) { return x^^2 * log(1/x); }
    auto i = integrateDE(&f, 0.0, 1.0);
    ---
*/
Result!Real integrateDE(Func, Real)(scope Func f, Real a, Real b,
    Real epsRel = cast(Real) 1e-6)
{
    Real result, error;
    intde(f, a, b, epsRel, &result, &error);
    if (error < 0) throw new IntegrationException(
        "Integration failed",
        "The integrand may have discontinuities, singularities or "
        ~"oscillatory behaviour in the integration interval.",
        result,
        real.nan);
    return typeof(return)(result, error);
}


unittest
{
    double f(double x) { return x^^2 * log(1/x); }
    auto result = integrateDE(&f, 0.0, 1.0, 1e-7);
    assert (isAccurate(result.value, result.error, 1.0/9, 1e-6));
}




/** Calculate the integral of $(D f(x)) over the infinite interval
    ($(D a),$(INFTY)) using double exponential integration.

    Example:
    ---
    double f(double x) { return (1 + 10*x)^^(-2) / sqrt(x); }
    auto i = integrateDEI(&f, 0.0);
    ---
*/
Result!Real integrateDEI(Func, Real)(scope Func f, Real a,
    Real epsRel = cast(Real) 1e-6)
{
    Real result, error;
    intdei!Real(f, a, epsRel, &result, &error);
    if (error < 0) throw new IntegrationException(
        "Integration failed",
        "The integrand may have discontinuities, singularities or "
        ~"oscillatory behaviour in the integration interval.",
        result,
        real.nan);
    return typeof(return)(result, error);
}


unittest
{
    double f(double x) { return (1 + 10*x)^^(-2) / sqrt(x); }
    auto result = integrateDEI(&f, 0.0);
    double expected = PI / (2 * sin(PI/2)) / sqrt(10.0);
    assert (matchDigits(result.value, expected, 6));
}




/** Calculate the integral of an oscillating function $(D f(x))
    over the infinite interval ($(D a),$(INFTY)) using double exponential
    integration.

    $(D f(x)) is assumed to take the form
    ---
    f(x) = g(x) * sin(omega * x + theta)
    ---
    as $(D x) goes to infinity ($(D theta) is not specified).

    Example:
    ---
    real f(real x) { return x <= 0 ? 0 : cos(x) * exp(-x/64) / sqrt(x); }
    auto i = integrateDEO(&f, 0.0L, 1.0L, 1e-18L);
    ---
*/
Result!Real integrateDEO(Func, Real)(scope Func f, Real a,
    Real omega, Real epsRel = cast(Real) 1e-6)
{
    Real result, error;
    intdeo(f, a, omega, epsRel, &result, &error);
    if (error < 0) throw new IntegrationException(
        "Integration failed",
        "The integrand may have discontinuities, singularities or "
        ~"oscillatory behaviour in the integration interval.",
        result,
        real.nan);
    return typeof(return)(result, error);
}


unittest
{
    real f(real x) { return x <= 0 ? 0 : cos(x) * exp(-x/64) / sqrt(x); }
    auto i = integrateDEO(&f, 0.0L, 1.0L, 1e-18L);
    real expect = sqrt(PI) * (1 + 1.0L/4096)^^(-0.25L) * cos(atan(64.0L)/2);
    assert (matchDigits(i.value, expect, 18));
}




/** Exception thrown when an integration routine fails. */
class IntegrationException : Exception
{
    /** The current estimates of the integration result and
        absolute error.

        Depending on the type of error, these may or may not be
        close to the correct values.  They are, however, very often
        useful for figuring out what has gone wrong.
    */
    immutable real resultEstimate;
    immutable real errorEstimate;       /// ditto


    /** This value will sometimes contain more detailed information
        about the error.
    */
    immutable string extra;


    this(string msg, real resultEstimate, real errorEstimate,
        string file = __FILE__, size_t line = __LINE__)
    {
        this(msg, null, resultEstimate, errorEstimate, file, line);
    }


    this(string msg, string extra, real resultEstimate, real errorEstimate,
        string file = __FILE__, size_t line = __LINE__)
    {
        super (msg, file, line);

        this.extra = extra;
        this.resultEstimate = resultEstimate;
        this.errorEstimate = errorEstimate;
    }


    override string toString()
    {
        return super.toString()
            ~ (extra.length > 0 ? "\n"~extra : "")
            ~ text("\nCurrent estimate: ", resultEstimate, " ± ", errorEstimate);
    }
}


/*  Verify the value of 'ier' returned by QUADPACK routines, and
    throw an appropriate exception if it is not zero.
*/
void checkQuadpackStatus(Real)(int ier, Real result, Real abserr,
    string file = __FILE__, size_t line = __LINE__)
{
    string msg, extra;
    switch(ier)
    {
        case 0: return;
        case 1:
            msg = "Reached the maximum number of iterations or subdivisions";
            extra = "If possible, try using a more special-purpose integrator.";
            break;
        case 2:
            msg = "Roundoff error detected";
            extra = "The occurrence of roundoff error is detected, which "
                ~ "prevents the requested tolerance from being achieved.  "
                ~ "The error may be under-estimated.";
            break;
        case 3:
            msg = "Bad integrand behaviour";
            extra = "Extremely bad integrand behaviour occurs at at some "
                ~ "points of the integration interval.";
            break;
        case 4:
            msg = "Roundoff error detected.";
            extra = "The algorithm does not converge.  Roundoff error is "
                ~ "detected in the extrapolation table.  It is assumed that "
                ~ "the requested tolerance cannot be achieved, and that the "
                ~ "returned result is the best that can be obtained.";
            break;
        case 5:
            msg = "Integral is divergent or converges too slowly";
            break;
        case 6:
            msg = "Invalid input";
            extra = "If you get this error message, it is a bug in SciD.  "
                ~ "Please report.";
            break;
        default:
            msg = "Invalid QUADPACK error code";
            extra = "If you get this error message, it is a bug in SciD.  "
                ~ "Please report.";
    }
    throw new IntegrationException(msg, extra, result, abserr, file, line);
}




/** Calculate the derivative of a function.

    This function uses Ridders' method of extrapolating the results
    of finite difference formulas for consecutively smaller step sizes,
    with an improved stopping criterion described in the Numerical Recipes
    books by Press et al.

    This method gives a much higher degree of accuracy in the answer
    compared to a single finite difference calculation, but requires
    more function evaluations; typically 6 to 12. The maximum number
    of function evaluations is $(D 2*tableauSize).

    Params:
        f = The function of which to take the derivative.
        x = The point at which to take the derivative.
        scale = A "characteristic scale" over which the function
            changes. (optional)
        tableauSize = The values of the consecutive approximations
            are stored in a $(D tableauSize)-by-$(D tableauSize) triangular
            matrix. Sometimes this number can be reduced to limit
            the possible number of function evaluations (see above)
            and thus greater speed, and sometimes it can be increased
            to allow for higher accuracy. (optional)

    Example:
    This example is from Ridders' paper:
    ---
    // Let's take the derivative of
    //                 x
    //                e
    //     f(x) = ----------
    //                     2
    //            sin x - x
    // at x=1.
    real f(real x) { return exp(x)/(sin(x)-x*x); }
    real dfdxAt1 = 140.73773557129660339;

    // scale=0.01 is appropriate for this function.
    auto r = diff(&f, 1.0, 0.01);
    assert (abs(r.value - dfdxAt1) <= r.error);

    // Note that if we use the default scale, it won't work:
    auto s = diff(&f, 1.0);
    writeln(s.error);           // prints "inf"
    ---

    References:
    $(UL
        $(LI
            C. J. F. Ridders,
            $(I Accurate computation of F'(x) and F'(x)F''(x)).
            Advances in Engineering Software, vol. 4 (1982), issue 2, p. 75.)
        $(LI
            W. H. Press, S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
            $(I Numerical Recipes in C++) (2nd ed.).
            Cambridge University Press, 2003.)
    )
*/
Result!real diff (Func)
    (scope Func f, real x, real scale = 1.0, size_t tableauSize = 10)
in
{
    assert (tableauSize > 0);
}
body
{
    mixin(newFrame);

    static assert (__traits(compiles, { real y; Func g; real z = g(y); }),
        "diff: Invalid function type: "~Func.stringof);

    // Set up Romberg tableau.
    auto workspace = newStack!real((tableauSize*(tableauSize+1))/2);
    auto tab =
        MatrixView!(real, Storage.Triangular)(workspace, tableauSize);

    // Let's keep the notation in order.
    real h = scale == 0.0 ? 1.0 : abs(scale);

    // Divide h by the factor FAC for each finite-difference
    // approximation. Ridders uses FAC=2.0 in his paper, but
    // the NR book uses FAC=1.4.
    // TODO: Experiment and see which value is best.
    enum real FAC = 2.0;
    enum real FACSQ = FAC*FAC;
    enum real INVFAC = 1.0/FAC;

    // From the NR book: Stop when the difference between consecutive
    // approximations is bigger than SAFE*error, where error is an
    // estimate of the absolute error in the current (best) approximation.
    enum int SAFE = 2;

    // First approximation: A_0
    tab[0,0] = (f(x+h) - f(x-h))/(2.0*h);

    real result, error = real.infinity;
    for (size_t n=1; n<tableauSize; n++)
    {
        // Decrease h.
        h *= INVFAC;

        // Compute A_n
        tab[0,n] = (f(x+h) - f(x-h))/(2.0*h);

        real facm = 1.0;
        for (size_t m=1; m<=n; m++)
        {
            facm *= FACSQ;

            // Compute B_(n-1), C_(n-2), ...
            real upLeft  = tab[m-1, n-1];
            real up      = tab[m-1, n];
            real current = (facm*up - upLeft)/(facm-1);
            tab[m,n] = current;

            // Calculate and check error.
            real currentError = max(abs(current - upLeft), abs(current - up));
            if (currentError <= error)
            {
                result = current;
                error = currentError;
            }
        }

        if (abs(tab[n,n]-tab[n-1,n-1]) >= SAFE*error)  break;
    }

    return Result!real(result, error);
}

unittest
{
    alias diff!(float function(float)) diff_float;
    alias diff!(double function(double)) diff_double;
    alias diff!(real function(real)) diff_real;

    real f(real x) { return exp(x)/(sin(x)-x*x); }
    auto r = diff(&f, 1.0, 0.01, 5);
    assert (approxEqual(r.value, 140.73773557129660339L,
        0.0L, min(r.error, sqrt(real.epsilon))));
}




/** Calculate the derivative of a function using a two-point
    formula, i.e. a forward- or backward-difference formula.

    This method only evaluates the function once (in addition
    to the function value provided by the user), and is therefore
    the fastest way to compute a derivative.  However, it is also
    the least accurate method, and should only be used if the
    function is very expensive to compute $(I and) you have
    already calculated $(D f(x)).

    Params:
        f = The function to differentiate.
        x = The point at which to take the derivative.
        fx = The function value at $(D x), i.e. $(D f(x))
        scale = A characteristic scale for $(D f).  When this
            is positive, the forward-difference formula
            is used, and when it is negative, the
            backward-difference formula is used.

    Returns:
        An approximation to the derivative of $(D f) at the point
        $(D x).  The relative error in the result is $(I at best)
        on the order of $(D sqrt(real.epsilon)).
        Usually it is much higher.
*/
real diff2(Func)(scope Func f, real x, real fx, real scale = 1.0)
in { assert (scale != 0); }
body
{
    enum sqrtEpsilon = sqrt(real.epsilon);
    immutable xph = x + sqrtEpsilon*scale;
    immutable h = xph - x;
    return (f(xph) - fx) / (xph - x);
}

unittest
{
    real f(real x) { return sin(x) * log(x); }
    real df(real x) { return cos(x) * log(x) + sin(x) / x; }

    foreach (x; iota(0.1, 1.0, 0.1))
    {
        assert (matchDigits(diff2(&f, x, f(x)), df(x), 8));
    }
}




/** Calculate the derivative of a function using a three-point
    formula, a.k.a. a central difference formula.

    The function is evaluated twice, at points just below and
    just above $(D x).

    Params:
        f = The function to differentiate.
        x = The point at which to take the derivative.
        scale = A characteristic scale for $(D f).

    Returns:
        An approximation to the derivative of $(D f) at the point
        $(D x).  The relative error in the result is at best
        on the order of $(D (real.epsilon)^(2/3)), roughly three
        orders of magnitude more accurate than $(D diff2()).
*/
real diff3(Func)(scope Func f, real x, real scale = 1.0)
in { assert (scale != 0); }
body
{
    immutable h = cbrt(real.epsilon) * scale;
    immutable xmh = x - h;
    immutable xph = x + h;
    return (f(xph) - f(xmh))/(xph - xmh);
}

unittest
{
    real f(real x) { return sin(x) * log(x); }
    real df(real x) { return cos(x) * log(x) + sin(x) / x; }

    foreach (x; iota(0.1, 1.0, 0.1))
    {
        assert (matchDigits(diff3(&f, x), df(x), 11));
    }
}




/** Calculate the Jacobian matrix associated with a set of $(I m) functions
    in $(I n) variables using a central-difference approximation to the
    Jacobian.

    This method requires 2$(I n) function evaluations, but is more
    accurate than the faster $(D jacobian2()) function. The relative
    accuracy in the result is, at best, on the order of
    $(D (real.epsilon)^(2/3)).

    Params:
        f = The set of functions. This is typically a function or delegate
            which takes a vector as input and returns a vector. If the
            function takes a second vector parameter, this is assumed to
            be a buffer for the return value, and will be used as such.

        m = The number of functions in $(D f).

        x = The point at which to take the derivative.

        scale = A "characteristic scale" over which the function changes.
            (optional)

        buffer = A buffer of length at least $(I m)*$(I n), for storing the calculated
            Jacobian matrix. (optional)

    Examples:
    ---
    // Let's say we want to find the Jacobian of the function
    //      f(x,y) = (xy, x-y)
    // at the point p = (1, 2).
    real[] p = [ 1.0, 2.0 ];

    // This is the simplest way to do it:
    real[] f(real[] a)
    {
        auto r = new real[2];
        r[0] = a[0] * a[1];
        r[1] = a[0] - a[1];
        return r;
    }

    auto j = jacobian(&f, 2, p);

    // However, if we need to perform this operation many times,
    // we may want to speed things up a bit. To avoid unnecessary
    // repeated allocations, we can add a buffer parameter to the
    // function:
    real[] g(real[] a, real[] r)
    {
        r[0] = a[0] * a[1];
        r[1] = a[0] - a[1];
        return r[0 .. 2];
    }

    auto jFaster = jacobian(&g, 2, p);
    ---

    See_also:
    $(UL $(LI Wikipedia:
        $(LINK2 http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant,Jacobian matrix and determinant).
    ))
*/
MatrixView!Real jacobian (Real, Func) (scope Func f, int m, Real[] x,
    real scale=1.0, Real[] buffer=null)
{
    mixin (newFrame);

    static assert (isFloatingPoint!Real,
        "jacobian: Not a floating-point type: "~Real.stringof);
    static assert (isVectorField!(Func, Real),
        "jacobian: Invalid function type ("~Func.stringof~"), or "
       ~"function type doesn't match parameter type ("~Real.stringof~")");
    assert (m > 0);

    immutable size_t n = x.length;
    buffer.length = m*n;
    auto jaco = MatrixView!Real(buffer, m, n);

    // Allocate workspace.
    auto workspace = newStack!Real(2*m);
    auto fxph = workspace[0 .. m];
    auto fxmh = workspace[m .. $];

    // Determine step length.
    scale = scale == 0.0 ? 1.0 : abs(scale);
    immutable real CBRT_EPSILON = cbrt(real.epsilon); // TODO: make enum
    real h = CBRT_EPSILON * scale;

    foreach (j; 0 .. n)
    {
        // We don't want to change x.
        real save = x[j];

        // Ensure 2h is exactly (x+h)-(x-h)
        real xph = save + h;
        real xmh = save - h;
        real twh = xph - xmh;

        // Evaluate function at both points, provide buffer if possible.
        x[j] = xph;
        static if (isBufferVectorField!(Func,Real))  fxph = f(x, fxph);
        else  fxph = f(x);

        x[j] = xmh;
        static if (isBufferVectorField!(Func,Real))  fxmh = f(x, fxmh);
        else  fxmh = f(x);

        // Restore x.
        x[j] = save;

        // Calculate one column of the Jacobian.
        foreach (i; 0 .. m)
            jaco[i, j] = (fxph[i] - fxmh[i])/twh;
    }

    return jaco;
}

unittest
{
    alias jacobian!(float, float[] function(float[])) jacobian_float;
    alias jacobian!(double, double[] function(double[])) jacobian_double;
    alias jacobian!(real, real[] function(real[])) jacobian_real;

    // Note: This is the example from the doc comment.
    // There are more unittests for jacobian() below jacobian2().
    double[] f(double[] a, double[] r=null)
    {
        if (r.length < 2)  r = new double[2];
        auto x = a[0], y = a[1];

        r[0] = x*y;
        r[1] = x-y;
        return r[0 .. 2];
    }

    double[] x = [ 1.0, 2.0 ];
    double[4] buffer;
    auto j1 = jacobian(&f, 2, x);
    auto j2 = jacobian(&f, 2, x, 1.0, buffer);

    double e = 1e-6;
    assert (approxEqual(j1[0,0], 2.0, e) && approxEqual(j1[0,1],  1.0, e)
      &&  approxEqual(j1[1,0], 1.0, e) && approxEqual(j1[1,1], -1.0, e));

    assert (approxEqual(j2[0,0], 2.0, e) && approxEqual(j2[0,1],  1.0, e)
      &&  approxEqual(j2[1,0], 1.0, e) && approxEqual(j2[1,1], -1.0, e));
}




/** Calculate the Jacobian using forward-/backward-difference
    methods, also known as 2-point formulas.

    This is less accurate than the central-difference method
    used in the $(D jacobian()) function, but requires only half
    as many function evaluations. The relative accuracy is,
    at best, on the order of $(D sqrt(real.epsilon)).

    This function is used more or less like $(D jacobian()). The
    differences lie in the parameters described below, which
    are all optional.

    Params:
        scale = $(D abs(scale)) is the characteristic scale of the
            function, and the sign of scale determines which
            differentiation method is used. If $(D scale) is
            negative, the backward difference method is used, and
            if it is positive, forward differences are used.
            (optional)

        fx = The result of evaluating the function at the point $(D x).
            Providing this saves one function evaluation if you
            for some reason have already calculated the value. (optional)

    Examples:
    ---
    // Continuing with the example for the jacobian() function,
    // we make it even faster by using a 2-point formula.
    auto gp = g(p);
    auto jFastest = jacobian2(&g, p, 1.0, gp, buffer, workspace);
    ---
*/
MatrixView!Real jacobian2 (Real, Func) (scope Func f, int m, Real[] x,
    real scale=1.0, Real[] fx=null, Real[] buffer=null)
{
    mixin (newFrame);

    static assert (isFloatingPoint!Real,
        "jacobian2: Not a floating-point type: "~Real.stringof);
    static assert (isVectorField!(Func, Real),
        "jacobian2: Invalid function type ("~Func.stringof~"), or "
       ~"function type doesn't match parameter type ("~Real.stringof~")");
    assert (m > 0);

    immutable size_t n = x.length;
    buffer.length = m*n;
    auto jaco = MatrixView!Real(buffer, m, n);

    // Allocate workspace.
    auto fxph = newStack!Real(m);

    // Calculate f(x) if it's not provided.
    if (fx.length == 0)
    {
        static if (isBufferVectorField!(Func, Real))
        {
            fx = newStack!Real(m);
            fx = f(x, fx);
        }
        else fx = f(x);
    }
    else assert (fx.length == m);

    // Determine step length.
    scale = scale == 0.0 ? 1.0 : abs(scale);
    real h = sqrt(real.epsilon) * scale;

    for (int j=0; j<n; j++)
    {
        // We don't want to change x.
        real save = x[j];

        // Ensure h is exactly (x+h)-x
        real xph = save + h;
        real hExact = xph - save;

        // Evaluate function at x+h, provide buffer if possible.
        x[j] = xph;
        static if (isBufferVectorField!(Func,Real))  fxph = f(x, fxph);
        else  fxph = f(x);

        // Restore x.
        x[j] = save;

        // Calculate one column of the Jacobian.
        for (int i=0; i<m; i++)
            jaco[i, j] = (fxph[i] - fx[i])/hExact;
    }

    return jaco;
}

unittest
{
    alias jacobian2!(float, float[] function(float[])) jacobian2_float;
    alias jacobian2!(double, double[] function(double[])) jacobian2_double;
    alias jacobian2!(real, real[] function(real[])) jacobian2_real;

    real[] f(real[] x, real[] buf=null)
    {
        if (buf.length < 4)  buf = new real[4];
        real x1 = x[0], x2 = x[1], x3 = x[2];

        buf[0] = x1;
        buf[1] = 5*x3;
        buf[2] = 4*x2*x2 - 2*x3;
        buf[3] = x3*sin(x1);

        return buf[0 .. 4];
    }

    real[] p = [0.0L, 1, 2];
    real[] answer = [1.0L, 0,      0, p[2]*cos(p[0]),
                        0, 0, 8*p[1],              0,
                        0, 5,     -2,      sin(p[0])];
    real[12] buffer;


    // Central difference
    real e = 1e-10;
    auto j = jacobian(&f, 4, p);
    assert (approxEqual(j.array, answer, e, e));

    j = jacobian(&f, 4, p, 1.0, buffer);
    assert (approxEqual(j.array, answer, e, e));


    // The next ones are less accurate.
    e = 1e-6;

    // Forward difference
    j = jacobian2(&f, 4, p);
    assert (approxEqual(j.array, answer, e, e));

    // Backward difference
    j = jacobian2(&f, 4, p, -1.0, f(p), buffer);
    assert (approxEqual(j.array, answer, e, e));
}




/** Calculate the _gradient of a function of several variables.

    This function calculates a central-difference approximation to
    the _gradient of a function $(D f). The error in the result is, at best,
    on the order of $(D sqrt(real.epsilon)). The function $(D f) is evaluated
    2$(I n) times, where $(I n) is the length of the vector $(D x).

    Params:
        f = The function of which to find the _gradient.
        x = The point at which to find the _gradient.
        scale = A "characteristic scale" over which the function
            changes significantly. (optional)
        buffer = A buffer of the same length as $(D x), for the returned
            _gradient vector. (optional)

    Example:
    ---
    // Let's find the gradient of f(x,y) = x exp(y) at
    // the point p = (2,1).
    real f(real[] x)  { return x[0] * exp(x[1]); }
    real[] p = [2.0L, 1.0L];

    auto g = gradient(&f, p);
    ---

    See_also:
    $(UL $(LI Wikipedia:
        $(LINK2 http://en.wikipedia.org/wiki/Gradient,Gradient).
    ))
*/
Real[] gradient (Real, Func)
    (scope Func f, Real[] x, real scale=1.0, Real[] buffer=null)
{
    static assert (isFloatingPoint!Real,
        "gradient: Not a floating-point type: "~Real.stringof);
    static assert (__traits(compiles, { Real[] u; Func g; Real v = g(u); }),
        "gradient: Invalid function type ("~Func.stringof~"), or "
       ~"function type doesn't match parameter type ("~Real.stringof~")");

    immutable size_t n = x.length;
    buffer.length = n;

    // Determine step length.
    scale = scale == 0.0 ? 1.0 : abs(scale);
    real h = cbrt(real.epsilon) * scale;

    foreach (i; 0 .. n)
    {
        // We don't want to modify the vector x.
        real save = x[i];

        // Ensure the difference between x+h and x-h is exactly 2h.
        real xph = save + h;
        real xmh = save - h;
        real twh = xph - xmh;

        // Calculate central difference and divide by 2h.
        x[i] = xph;
        real sum = f(x);
        x[i] = xmh;
        sum -= f(x);
        buffer[i] = sum/twh;

        // Restore vector.
        x[i] = save;
    }

    return buffer[0 .. n];
}

unittest
{
    alias gradient!(float, float function(float[])) gradient_float;
    alias gradient!(double, double function(double[])) gradient_double;
    alias gradient!(real, real function(real[])) gradient_real;

    real f(real[] x)  { return x[0] * exp(x[1]); }
    real[] p = [2.0L, 1.0L];

    auto g = gradient(&f, p);

    real[] fgrad(real[] x) { return [exp(x[1]), x[0]*exp(x[1])]; }
    auto gtest = fgrad(p);
    assert (approxEqual(g, gtest, 1e-10, 1e-10));
}



/** Calculate the gradient of a function of several variables
    using Ridders' method.

    This function uses $(D diff()) to calculate the derivative
    in each direction. It is therefore more accurate, but slower, than
    the $(D gradient()) function. The function $(D f) is typically evaluated
    between 6$(I n) and 12$(I n) times, where $(I n) is the length of $(D x).
    See the documentation for $(D diff()) for more information on this method.

    This function is used in the same way as $(D gradient()).
*/
Real[] gradientR (Real, Func)
    (scope Func f, Real[] x, real scale=1.0, Real[] buffer=null)
{
    static assert (isFloatingPoint!Real,
        "gradientR: Not a floating-point type: "~Real.stringof);
    static assert (__traits(compiles, { Real[] u; Func g; Real v = g(u); }),
        "gradientR: Invalid function type ("~Func.stringof~"), or "
       ~"function type doesn't match parameter type ("~Real.stringof~")");

    immutable size_t n = x.length;
    buffer.length = n;

    // Calculate the derivative in each direction using Ridders' method.
    foreach (i; 0 .. n)
    {
        real save = x[i];

        // Provide a real->real proxy function.
        real g(real y)
        {
            x[i] = y;
            return f(x);
        }

        // Use Ridders' method to differentiate.
        buffer[i] = diff(&g, save, scale, 10);
        x[i] = save;
    }

    return buffer;
}

unittest
{
    alias gradientR!(float, float function(float[])) gradientR_float;
    alias gradientR!(double, double function(double[])) gradientR_double;
    alias gradientR!(real, real function(real[])) gradientR_real;

    real f(real[] v)
    {
        real x = v[0], y = v[1], z = v[2];
        return exp(x) + log(y) + z*z;
    }

    real[] df(real[] v)
    {
        real x = v[0], y = v[1], z = v[2];
        real[] g = new real[3];
        g[0] = exp(x);
        g[1] = 1/y;
        g[2] = 2*z;
        return g;
    }

    real[] p = [ 0.0L, 2.0, 4.0 ];
    real[] ans = df(p);
    auto g = gradient(&f, p);

    assert (approxEqual(g, ans, sqrt(real.epsilon)));
}




/** Calculate the Hessian matrix of a function of several variables
    using a central-difference approximation.

    This function stores its results in an $(I n)-by-$(I n) symmetric matrix,
    where $(I n) is the number of variables (i.e. the length of $(D x)).
    The function $(D f) is evaluated 1+2$(I n)$(SUP 2) times.

    Params:
        f = The function of which to calculate the Hessian.
        x = The point at which to calculate the Hessian.
        scale = A "characteristic scale" over which the function
            changes significantly. (optional)
        buffer = A buffer of size at least $(I n)($(I n)+1)/2, for storing
            the Hessian matrix. (optional)

    Example:
    ---
    // Above, we found the gradient of f(x,y) at the point p.
    // Finding the Hessian matrix is just as simple:
    auto h = hessian(&f, p);
    ---

    See_also:
    $(UL $(LI Wikipedia:
        $(LINK2 http://en.wikipedia.org/wiki/Hessian_matrix,Hessian matrix).
    ))
*/
MatrixView!(Real, Storage.Symmetric) hessian(Real, Func)
    (scope Func f, Real[] x, real scale=1.0, Real[] buffer=null)
{
    static assert (isFloatingPoint!Real,
        "hessian: Not a floating-point type: "~Real.stringof);
    static assert (__traits(compiles, { Real[] u; Func g; Real v = g(u); }),
        "hessian: Invalid function type ("~Func.stringof~"), or "
       ~"function type doesn't match parameter type ("~Real.stringof~")");

    immutable size_t n = x.length;

    buffer.length = (n*n + n)/2;
    auto hess = MatrixView!(Real, Storage.Symmetric)(buffer, n, n);

    // Determine step length.
    enum real SQRT_EPSILON = sqrt(real.epsilon);
    enum real QDRT_EPSILON = sqrt(SQRT_EPSILON);
    scale = scale == 0.0 ? 1.0 : abs(scale);
    real hOffdiag = QDRT_EPSILON*scale;
    real hDiag = cbrt(real.epsilon*(scale^^4));

    // We need the function at x for the diagonal elements.
    real m2fx = -2.0*f(x);

    // Loop over columns.
    foreach (j; 0 .. n)
    {
        // Loop over rows, calculate offdiagonal elements.
        foreach (i; 0 .. j)
        {
            // Don't want to change x.
            real savex = x[i];
            real savey = x[j];

            // Make hx and hy exactly representable numbers.
            real xph = savex + hOffdiag;
            real xmh = savex - hOffdiag;
            real twhx = xph - xmh;

            real yph = savey + hOffdiag;
            real ymh = savey - hOffdiag;
            real twhy = yph - ymh;

            // Sum function value at all four points.
            x[i] = xph;
            x[j] = yph;
            real sum = f(x);

            x[j] = ymh;
            sum -= f(x);

            x[i] = xmh;
            sum += f(x);

            x[j] = yph;
            sum -= f(x);

            // Restore x and save one element of the Hessian.
            x[i] = savex;
            x[j] = savey;
            hess[i, j] = sum/(twhx*twhy);
        }

        // Calculate the diagonal element.
        {
            // Don't want to change x.
            real save = x[j];

            // Make h an exactly representable number.
            real xph = save + hDiag;
            real h = xph - save;

            // We have already calculated -2f(x).
            real sum = m2fx;

            x[j] = xph;
            sum += f(x);

            x[j] = save - hDiag;
            sum += f(x);

            // Reset x and calculate diagonal element.
            x[j] = save;
            hess[j, j] = sum/(h*h);
        }
    }

    return hess;
}

unittest
{
    alias hessian!(float, float function(float[])) hessian_float;
    alias hessian!(double, double function(double[])) hessian_double;
    alias hessian!(real, real function(real[])) hessian_real;

    real f(real[] x) { return cos(x[0]+x[1]) + sin(x[0]*x[1]); }
    real fxx(real[] x) { return -cos(x[0]+x[1]) - x[1]*x[1]*sin(x[0]*x[1]); }
    real fxy(real[] x) { return -cos(x[0]+x[1]) - x[0]*x[1]*sin(x[0]*x[1])
        + cos(x[0]*x[1]); }
    real fyx(real[] x) { return fxy(x); }
    real fyy(real[] x) { return -cos(x[0]+x[1]) - x[0]*x[0]*sin(x[0]*x[1]); }

    real[] x = [ 10., 100.0 ];
    real[3] buffer;

    auto h = hessian(&f, x, PI*0.25, buffer);

    assert (approxEqual(h[0,0], fxx(x), 1e-6L));
    assert (approxEqual(h[0,1], fxy(x), 1e-6L));
    assert (approxEqual(h[1,0], fyx(x), 1e-6L));
    assert (approxEqual(h[1,1], fyy(x), 1e-6L));
}
