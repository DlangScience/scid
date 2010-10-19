/** Functions for numerical integration (quadrature) and differentiation.

    It is recommended to read the following before using any of the
    functions in this module, as it discusses issues that are common
    to several of them.
    
    Integration:

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
    use the general-purpose integrate() function.  It tries to select
    the method which is most likely to succeed for the given interval.

    All the integration functions have a set of input parameters in common:
    The function to integrate, the limits of the integral, and some accuracy
    requirements.  All of these are described in the documentation for
    integrate().

    Most of the quadrature routines in SciD are ported from from
    $(LINK2 http://www.netlib.org/quadpack,QUADPACK)
    which is a well-known FORTRAN package for numerical evaluation of
    integrals using
    $(LINK2 http://en.wikipedia.org/wiki/Gaussian_quadrature,Gaussian quadrature).
    Gaussian quadrature is
    $(LINK2 http://mathworld.wolfram.com/NumericalIntegration.html,considered)
    to be the most accurate method available for numerical integration of
    smooth functions.  The QUADPACK integrators are named integrateQ*() below.

    In addition, SciD has support for
    $(LINK2 http://mathworld.wolfram.com/DoubleExponentialIntegration.html,double exponential integration).
    The DE algorithms are ported from Takuya Ooura's
    $(LINK2 http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html,DE-Quadrature package),
    and are available in this module through the integrateDE*() functions.

    Note that SciD contains complete ports of QUADPACK and DE-Quadrature,
    but not all the methods are available through this module yet.  This
    will happen, but until then they can be accessed directly from the
    scid.ports.quadpack and scid.ports.intde packages, albeit with rather
    ugly C/FORTRAN-style interfaces.

    Differentiation:

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
    f(x) is already known it requires one more function evaluation compared
    to the forward-/backward-difference formulas. This is of particular
    importance when one needs to evaluate several derivatives, as is the
    case with e.g. gradients and Jacobians.

    Some functions in this module, like diff(), use another method by Ridders,
    which extrapolates the results of finite-difference formulas like the
    above to make the approximation more accurate – usually a lot more so.
    This requires several function evaluations and is therefore also a lot
    slower than simple finite differences.

    All the differentiation methods in this module take an optional
    scale parameter. This is a "characteristic scale" for the function,
    i.e. a scale over which the function changes significantly. It is
    worth experimenting a bit with this parameter, as it can have a
    drastic impact on accuracy. (For an extreme example, see the "Example"
    section for the diff() function below.) For practical reasons, 1.0 is
    chosen as the default scale in this module. Another, probably more
    common, choice that sometimes works well is to set scale=x, where x
    is the point at which the derivative is taken.


    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009-2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.calculus;


import std.algorithm: max;
import std.math;
import std.traits;

import scid.core.memory;
import scid.core.traits;
import scid.internal.calculus.integrate_qng;
import scid.ports.intde.intde1;
import scid.ports.quadpack.qage;
import scid.ports.quadpack.qagie;
import scid.ports.quadpack.qagse;
import scid.ports.quadpack.qawce;
import scid.exception;
import scid.matrix;
import scid.types;

version(unittest) import scid.core.testing;




/** Integrate a function from a to b.

    This module contains a lot of different integration routines, and
    integrate() aims to be a simple interface to several of
    these. Currently, it only calls the integrateQAGS() and
    integrateQAGI() routines, but this will change in the future.
    In particular, it should take into account whether a function
    is oscillatory, if it has discontinuities and singularities in
    the integration interval, and so on.

    Params:
        f       = The function to integrate.  This can be given as a
                    function pointer, a delegate, or a struct/class
                    which defines opCall().
        a       = The lower limit of integration.
        b       = The upper limit of integration.
        epsRel  = (optional) The requested relative accuracy.
        epsAbs  = (optional) The requested absolute accuracy.

    Note:
    This function is different from the other integrateXXX() routines
    in this module in that it allows a and b to take on the special
    values Real.infinity and -Real.infinity, where Real is the
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
Result!Real integrate (Func, Real)(Func f, Real a, Real b,
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
    check (approxEqual(logsqrtResult.value, -4.0L, 0.0L, logsqrtResult.error));
    check (logsqrtResult.error < abs(logsqrtResult.value*1e-10L));
        

    // Tests for when both limits are infinite.
    real gauss(real x)
    {
        return exp(-x*x);
    }

    auto gaussResult1 =
        integrate(&gauss, -real.infinity, real.infinity, 1e-10L);
    check (approxEqual(gaussResult1.value, sqrt(PI), 0.0L, gaussResult1.error));
    check (gaussResult1.error < abs(gaussResult1.value*1e-10));

    auto gaussResult2 =
        integrate(&gauss, real.infinity, -real.infinity, 1e-10L);
    check (approxEqual(gaussResult2.value,-sqrt(PI), 0.0L, gaussResult2.error));

    auto gaussResult3 =
        integrate(&gauss, -real.infinity, -real.infinity, 1e-10L);
    check (gaussResult3.value == 0.0);

    auto gaussResult4 =
        integrate(&gauss, real.infinity, real.infinity, 1e-10L);
    check (gaussResult4.value == 0.0);


    // Tests for when one limit is infinite.
    real u(real x)  { return log( x)/(1 + 100*x*x); }
    real l(real x)  { return log(-x)/(1 + 100*x*x); }
    
    real ulExact = -PI * log(10.0)/20.0;
    auto uResult = integrate(&u, 0.0L, real.infinity, 0.001L);
    auto muResult = integrate(&u, real.infinity, 0.0L, 0.001L);
    auto lResult = integrate(&l, -real.infinity, 0.0L, 0.001L);
    auto mlResult = integrate(&l, 0.0L, -real.infinity, 0.001L);

    check (approxEqual(uResult.value, ulExact, 0.0L, uResult.error));
    check (approxEqual(lResult.value, ulExact, 0.0L, lResult.error));
    check (uResult.value == -muResult.value);
    check (lResult.value == -mlResult.value);
}




/** Calculate the integral of f over the finite interval (a,b) using a simple
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
Result!Real integrateQNG(Func, Real)(Func f, Real a, Real b,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
{
    return scid.internal.calculus.integrate_qng.qng(f, a, b, epsRel, epsAbs);
}


unittest
{
    double f(double x) { return x^^2 * log(1/x); }
    auto result = integrateQNG(&f, 0.0, 1.0);
    check (abs(result.value - 1.0/9) <= result.error);
}




// Integration rules for integrateQAG()
enum GaussKronrod { rule15 = 1, rule21, rule31, rule41, rule51, rule61 }




/** Calculate the integral of f over the finite interval (a,b) using
    a simple globally adaptive integrator.

    This method is suitable for functions without singularities
    or discontinuities, which are too difficult for integrateQNG(),
    and, in particular for functions with oscillating behaviour of a
    non-specific type.

    It is possible to choose between 6 pairs of
    $(LINK2 http://en.wikipedia.org/wiki/Gauss-Kronrod_quadrature_formula,
    Gauss-Kronrod quadrature formulae)
    for the rule evaluation component.  The pairs of high
    degree of precision are suitable for handling integration
    difficulties due to a strongly oscillating integrand, whereas
    the lower-order rules are more efficient for well-behaved integrands.
    The rule parameter may take on the following values, corresponding to
    15-, 21-, 31-, 41-, 51-, and 61-point Gauss-Kronrod rules:
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
    (Func f, Real a, Real b, GaussKronrod rule = GaussKronrod.rule31,
     Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
{
    Real result, abserr;
    int neval, ier, last;

    enum int limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord;

    qage(f, a, b, epsAbs, epsRel, rule, limit, result, abserr,
        neval, ier, alist.ptr, blist.ptr, rlist.ptr, elist.ptr,
        iord.ptr, last);

    switch (ier)
    {
    case 0:     return typeof(return)(result, abserr);

    case 1:     enforceNE(false, NE.Limit);
    case 2:     enforceNE(false, NE.Roundoff);
    case 3:     enforceNE(false, NE.Behaviour);
    case 6:     enforceNE(false, NE.InvalidInput);
    default:    assert(0);
    }
}


unittest
{
    double f(double x) { return cos(100 * sin(x)); }
    auto i = integrateQAG(&f, 0.0, cast(double) PI, GaussKronrod.rule61);
    check (isAccurate(i.value, i.error, 0.062787400491492695655, 1e-6));
}




/** Calculate the integral of f over the finite interval (a,b) using a
    general-purpose integration algorithm.

    integrateQAGS() is an integrator based on globally adaptive interval
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
Result!Real integrateQAGS(Func, Real)(Func f, Real a, Real b,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
{
    Real result, abserr;
    int neval, ier;
    int last;

    enum int limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord;

    qagse(f, a, b, epsAbs, epsRel, limit, result, abserr, neval, ier,
        alist.ptr, blist.ptr, rlist.ptr, elist.ptr, iord.ptr, last);

    switch (ier)
    {
    case 0:     return typeof(return)(result, abserr);

    case 1:     enforceNE(false, NE.Limit);
    case 2:     enforceNE(false, NE.Roundoff);
    case 3:     enforceNE(false, NE.Behaviour);
    case 4:     enforceNE(false, NE.Accuracy);
    case 5:     enforceNE(false, NE.Convergence);
    case 6:     enforceNE(false, NE.InvalidInput);
    default:    assert(0);
    }
}


unittest
{
    double f(double x) { return 1/sqrt(abs(x-1.0/3)); }
    auto i = integrateQAGS(&f, 0.0, 1.0);
    double expected = 2 * (sqrt(2.0/3) + sqrt(1.0/3));
    check (isAccurate(i.value, i.error, expected, 1e-8));
}




// Integration ranges for integrateQAGI()
enum Infinite { upper = 1, lower = -1 }




/** Calculate the integral of f over an infinite interval.

    The infinite range is mapped onto a finite interval and subsequently
    the same strategy as in integrateQAGS() is applied.

    To integrate f over the interval (-infinity,infinity), use the
    first form.  To integrate f over the interval (-infinity,a) or
    (a,infinity) use the second form with inf=Infinite.lower or
    inf=Infinite.upper, respectively.

    Example:
    ---
    // Slowly convergent integral over infinite interval,
    // integrand with endpoint singularity.
    double f(double x) { return (1 + 10*x)^^(-2) / sqrt(x); }
    auto i = integrateQAGI(&f, 0.0, Infinite.upper, 1e-8);
    ---
*/
Result!Real integrateQAGI(Func, Real)
    (Func f, Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
{
    return integrateQAGI(f, Real.init, cast(Infinite) 2, epsRel, epsAbs);
}


/// ditto
Result!Real integrateQAGI(Func, Real)(Func f, Real a, Infinite inf,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
{
    Real result, abserr;
    int neval, ier, last;

    enum int limit = 500;
    Real[limit] alist, blist, rlist, elist;
    int[limit] iord;

    qagie(f, a, inf, epsAbs, epsRel, limit, result, abserr, neval, ier,
        alist.ptr, blist.ptr, rlist.ptr, elist.ptr, iord.ptr, last);

    NE errCode;
    switch (ier)
    {
        case 0: return Result!(Real)(result, abserr);
        case 1: errCode = NE.Limit; break;
        case 2: errCode = NE.Roundoff; break;
        case 3: errCode = NE.Behaviour; break;
        case 4: errCode = NE.Accuracy; break;
        case 5: errCode = NE.Convergence; break;
        case 6: errCode = NE.InvalidInput; break;
        default: assert(0);
    }
    enforceNE(false, errCode);
    assert(0);
}


unittest
{
    double f(double x) { return (1 + 10*x)^^(-2) / sqrt(x); }
    auto i = integrateQAGI(&f, 0.0, Infinite.upper);
    double expected = PI / (2 * sin(PI/2)) / sqrt(10.0);
    check(isAccurate(i.value, i.error, expected, 1e-6));
}




/** Calculate a Cauchy principal value integral.

    Use this to calculate the integral of f(x)/(x-c) over the finite
    interval (a,b), where f(x) is smooth on the entire interval and
    c is not one of the endpoints.

    The strategy is globally adaptive.  Modified Clenshaw-Curtis
    integration is used on those intervals containing the point
    x = c.

    Example:
    ---
    // Integrate cos(x-1)/(x-1) over the interval (0,3)
    real f(real x) { return cos(x-1); }
    auto i = integrateQAWC(&f, 0.0L, 3.0L, 1.0L, 1e-15L);
    ---
*/
Result!Real integrateQAWC(Func, Real)(Func f, Real a, Real b, Real c,
    Real epsRel = cast(Real) 1e-6, Real epsAbs = cast(Real) 0)
in
{
    assert (epsAbs >= 0 && epsRel >= 0, "Requested accuracy is negative.");
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

    NE errCode;
    switch (ier)
    {
        case 0: return typeof(return)(result, abserr);
        case 1: errCode = NE.Limit; break;
        case 2: errCode = NE.Roundoff; break;
        case 3: errCode = NE.Behaviour; break;
        case 6: errCode = NE.InvalidInput; break;
        default: assert(0);
    }
    enforceNE(false, errCode);
    assert(0);
}


unittest
{
    // This is just a simple test of basic functionality.
    // More thorough unittests are in the scid.ports.quadpack.qawce module.
    enum eps = 1e-15L;
    enum expect = 0.085576905873896861036L;
    static real f(real x) { return cos(x-1); }
    auto i = integrateQAWC(&f, 0.0L, 3.0L, 1.0L, eps);
    check(abs(i.value-expect) < eps*expect && i.error <= eps);
}




/** Calculate the integral of f over the finite interval (a,b) using
    double exponential integration.

    Example:
    ---
    double f(double x) { return x^^2 * log(1/x); }
    auto i = integrateDE(&f, 0.0, 1.0);
    ---
*/
Result!Real integrateDE(Func, Real)(Func f, Real a, Real b,
    Real epsRel = cast(Real) 1e-6)
{
    Real result, error;
    intde(f, a, b, epsRel, &result, &error);
    enforceNE(error >= 0, NE.Convergence);
    return typeof(return)(result, error);
}


unittest
{
    double f(double x) { return x^^2 * log(1/x); }
    auto result = integrateDE(&f, 0.0, 1.0, 1e-7);
    check (isAccurate(result.value, result.error, 1.0/9, 1e-6));
}




/** Calculate the derivative of a function.

    This function uses Ridders' method of extrapolating the results
    of finite difference formulas for consecutively smaller step sizes,
    with an improved stopping criterion described in the Numerical Recipes
    books by Press et al.

    This method gives a much higher degree of accuracy in the answer
    compared to a single finite difference calculation, but requires
    more function evaluations; typically 6 to 12. The maximum number
    of function evaluations is 2*tableauSize.

    Params:
        f = The function of which to take the derivative.
        x = The point at which to take the derivative.
        scale = A "characteristic scale" over which the function
            changes. (optional)
        tableauSize = The values of the consecutive approximations
            are stored in a tableauSize-by-tableauSize triangular
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
    (Func f, real x, real scale = 1.0, size_t tableauSize = 10)
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
    check (approxEqual(r.value, 140.73773557129660339L,
        0.0L, min(r.error, sqrt(real.epsilon))));
}




/** Calculate the Jacobian matrix associated with a set of m functions
    in n variables using a central-difference approximation to the
    Jacobian.
    
    This method requires 2n function evaluations, but is more
    accurate than the faster jacobian2() function. The relative
    accuracy in the result is, at best, on the order of
    (real.epsilon)^(2/3).
    
    Params:
        f = The set of functions. This is typically a function or delegate
            which takes a vector as input and returns a vector. If the
            function takes a second vector parameter, this is assumed to
            be a buffer for the return value, and will be used as such.

        x = The point at which to take the derivative.

        scale = A "characteristic scale" over which the function changes.
            (optional)

        m = The number of functions in f. If this is negative, as it
            is by default, the function will be called once just to
            determine the length of the returned vector. Providing this
            number is therefore a simple way of speeding up this routine.
            (optional)

        buffer = A buffer of length at least m*n, for storing the calculated
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

    auto j = jacobian(&f, p);

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

    auto jFaster = jacobian(&g, p);
    ---

    See_also:
    $(UL $(LI Wikipedia:
        $(LINK2 http://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant,Jacobian matrix and determinant).
    ))
*/
MatrixView!Real jacobian (Real, Func) (Func f, Real[] x, real scale=1.0,
    int m=-1, Real[] buffer=null)
{
    mixin (newFrame);

    static assert (isFloatingPoint!Real,
        "jacobian: Not a floating-point type: "~Real.stringof);
    static assert (isVectorField!(Func, Real),
        "jacobian: Invalid function type ("~Func.stringof~"), or "
       ~"function type doesn't match parameter type ("~Real.stringof~")");

    // Do we have to evaluate the function once just to determine how
    // long the f-vector is? (That would be stupid.)
    if (m < 0)  m = f(x).length;

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
    auto j1 = jacobian(&f, x);
    auto j2 = jacobian(&f, x, 1.0, 2, buffer);
    
    double e = 1e-6;
    check (approxEqual(j1[0,0], 2.0, e) && approxEqual(j1[0,1],  1.0, e)
      &&  approxEqual(j1[1,0], 1.0, e) && approxEqual(j1[1,1], -1.0, e));

    check (approxEqual(j2[0,0], 2.0, e) && approxEqual(j2[0,1],  1.0, e)
      &&  approxEqual(j2[1,0], 1.0, e) && approxEqual(j2[1,1], -1.0, e));
}




/** Calculate the Jacobian using forward-/backward-difference
    methods, also known as 2-point formulas.
    
    This is less accurate than the central-difference method
    used in the jacobian() function, but requires only half
    as many function evaluations. The relative accuracy is,
    at best, on the order of sqrt(real.epsilon).

    This function is used more or less like jacobian(). The
    differences lie in the parameters described below, which
    are all optional.

    Params:
        scale = abs(scale) is the characteristic scale of the
            function, and the sign of scale determines which
            differentiation method is used. If scale is
            negative, the backward difference method is used, and
            if scale is positive, forward differences are used.
            (optional)

        fx = The result of evaluating the function at the point x.
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
MatrixView!Real jacobian2 (Real, Func) (Func f, Real[] x, real scale=1.0,
    Real[] fx=null, Real[] buffer=null)
{
    mixin (newFrame);

    static assert (isFloatingPoint!Real,
        "jacobian2: Not a floating-point type: "~Real.stringof);
    static assert (isVectorField!(Func, Real),
        "jacobian2: Invalid function type ("~Func.stringof~"), or "
       ~"function type doesn't match parameter type ("~Real.stringof~")");

    // Calculate f(x) if it's not provided.
    if (fx.length == 0)  fx = f(x);

    // Determine dimensions.
    immutable size_t m = fx.length;
    immutable size_t n = x.length;
    
    buffer.length = m*n;
    auto jaco = MatrixView!Real(buffer, m, n);

    // Allocate workspace.
    auto fxph = newStack!Real(m);

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
    auto j = jacobian(&f, p);
    check (approxEqual(j.array, answer, e, e));

    j = jacobian(&f, p, 1.0, 4, buffer);
    check (approxEqual(j.array, answer, e, e));


    // The next ones are less accurate.
    e = 1e-6;

    // Forward difference
    j = jacobian2(&f, p);
    check (approxEqual(j.array, answer, e, e));

    // Backward difference
    j = jacobian2(&f, p, -1.0, f(p), buffer);
    check (approxEqual(j.array, answer, e, e));

} 




/** Calculate the gradient of a function of several variables.

    This function calculates a central-difference approximation to
    the gradient of a function f. The error in the result is, at best,
    on the order of sqrt(real.epsilon). The function f is evaluated 2n
    times, where n is the length of the vector x.

    Params:
        f = The function of which to find the gradient.
        x = The point at which to find the gradient.
        scale = A "characteristic scale" over which the function
            changes significantly. (optional)
        buffer = A buffer of the same length as x, for the returned
            gradient vector. (optional)

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
    (Func f, Real[] x, real scale=1.0, Real[] buffer=null)
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
    check (approxEqual(g, gtest, 1e-10, 1e-10));
}



/** Calculate the gradient of a function of several variables
    using Ridders' method.

    This function uses diff() to calculate the derivative
    in each direction. It is therefore more accurate, but slower, than
    the gradient() function. The function f is typically evaluated
    between 6n and 12n times, where n is the length of x.
    See the documentation for diff() for more information on this method.

    This function is used in the same way as gradient().
*/
Real[] gradientR (Real, Func)
    (Func f, Real[] x, real scale=1.0, Real[] buffer=null)
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

    check (approxEqual(g, ans, sqrt(real.epsilon)));
}




/** Calculate the Hessian matrix of a function of several variables
    using a central-difference approximation.

    This function stores its results in an n-by-n symmetric matrix,
    where n is the number of variables (i.e. the length of x).
    The function f is evaluated 1+2n^2 times.

    Params:
        f = The function of which to calculate the Hessian.
        x = The point at which to calculate the Hessian.
        scale = A "characteristic scale" over which the function
            changes significantly. (optional)
        buffer = A buffer of size at least n(n+1)/2, for storing
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
    (Func f, Real[] x, real scale=1.0, Real[] buffer=null)
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

    check (approxEqual(h[0,0], fxx(x), 1e-6L));
    check (approxEqual(h[0,1], fxy(x), 1e-6L));
    check (approxEqual(h[1,0], fyx(x), 1e-6L));
    check (approxEqual(h[1,1], fyy(x), 1e-6L));
}
