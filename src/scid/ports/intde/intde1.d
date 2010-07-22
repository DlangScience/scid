/** D port of intde1.c from Takuya OOURA's DE-Quadrature package
    ($(LINK http://www.kurims.kyoto-u.ac.jp/~ooura/)).

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.intde.intde1;


import std.math;




/** I = integral of f(x) over (a,b)

    Params:
        f         = integrand f(x)
        a         = lower limit of integration
        b         = upper limit of integration
        eps       = relative error requested
        i         = approximation to the integral
        err       = estimate of the absolute error

    Remarks:
    <pre>
        function
            f(x) needs to be analytic over (a,b).
        relative error
            eps is relative error requested excluding 
            cancellation of significant digits.
            i.e. eps means : (absolute error) / 
                             (integral_a^b |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination (m >= mmax).
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has 
                              discontinuous points or sharp 
                              peaks over (a,b).
                              you must divide the interval 
                              (a,b) at this points.
                           2. relative error of f(x) is 
                              greater than eps.
                           3. f(x) has oscillatory factor 
                              and frequency of the oscillation 
                              is very high.
        </pre>
*/
void intde(Real, Func)(Func f, Real a, Real b, Real eps, Real *i, Real *err)
{
    /* ---- adjustable parameter ---- */
    enum MMAX = 256;
    enum : Real
    {
        EFS = 0.1,
        HOFF = 8.5
    }
    /* ------------------------------ */
    int m;
    Real epsln, epsh, h0, ehp, ehm, epst, ba, ir, h, iback, 
        irback, t, ep, em, xw, xa, wg, fa, fb, errt, errh, errd;
    
    epsln = 1 - log(EFS * eps);
    epsh = sqrt(EFS * eps);
    h0 = HOFF / epsln;
    ehp = exp(h0);
    ehm = 1 / ehp;
    epst = exp(-ehm * epsln);
    ba = b - a;
    ir = f((a + b) * 0.5) * (ba * 0.25);
    *i = ir * PI;
    *err = fabs(*i) * epst;
    h = 2 * h0;
    m = 1;
    do {
        iback = *i;
        irback = ir;
        t = h * 0.5;
        do {
            em = exp(t);
            ep = PI_2 * em;
            em = PI_2 / em;
            do {
                xw = 1 / (1 + exp(ep - em));
                xa = ba * xw;
                wg = xa * (1 - xw);
                fa = f(a + xa) * wg;
                fb = f(b - xa) * wg;
                ir += fa + fb;
                *i += (fa + fb) * (ep + em);
                errt = (fabs(fa) + fabs(fb)) * (ep + em);
                if (m == 1) *err += errt * epst;
                ep *= ehp;
                em *= ehm;
            } while (errt > *err || xw > epsh);
            t += h;
        } while (t < h0);
        if (m == 1) {
            errh = (*err / epst) * epsh * h0;
            errd = 1 + 2 * errh;
        } else {
            errd = h * (fabs(*i - 2 * iback) + 4 * fabs(ir - 2 * irback));
        }
        h *= 0.5;
        m *= 2;
    } while (errd > errh && m < MMAX);
    *i *= h;
    if (errd > errh) {
        *err = -errd * m;
    } else {
        *err = errh * epsh * m / (2 * EFS);
    }
}


unittest
{
    real f1(real x) { return 1.0/sqrt(x); }
    real f2(real x) { return sqrt(4 - x*x); }

    real i, err;
    real eps = 1e-15;

    intde(&f1, 0.0L, 1.0L, eps, &i, &err);
    assert (approxEqual(i, 2.0, eps));

    intde(&f2, 0.0L, 2.0L, eps, &i, &err);
    assert (approxEqual(i, PI, eps));
}

    


/** I = integral of f(x) over (a,infinity), 
    f(x) has not oscillatory factor.

    Params:
        f         = integrand f(x)
        a         = lower limit of integration
        eps       = relative error requested
        i         = approximation to the integral
        err       = estimate of the absolute error
    
    Remarks:
    <pre>
        function
            f(x) needs to be analytic over (a,infinity).
        relative error
            eps is relative error requested excluding 
            cancellation of significant digits.
            i.e. eps means : (absolute error) / 
                             (integral_a^infinity |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination (m >= mmax).
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has 
                              discontinuous points or sharp 
                              peaks over (a,infinity).
                              you must divide the interval 
                              (a,infinity) at this points.
                           2. relative error of f(x) is 
                              greater than eps.
                           3. f(x) has oscillatory factor 
                              and decay of f(x) is very slow 
                              as x -> infinity.
    </pre>
*/
void intdei(Real, Func)(Func f, Real a, Real eps, Real *i, Real *err)
{
    /* ---- adjustable parameter ---- */
    enum  MMAX = 256;
    enum : Real
    {
        EFS = 0.1,
        HOFF = 11.0
    }
    /* ------------------------------ */
    int m;
    Real epsln, epsh, h0, ehp, ehm, epst, ir, h, iback, irback, 
        t, ep, em, xp, xm, fp, fm, errt, errh, errd;
    
    epsln = 1 - log(EFS * eps);
    epsh = sqrt(EFS * eps);
    h0 = HOFF / epsln;
    ehp = exp(h0);
    ehm = 1 / ehp;
    epst = exp(-ehm * epsln);
    ir = f(a + 1);
    *i = ir * PI_2;
    *err = fabs(*i) * epst;
    h = 2 * h0;
    m = 1;
    do {
        iback = *i;
        irback = ir;
        t = h * 0.5;
        do {
            em = exp(t);
            ep = PI_4 * em;
            em = PI_4 / em;
            do {
                xp = exp(ep - em);
                xm = 1 / xp;
                fp = f(a + xp) * xp;
                fm = f(a + xm) * xm;
                ir += fp + fm;
                *i += (fp + fm) * (ep + em);
                errt = (fabs(fp) + fabs(fm)) * (ep + em);
                if (m == 1) *err += errt * epst;
                ep *= ehp;
                em *= ehm;
            } while (errt > *err || xm > epsh);
            t += h;
        } while (t < h0);
        if (m == 1) {
            errh = (*err / epst) * epsh * h0;
            errd = 1 + 2 * errh;
        } else {
            errd = h * (fabs(*i - 2 * iback) + 4 * fabs(ir - 2 * irback));
        }
        h *= 0.5;
        m *= 2;
    } while (errd > errh && m < MMAX);
    *i *= h;
    if (errd > errh) {
        *err = -errd * m;
    } else {
        *err = errh * epsh * m / (2 * EFS);
    }
}


unittest
{
    real f3(real x) { return 1.0/(1 + x*x); }
    real f4(real x) { return exp(-x)/sqrt(x); }

    real i, err;
    real eps = 1e-15;

    intdei(&f3, 0.0L, eps, &i, &err);
    assert (approxEqual(i, PI_2, eps));

    intdei(&f4, 0.0L, eps, &i, &err);
    assert (approxEqual(i, sqrt(PI), eps));
}




/** I = integral of f(x) over (a,infinity), 
    f(x) has oscillatory factor :
    ---
    f(x) = g(x) * sin(omega * x + theta) as x -> infinity
    ---

    Params:
        f         = integrand f(x)
        a         = lower limit of integration
        omega     = frequency of oscillation
        eps       = relative error requested
        i         = approximation to the integral
        err       = estimate of the absolute error

    Remarks:
    <pre>
        function
            f(x) needs to be analytic over (a,infinity).
        relative error
            eps is relative error requested excluding 
            cancellation of significant digits.
            i.e. eps means : (absolute error) / 
                             (integral_a^R |f(x)| dx).
            eps does not mean : (absolute error) / I.
        error message
            err >= 0 : normal termination.
            err < 0  : abnormal termination (m >= mmax).
                       i.e. convergent error is detected :
                           1. f(x) or (d/dx)^n f(x) has 
                              discontinuous points or sharp 
                              peaks over (a,infinity).
                              you must divide the interval 
                              (a,infinity) at this points.
                           2. relative error of f(x) is 
                              greater than eps.
    </pre>
*/
void intdeo(Real, Func)(Func f, Real a, Real omega, Real eps,
    Real *i, Real *err)
{
    /* ---- adjustable parameter ---- */
    enum
    {
        MMAX = 256,
        LMAX = 5
    }

    enum : Real
    {
        EFS = 0.1,
        ENOFF = 0.40,
        PQOFF = 2.9,
        PPOFF = -0.72
    }
    /* ------------------------------ */
    int n, m, l, k;
    Real epsln, epsh, frq4, per2, pp, pq, ehp, ehm, ir, h, iback, 
        irback, t, ep, em, tk, xw, wg, xa, fp, fm, errh, tn, errd;
    
    epsln = 1 - log(EFS * eps);
    epsh = sqrt(EFS * eps);
    n = cast(int) (ENOFF * epsln);
    frq4 = fabs(omega) * M_2_PI;
    per2 = PI / fabs(omega);
    pq = PQOFF / epsln;
    pp = PPOFF - log(pq * pq * frq4);
    ehp = exp(2 * pq);
    ehm = 1 / ehp;
    xw = exp(pp - PI_2);
    *i = f(a + sqrt(xw * (per2 * 0.5)));
    ir = *i * xw;
    *i *= per2 * 0.5;
    *err = fabs(*i);
    h = 2;
    m = 1;
    do {
        iback = *i;
        irback = ir;
        t = h * 0.5;
        do {
            em = exp(2 * pq * t);
            ep = PI_4 * em;
            em = PI_4 / em;
            tk = t;
            do {
                xw = exp(pp - ep - em);
                wg = sqrt(frq4 * xw + tk * tk);
                xa = xw / (tk + wg);
                wg = (pq * xw * (ep - em) + xa) / wg;
                fm = f(a + xa);
                fp = f(a + xa + per2 * tk);
                ir += (fp + fm) * xw;
                fm *= wg;
                fp *= per2 - wg;
                *i += fp + fm;
                if (m == 1) *err += fabs(fp) + fabs(fm);
                ep *= ehp;
                em *= ehm;
                tk += 1;
            } while (ep < epsln);
            if (m == 1) {
                errh = *err * epsh;
                *err *= eps;
            }
            tn = tk;
            while (fabs(fm) > *err) {
                xw = exp(pp - ep - em);
                xa = xw / tk * 0.5;
                wg = xa * (1 / tk + 2 * pq * (ep - em));
                fm = f(a + xa);
                ir += fm * xw;
                fm *= wg;
                *i += fm;
                ep *= ehp;
                em *= ehm;
                tk += 1;
            }
            fm = f(a + per2 * tn);
            em = per2 * fm;
            *i += em;
            if (fabs(fp) > *err || fabs(em) > *err) {
                l = 0;
                for (;;) {
                    l++;
                    tn += n;
                    em = fm;
                    fm = f(a + per2 * tn);
                    xa = fm;
                    ep = fm;
                    em += fm;
                    xw = 1;
                    wg = 1;
                    for (k = 1; k <= n - 1; k++) {
                        xw = xw * (n + 1 - k) / k;
                        wg += xw;
                        fp = f(a + per2 * (tn - k));
                        xa += fp;
                        ep += fp * wg;
                        em += fp * xw;
                    }
                    wg = per2 * n / (wg * n + xw);
                    em = wg * fabs(em);
                    if (em <= *err || l >= LMAX) break;
                    *i += per2 * xa;
                }
                *i += wg * ep;
                if (em > *err) *err = em;
            }
            t += h;
        } while (t < 1);
        if (m == 1) {
            errd = 1 + 2 * errh;
        } else {
            errd = h * (fabs(*i - 2 * iback) + pq * fabs(ir - 2 * irback));
        }
        h *= 0.5;
        m *= 2;
    } while (errd > errh && m < MMAX);
    *i *= h;
    if (errd > errh) {
        *err = -errd;
    } else {
        *err *= m * 0.5;
    }
}


unittest
{
    real f5(real x) { return sin(x)/x; }
    real f6(real x) { return cos(x)/sqrt(x); }

    real i, err;
    real eps = 1e-15;

    intdeo(&f5, 0.0L, 1.0L, eps, &i, &err);
    assert (approxEqual(i, PI_2, eps));

    intdeo(&f6, 0.0L, 1.0L, eps, &i, &err);
    assert (approxEqual(i, sqrt(PI_2), eps));
}
