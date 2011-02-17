/** D port of intde2.c from Takuya OOURA's DE-Quadrature package
    ($(LINK http://www.kurims.kyoto-u.ac.jp/~ooura/)).

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.intde.intde2;


import std.math;




/** I = integral of f(x) over (a,b)
    
    Examples:
        ---
        auto aw = new real[8000];
        intdeini(aw.length, tiny, eps, aw.ptr);  // initialization of aw
        ...
        intde(f, a, b, aw.ptr, &i, &err);
        ---

    Params:
        lenaw     = length of aw
        tiny      = minimum value that 1/tiny does not 
                    overflow
        eps       = relative error requested
        aw        = points and weights of the quadrature 
                    formula, aw[0...lenaw-1]
        f         = integrand f(x)
        a         = lower limit of integration
        b         = upper limit of integration
        i         = approximation to the integral
        err       = estimate of the absolute error

    Remarks:
    <pre>
        initial parameters
            lenaw > 1000, 
            IEEE double :
                lenaw = 8000;
                tiny = 1.0e-307;
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
            err < 0  : abnormal termination.
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
void intdeini(Real)(size_t lenaw, Real tiny, Real eps, Real *aw)
{
    assert (lenaw > 1000);
    /* ---- adjustable parameter ---- */
    enum : Real { EFS = 0.1, HOFF = 8.5 }
    /* ------------------------------ */
    int nk, k, j;
    Real tinyln, epsln, h0, ehp, ehm, h, t, ep, em, xw, wg;
    
    tinyln = -log(tiny);
    epsln = 1 - log(EFS * eps);
    h0 = HOFF / epsln;
    ehp = exp(h0);
    ehm = 1 / ehp;
    aw[2] = eps;
    aw[3] = exp(-ehm * epsln);
    aw[4] = sqrt(EFS * eps);
    enum NOFF = 5;
    aw[NOFF] = 0.5;
    aw[NOFF + 1] = h0;
    aw[NOFF + 2] = PI_4 * h0;
    h = 2;
    nk = 0;
    k = NOFF + 3;
    do {
        t = h * 0.5;
        do {
            em = exp(h0 * t);
            ep = PI_2 * em;
            em = PI_2 / em;
            j = k;
            do {
                xw = 1 / (1 + exp(ep - em));
                wg = xw * (1 - xw) * h0;
                aw[j] = xw;
                aw[j + 1] = wg * 4;
                aw[j + 2] = wg * (ep + em);
                ep *= ehp;
                em *= ehm;
                j += 3;
            } while (ep < tinyln && j <= lenaw - 3);
            t += h;
            k += nk;
        } while (t < 1);
        h *= 0.5;
        if (nk == 0) {
            if (j > lenaw - 6) j -= 3;
            nk = j - NOFF;
            k += nk;
            aw[1] = nk;
        }
    } while (2 * k - NOFF - 3 <= lenaw);
    aw[0] = k - 3;
}


/// ditto
void intde(Real, Func)(Func f, Real a, Real b, Real *aw, 
    Real *i, Real *err)
{
    int lenawm, nk, k, j, jtmp, jm, m, klim;
    Real epsh, ba, ir, xa, fa, fb, errt, errh, errd, h, iback, irback;
    
    enum NOFF = 5;
    lenawm = cast(int) (aw[0] + 0.5);
    nk = cast(int) (aw[1] + 0.5);
    epsh = aw[4];
    ba = b - a;
    *i = f((a + b) * aw[NOFF]);
    ir = *i * aw[NOFF + 1];
    *i *= aw[NOFF + 2];
    *err = fabs(*i);
    k = nk + NOFF;
    j = NOFF;
    do {
        j += 3;
        xa = ba * aw[j];
        fa = f(a + xa);
        fb = f(b - xa);
        ir += (fa + fb) * aw[j + 1];
        fa *= aw[j + 2];
        fb *= aw[j + 2];
        *i += fa + fb;
        *err += fabs(fa) + fabs(fb);
    } while (aw[j] > epsh && j < k);
    errt = *err * aw[3];
    errh = *err * epsh;
    errd = 1 + 2 * errh;
    jtmp = j;
    while (fabs(fa) > errt && j < k) {
        j += 3;
        fa = f(a + ba * aw[j]);
        ir += fa * aw[j + 1];
        fa *= aw[j + 2];
        *i += fa;
    }
    jm = j;
    j = jtmp;
    while (fabs(fb) > errt && j < k) {
        j += 3;
        fb = f(b - ba * aw[j]);
        ir += fb * aw[j + 1];
        fb *= aw[j + 2];
        *i += fb;
    }
    if (j < jm) jm = j;
    jm -= NOFF + 3;
    h = 1;
    m = 1;
    klim = k + nk;
    while (errd > errh && klim <= lenawm) {
        iback = *i;
        irback = ir;
        do {
            jtmp = k + jm;
            for (j = k + 3; j <= jtmp; j += 3) {
                xa = ba * aw[j];
                fa = f(a + xa);
                fb = f(b - xa);
                ir += (fa + fb) * aw[j + 1];
                *i += (fa + fb) * aw[j + 2];
            }
            k += nk;
            j = jtmp;
            do {
                j += 3;
                fa = f(a + ba * aw[j]);
                ir += fa * aw[j + 1];
                fa *= aw[j + 2];
                *i += fa;
            } while (fabs(fa) > errt && j < k);
            j = jtmp;
            do {
                j += 3;
                fb = f(b - ba * aw[j]);
                ir += fb * aw[j + 1];
                fb *= aw[j + 2];
                *i += fb;
            } while (fabs(fb) > errt && j < k);
        } while (k < klim);
        errd = h * (fabs(*i - 2 * iback) + fabs(ir - 2 * irback));
        h *= 0.5;
        m *= 2;
        klim = 2 * klim - NOFF;
    }
    *i *= h * ba;
    if (errd > errh) {
        *err = -errd * (m * fabs(ba));
    } else {
        *err = *err * aw[2] * (m * fabs(ba));
    }
}


unittest
{
    real f1(real x) { return 1/sqrt(x); }
    real f2(real x) { return sqrt(4 - x*x); }

    real[] aw = new real[8000];
    real i, err;
    real tiny = 1e-307;
    real eps = 1e-15;

    intdeini(aw.length, tiny, eps, aw.ptr);

    intde(&f1, 0.0L, 1.0L, aw.ptr, &i, &err);
    assert (approxEqual(i, 2.0, eps));

    intde(&f2, 0.0L, 2.0L, aw.ptr, &i, &err);
    assert (approxEqual(i, PI, eps));
}




/** I = integral of f(x) over (a,infinity), 
    f(x) has not oscillatory factor.

    Examples:
        ---
        auto aw = new real[8000];
        intdeiini(aw.length, tiny, eps, aw.ptr);  // initialization of aw
        ...
        intdei(f, a, aw.ptr, &i, &err);
        ---

    Params:
        lenaw     = length of aw
        tiny      = minimum value that 1/tiny does not 
                    overflow
        eps       = relative error requested
        aw        = points and weights of the quadrature 
                    formula, aw[0...lenaw-1]
        f         = integrand f(x)
        a         = lower limit of integration
        i         = approximation to the integral
        err       = estimate of the absolute error

    Remarks:
    <pre>
        initial parameters
            lenaw > 1000, 
            IEEE double :
                lenaw = 8000;
                tiny = 1.0e-307;
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
            err < 0  : abnormal termination.
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
void intdeiini(Real)(size_t lenaw, Real tiny, Real eps, Real *aw)
{
    assert (lenaw > 1000);
    /* ---- adjustable parameter ---- */
    enum : Real { EFS = 0.1, HOFF = 11.0 }
    /* ------------------------------ */
    int nk, k, j;
    Real tinyln, epsln, h0, ehp, ehm, h, t, ep, em, xp, xm, 
        wp, wm;
    
    tinyln = -log(tiny);
    epsln = 1 - log(EFS * eps);
    h0 = HOFF / epsln;
    ehp = exp(h0);
    ehm = 1 / ehp;
    aw[2] = eps;
    aw[3] = exp(-ehm * epsln);
    aw[4] = sqrt(EFS * eps);
    enum NOFF = 5;
    aw[NOFF] = 1;
    aw[NOFF + 1] = 4 * h0;
    aw[NOFF + 2] = PI_2 * h0;
    h = 2;
    nk = 0;
    k = NOFF + 6;
    do {
        t = h * 0.5;
        do {
            em = exp(h0 * t);
            ep = PI_4 * em;
            em = PI_4 / em;
            j = k;
            do {
                xp = exp(ep - em);
                xm = 1 / xp;
                wp = xp * ((ep + em) * h0);
                wm = xm * ((ep + em) * h0);
                aw[j] = xm;
                aw[j + 1] = xp;
                aw[j + 2] = xm * (4 * h0);
                aw[j + 3] = xp * (4 * h0);
                aw[j + 4] = wm;
                aw[j + 5] = wp;
                ep *= ehp;
                em *= ehm;
                j += 6;
            } while (ep < tinyln && j <= lenaw - 6);
            t += h;
            k += nk;
        } while (t < 1);
        h *= 0.5;
        if (nk == 0) {
            if (j > lenaw - 12) j -= 6;
            nk = j - NOFF;
            k += nk;
            aw[1] = nk;
        }
    } while (2 * k - NOFF - 6 <= lenaw);
    aw[0] = k - 6;
}


/// ditto
void intdei(Real, Func)(Func f, Real a, Real *aw, Real *i, 
    Real *err)
{
    int lenawm, nk, k, j, jtmp, jm, m, klim;
    Real epsh, ir, fp, fm, errt, errh, errd, h, iback, irback;
    
    enum NOFF = 5;
    lenawm = cast(int) (aw[0] + 0.5);
    nk = cast(int) (aw[1] + 0.5);
    epsh = aw[4];
    *i = f(a + aw[NOFF]);
    ir = *i * aw[NOFF + 1];
    *i *= aw[NOFF + 2];
    *err = fabs(*i);
    k = nk + NOFF;
    j = NOFF;
    do {
        j += 6;
        fm = f(a + aw[j]);
        fp = f(a + aw[j + 1]);
        ir += fm * aw[j + 2] + fp * aw[j + 3];
        fm *= aw[j + 4];
        fp *= aw[j + 5];
        *i += fm + fp;
        *err += fabs(fm) + fabs(fp);
    } while (aw[j] > epsh && j < k);
    errt = *err * aw[3];
    errh = *err * epsh;
    errd = 1 + 2 * errh;
    jtmp = j;
    while (fabs(fm) > errt && j < k) {
        j += 6;
        fm = f(a + aw[j]);
        ir += fm * aw[j + 2];
        fm *= aw[j + 4];
        *i += fm;
    }
    jm = j;
    j = jtmp;
    while (fabs(fp) > errt && j < k) {
        j += 6;
        fp = f(a + aw[j + 1]);
        ir += fp * aw[j + 3];
        fp *= aw[j + 5];
        *i += fp;
    }
    if (j < jm) jm = j;
    jm -= NOFF + 6;
    h = 1;
    m = 1;
    klim = k + nk;
    while (errd > errh && klim <= lenawm) {
        iback = *i;
        irback = ir;
        do {
            jtmp = k + jm;
            for (j = k + 6; j <= jtmp; j += 6) {
                fm = f(a + aw[j]);
                fp = f(a + aw[j + 1]);
                ir += fm * aw[j + 2] + fp * aw[j + 3];
                *i += fm * aw[j + 4] + fp * aw[j + 5];
            }
            k += nk;
            j = jtmp;
            do {
                j += 6;
                fm = f(a + aw[j]);
                ir += fm * aw[j + 2];
                fm *= aw[j + 4];
                *i += fm;
            } while (fabs(fm) > errt && j < k);
            j = jtmp;
            do {
                j += 6;
                fp = f(a + aw[j + 1]);
                ir += fp * aw[j + 3];
                fp *= aw[j + 5];
                *i += fp;
            } while (fabs(fp) > errt && j < k);
        } while (k < klim);
        errd = h * (fabs(*i - 2 * iback) + fabs(ir - 2 * irback));
        h *= 0.5;
        m *= 2;
        klim = 2 * klim - NOFF;
    }
    *i *= h;
    if (errd > errh) {
        *err = -errd * m;
    } else {
        *err *= aw[2] * m;
    }
}


unittest
{
    real f3(real x) { return 1.0/(1 + x*x); }
    real f4(real x) { return exp(-x)/sqrt(x); }

    real[] aw = new real[8000];
    real i, err;
    real tiny = 1e-307;
    real eps = 1e-15;

    intdeiini(aw.length, tiny, eps, aw.ptr);

    intdei(&f3, 0.0L, aw.ptr, &i, &err);
    assert (approxEqual(i, PI_2, eps));

    intdei(&f4, 0.0L, aw.ptr, &i, &err);
    assert (approxEqual(i, sqrt(PI), eps));
}




/** I = integral of f(x) over (a,infinity), 
    f(x) has oscillatory factor :
    ---
    f(x) = g(x) * sin(omega * x + theta) as x -> infinity
    ---

    Examples:
        ---
        auto aw = new real[8000];
        intdeoini(aw.length, tiny, eps, aw.ptr);  // initialization of aw
        ...
        intdeo(f, a, omega, aw.ptr, &i, &err);
        ---

    Params:
        lenaw     = length of aw
        tiny      = minimum value that 1/tiny does not 
                    overflow
        eps       = relative error requested
        aw        = points and weights of the quadrature 
                    formula, aw[0...lenaw-1]
        f         = integrand f(x)
        a         = lower limit of integration
        omega     = frequency of oscillation
        i         = approximation to the integral
        err       = estimate of the absolute error

    Remarks:
    <pre>
        initial parameters
            lenaw > 1000, 
            IEEE double :
                lenaw = 8000;
                tiny = 1.0e-307;
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
            err < 0  : abnormal termination.
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

void intdeoini(Real)(size_t lenaw, Real tiny, Real eps, Real *aw)
{
    assert (lenaw > 1000);
    /* ---- adjustable parameter ---- */
    enum LMAX = 5;
    enum : Real { EFS = 0.1, ENOFF = 0.40, PQOFF = 2.9, PPOFF = -0.72 }
    /* ------------------------------ */
    int noff0, nk0, noff, k, nk, j;
    Real tinyln, epsln, frq4, per2, pp, pq, ehp, ehm, h, t, 
        ep, em, tk, xw, wg, xa;
    
    tinyln = -log(tiny);
    epsln = 1 - log(EFS * eps);
    frq4 = M_2_PI;
    per2 = PI;
    pq = PQOFF / epsln;
    pp = PPOFF - log(pq * pq * frq4);
    ehp = exp(2 * pq);
    ehm = 1 / ehp;
    aw[3] = LMAX;
    aw[4] = eps;
    aw[5] = sqrt(EFS * eps);
    noff0 = 6;
    nk0 = 1 + cast(int) (ENOFF * epsln);
    aw[1] = nk0;
    noff = 2 * nk0 + noff0;
    wg = 0;
    xw = 1;
    for (k = 1; k <= nk0; k++) {
        wg += xw;
        aw[noff - 2 * k] = wg;
        aw[noff - 2 * k + 1] = xw;
        xw = xw * (nk0 - k) / k;
    }
    wg = per2 / wg;
    for (k = noff0; k <= noff - 2; k += 2) {
        aw[k] *= wg;
        aw[k + 1] *= wg;
    }
    xw = exp(pp - PI_2);
    aw[noff] = sqrt(xw * (per2 * 0.5));
    aw[noff + 1] = xw * pq;
    aw[noff + 2] = per2 * 0.5;
    h = 2;
    nk = 0;
    k = noff + 3;
    do {
        t = h * 0.5;
        do {
            em = exp(2 * pq * t);
            ep = PI_4 * em;
            em = PI_4 / em;
            tk = t;
            j = k;
            do {
                xw = exp(pp - ep - em);
                wg = sqrt(frq4 * xw + tk * tk);
                xa = xw / (tk + wg);
                wg = (pq * xw * (ep - em) + xa) / wg;
                aw[j] = xa;
                aw[j + 1] = xw * pq;
                aw[j + 2] = wg;
                ep *= ehp;
                em *= ehm;
                tk += 1;
                j += 3;
            } while (ep < tinyln && j <= lenaw - 3);
            t += h;
            k += nk;
        } while (t < 1);
        h *= 0.5;
        if (nk == 0) {
            if (j > lenaw - 6) j -= 3;
            nk = j - noff;
            k += nk;
            aw[2] = nk;
        }
    } while (2 * k - noff - 3 <= lenaw);
    aw[0] = k - 3;
}


/// ditto
void intdeo(Real, Func)(Func f, Real a, Real omega, Real *aw, 
    Real *i, Real *err)
{
    int lenawm, nk0, noff0, nk, noff, lmax, m, k, j, jm, l;
    Real eps, per, perw, w02, ir, h, iback, irback, t, tk, 
        xa, fm, fp, errh, s0, s1, s2, errd;
    
    lenawm = cast(int) (aw[0] + 0.5);
    nk0 = cast(int) (aw[1] + 0.5);
    noff0 = 6;
    nk = cast(int) (aw[2] + 0.5);
    noff = 2 * nk0 + noff0;
    lmax = cast(int) (aw[3] + 0.5);
    eps = aw[4];
    per = 1 / fabs(omega);
    w02 = 2 * aw[noff + 2];
    perw = per * w02;
    *i = f(a + aw[noff] * per);
    ir = *i * aw[noff + 1];
    *i *= aw[noff + 2];
    *err = fabs(*i);
    h = 2;
    m = 1;
    k = noff;
    do {
        iback = *i;
        irback = ir;
        t = h * 0.5;
        do {
            if (k == noff) {
                tk = 1;
                k += nk;
                j = noff;
                do {
                    j += 3;
                    xa = per * aw[j];
                    fm = f(a + xa);
                    fp = f(a + xa + perw * tk);
                    ir += (fm + fp) * aw[j + 1];
                    fm *= aw[j + 2];
                    fp *= w02 - aw[j + 2];
                    *i += fm + fp;
                    *err += fabs(fm) + fabs(fp);
                    tk += 1;
                } while (aw[j] > eps && j < k);
                errh = *err * aw[5];
                *err *= eps;
                jm = j - noff;
            } else {
                tk = t;
                for (j = k + 3; j <= k + jm; j += 3) {
                    xa = per * aw[j];
                    fm = f(a + xa);
                    fp = f(a + xa + perw * tk);
                    ir += (fm + fp) * aw[j + 1];
                    fm *= aw[j + 2];
                    fp *= w02 - aw[j + 2];
                    *i += fm + fp;
                    tk += 1;
                }
                j = k + jm;
                k += nk;
            }
            while (fabs(fm) > *err && j < k) {
                j += 3;
                fm = f(a + per * aw[j]);
                ir += fm * aw[j + 1];
                fm *= aw[j + 2];
                *i += fm;
            }
            fm = f(a + perw * tk);
            s2 = w02 * fm;
            *i += s2;
            if (fabs(fp) > *err || fabs(s2) > *err) {
                l = 0;
                for (;;) {
                    l++;
                    s0 = 0;
                    s1 = 0;
                    s2 = fm * aw[noff0 + 1];
                    for (j = noff0 + 2; j <= noff - 2; j += 2) {
                        tk += 1;
                        fm = f(a + perw * tk);
                        s0 += fm;
                        s1 += fm * aw[j];
                        s2 += fm * aw[j + 1];
                    }
                    if (s2 <= *err || l >= lmax) break;
                    *i += w02 * s0;
                }
                *i += s1;
                if (s2 > *err) *err = s2;
            }
            t += h;
        } while (t < 1);
        if (m == 1) {
            errd = 1 + 2 * errh;
        } else {
            errd = h * (fabs(*i - 2 * iback) + fabs(ir - 2 * irback));
        }
        h *= 0.5;
        m *= 2;
    } while (errd > errh && 2 * k - noff <= lenawm);
    *i *= h * per;
    if (errd > errh) {
        *err = -errd * per;
    } else {
        *err *= per * m * 0.5;
    }
}


unittest
{
    real f5(real x) { return sin(x)/x; }
    real f6(real x) { return cos(x)/sqrt(x); }

    real[] aw = new real[8000];
    real i, err;
    real tiny = 1e-307;
    real eps = 1e-15;

    intdeoini(aw.length, tiny, eps, aw.ptr);

    intdeo(&f5, 0.0L, 1.0L, aw.ptr, &i, &err);
    assert (approxEqual(i, PI_2, eps));

    intdeo(&f6, 0.0L, 1.0L, aw.ptr, &i, &err);
    assert (approxEqual(i, sqrt(PI_2), eps));
}
