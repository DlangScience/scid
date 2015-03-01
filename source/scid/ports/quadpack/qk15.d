// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.
// An idiomatic D port can be found in scid.internal.calculus.integrate_qk.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qk15;


import std.algorithm: max, min;
import std.math;

import scid.core.fortran;




///
void qk15(Real, Func)(Func f, Real a, Real b, out Real result, out Real abserr,
    out Real resabs, out Real resasc)
{
//***begin prologue  dqk15
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a1a2
//***keywords  15-point gauss-kronrod rules
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div - k.u.leuven
//***purpose  to compute i = integral of f over (a,b), with error
//                           estimate
//                       j = integral of abs(f) over (a,b)
//***description
//
//           integration rules
//           standard fortran subroutine
//           double precision version
//
//           parameters
//            on entry
//              f      - double precision
//                       function subprogram defining the integrand
//                       function f(x). the actual name for f needs to be
//                       declared e x t e r n a l in the calling program.
//
//              a      - double precision
//                       lower limit of integration
//
//              b      - double precision
//                       upper limit of integration
//
//            on return
//              result - double precision
//                       approximation to the integral i
//                       result is computed by applying the 15-point
//                       kronrod rule (resk) obtained by optimal addition
//                       of abscissae to the7-point gauss rule(resg).
//
//              abserr - double precision
//                       estimate of the modulus of the absolute error,
//                       which should not exceed abs(i-result)
//
//              resabs - double precision
//                       approximation to the integral j
//
//              resasc - double precision
//                       approximation to the integral of abs(f-i/(b-a))
//                       over (a,b)
//
//***references  (none)
//***routines called  d1mach
//***end prologue  dqk15
//
      Real absc,centr,dhlgth,
        epmach,fc,fsum,fval1,fval2,hlgth,
        resg,resk,reskh,uflow;
      Real[7] fv1_, fv2_;
      int j=1,jtw=1,jtwm1=1;
//
//           the abscissae and weights are given for the interval (-1,1).
//           because of symmetry only the positive abscissae and their
//           corresponding weights are given.
//
//           xgk    - abscissae of the 15-point kronrod rule
//                    xgk(2), xgk(4), ...  abscissae of the 7-point
//                    gauss rule
//                    xgk(1), xgk(3), ...  abscissae which are optimally
//                    added to the 7-point gauss rule
//
//           wgk    - weights of the 15-point kronrod rule
//
//           wg     - weights of the 7-point gauss rule
//
//
// gauss quadrature weights and kronron quadrature abscissae and weights
// as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
// bell labs, nov. 1981.
//
      static immutable Real[4] wg_ = [
        0.1294849661_6886969327_0611432679_082,
        0.2797053914_8927666790_1467771423_780,
        0.3818300505_0511894495_0369775488_975,
        0.4179591836_7346938775_5102040816_327];
//
      static immutable Real[8] xgk_ = [
        0.9914553711_2081263920_6854697526_329,
        0.9491079123_4275852452_6189684047_851,
        0.8648644233_5976907278_9712788640_926,
        0.7415311855_9939443986_3864773280_788,
        0.5860872354_6769113029_4144838258_730,
        0.4058451513_7739716690_6606412076_961,
        0.2077849550_0789846760_0689403773_245,
        0.0000000000_0000000000_0000000000_000];
//
      static immutable Real[8] wgk_ = [
        0.0229353220_1052922496_3732008058_970,
        0.0630920926_2997855329_0700663189_204,
        0.1047900103_2225018383_9876322541_518,
        0.1406532597_1552591874_5189590510_238,
        0.1690047266_3926790282_6583426598_550,
        0.1903505780_6478540991_3256402421_014,
        0.2044329400_7529889241_4161999234_649,
        0.2094821410_8472782801_2999174891_714];
//
      auto fv1 = dimension(fv1_.ptr, 7);
      auto fv2 = dimension(fv2_.ptr, 7);
      auto wg = dimension(wg_.ptr, 4);
      auto wgk = dimension(wgk_.ptr, 8);
      auto xgk = dimension(xgk_.ptr, 8);
//
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the interval
//           hlgth  - half-length of the interval
//           absc   - abscissa
//           fval*  - function value
//           resg   - result of the 7-point gauss formula
//           resk   - result of the 15-point kronrod formula
//           reskh  - approximation to the mean value of f over (a,b),
//                    i.e. to i/(b-a)
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqk15
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
      centr = 0.5*(a+b);
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
//
//           compute the 15-point kronrod approximation to
//           the integral, and estimate the absolute error.
//
      fc = f(centr);
      resg = fc*wg[4];
      resk = fc*wgk[8];
      resabs = fabs(resk);
      for (j=1; j<=3; j++) { //do 10 j=1,3
        jtw = j*2;
        absc = hlgth*xgk[jtw];
        fval1 = f(centr-absc);
        fval2 = f(centr+absc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1+fval2;
        resg = resg+wg[j]*fsum;
        resk = resk+wgk[jtw]*fsum;
        resabs = resabs+wgk[jtw]*(fabs(fval1)+fabs(fval2));
 l10: ;}
      for (j=1; j<=4; j++) { //do 15 j = 1,4
        jtwm1 = j*2-1;
        absc = hlgth*xgk[jtwm1];
        fval1 = f(centr-absc);
        fval2 = f(centr+absc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1+fval2;
        resk = resk+wgk[jtwm1]*fsum;
        resabs = resabs+wgk[jtwm1]*(fabs(fval1)+fabs(fval2));
 l15: ;}
      reskh = resk*0.5;
      resasc = wgk[8]*fabs(fc-reskh);
      for (j=1; j<=7; j++) { //do 20 j=1,7
        resasc = resasc+wgk[j]*(fabs(fv1[j]-reskh)+fabs(fv2[j]-reskh));
 l20: ;}
      result = resk*hlgth;
      resabs = resabs*dhlgth;
      resasc = resasc*dhlgth;
      abserr = fabs((resk-resg)*hlgth);
      if(resasc != 0.0 && abserr != 0.0)
        abserr = resasc*min(0.1e1,((cast(Real)0.2e3)*abserr/resasc)^^(cast(Real)1.5));
      if(resabs > uflow/(0.5e2*epmach)) abserr = max
        ((epmach*0.5e2)*resabs,abserr);
      return;
}

version(unittest) import scid.core.testing;
unittest
{
    alias qk15!(float, float delegate(float)) fqk15;
    alias qk15!(double, double delegate(double)) dqk15;
    alias qk15!(double, double function(double)) dfqk15;
    alias qk15!(real, real delegate(real)) rqk15;

    double f(double x) { return x^^3; }
    double result, abserr, resabs, resasc;
    qk15(&f, 0.0, 1.0, result, abserr, resabs, resasc);
    assert (isAccurate(result, abserr, 0.25, 1e-6));
}
