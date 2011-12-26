// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qk15i;


import std.algorithm: max, min;
import std.math;

import scid.core.fortran;
import scid.core.meta;




///
void qk15i(Real, Func)(Func f, Real boun, int inf, Real a, Real b,
    out Real result, out Real abserr, out Real resabs, out Real resasc)
{
//***begin prologue  dqk15i
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a3a2,h2a4a2
//***keywords  15-point transformed gauss-kronrod rules
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the original (infinite integration range is mapped
//            onto the interval (0,1) and (a,b) is a part of (0,1).
//            it is the purpose to compute
//            i = integral of transformed integrand over (a,b),
//            j = integral of abs(transformed integrand) over (a,b).
//***description
//
//           integration rule
//           standard fortran subroutine
//           double precision version
//
//           parameters
//            on entry
//              f      - double precision
//                       fuction subprogram defining the integrand
//                       function f(x). the actual name for f needs to be
//                       declared e x t e r n a l in the calling program.
//
//              boun   - double precision
//                       finite bound of original integration
//                       range (set to zero if inf = +2)
//
//              inf    - integer
//                       if inf = -1, the original interval is
//                                   (-infinity,bound),
//                       if inf = +1, the original interval is
//                                   (bound,+infinity),
//                       if inf = +2, the original interval is
//                                   (-infinity,+infinity) and
//                       the integral is computed as the sum of two
//                       integrals, one over (-infinity,0) and one over
//                       (0,+infinity).
//
//              a      - double precision
//                       lower limit for integration over subrange
//                       of (0,1)
//
//              b      - double precision
//                       upper limit for integration over subrange
//                       of (0,1)
//
//            on return
//              result - double precision
//                       approximation to the integral i
//                       result is computed by applying the 15-point
//                       kronrod rule(resk) obtained by optimal addition
//                       of abscissae to the 7-point gauss rule(resg).
//
//              abserr - double precision
//                       estimate of the modulus of the absolute error,
//                       which should equal or exceed abs(i-result)
//
//              resabs - double precision
//                       approximation to the integral j
//
//              resasc - double precision
//                       approximation to the integral of
//                       abs((transformed integrand)-i/(b-a)) over (a,b)
//
//***references  (none)
//***routines called  d1mach
//***end prologue  dqk15i
//
      Real absc,absc1,absc2,centr,dabs,dinf,
        epmach,fc,fsum,fval1,fval2,hlgth,
        resg,resk,reskh,tabsc1,tabsc2,uflow;
      Real[7] fv1_, fv2_;
      int j;
//
//           the abscissae and weights are supplied for the interval
//           (-1,1).  because of symmetry only the positive abscissae and
//           their corresponding weights are given.
//
//           xgk    - abscissae of the 15-point kronrod rule
//                    xgk(2), xgk(4), ... abscissae of the 7-point
//                    gauss rule
//                    xgk(1), xgk(3), ...  abscissae which are optimally
//                    added to the 7-point gauss rule
//
//           wgk    - weights of the 15-point kronrod rule
//
//           wg     - weights of the 7-point gauss rule, corresponding
//                    to the abscissae xgk(2), xgk(4), ...
//                    wg(1), wg(3), ... are set to zero.
//
      static immutable Real[8] wg_ = [
        0.0,
        0.1294849661_6886969327_0611432679_082,
        0.0,
        0.2797053914_8927666790_1467771423_780,
        0.0,
        0.3818300505_0511894495_0369775488_975,
        0.0,
        0.4179591836_7346938775_5102040816_327
      ];
//
      static immutable Real[8] xgk_ = [
        0.9914553711_2081263920_6854697526_329,
        0.9491079123_4275852452_6189684047_851,
        0.8648644233_5976907278_9712788640_926,
        0.7415311855_9939443986_3864773280_788,
        0.5860872354_6769113029_4144838258_730,
        0.4058451513_7739716690_6606412076_961,
        0.2077849550_0789846760_0689403773_245,
        0.0000000000_0000000000_0000000000_000
      ];
//
      static immutable Real[8] wgk_ = [
        0.0229353220_1052922496_3732008058_970,
        0.0630920926_2997855329_0700663189_204,
        0.1047900103_2225018383_9876322541_518,
        0.1406532597_1552591874_5189590510_238,
        0.1690047266_3926790282_6583426598_550,
        0.1903505780_6478540991_3256402421_014,
        0.2044329400_7529889241_4161999234_649,
        0.2094821410_8472782801_2999174891_714
      ];
//
      auto fv1 = dimension(fv1_.ptr, 7);
      auto fv2 = dimension(fv2_.ptr, 7);
      auto xgk = dimension(xgk_.ptr, 8);
      auto wgk = dimension(wgk_.ptr, 8);
      auto wg  = dimension(wg_.ptr, 8);
//
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the interval
//           hlgth  - half-length of the interval
//           absc*  - abscissa
//           tabsc* - transformed abscissa
//           fval*  - function value
//           resg   - result of the 7-point gauss formula
//           resk   - result of the 15-point kronrod formula
//           reskh  - approximation to the mean value of the transformed
//                    integrand over (a,b), i.e. to i/(b-a)
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqk15i
      epmach = Real.epsilon;
      uflow = Real.min_normal;
      dinf = min(1,inf);
//
      centr = 0.5*(a+b);
      hlgth = 0.5*(b-a);
      tabsc1 = boun+dinf*(0.1e1-centr)/centr;
      fval1 = f(tabsc1);
      if(inf == 2) fval1 = fval1+f(-tabsc1);
      fc = (fval1/centr)/centr;
//
//           compute the 15-point kronrod approximation to
//           the integral, and estimate the error.
//
      resg = wg[8]*fc;
      resk = wgk[8]*fc;
      resabs = fabs(resk);
      for (j=1; j<=7; j++) { // end: 10
        absc = hlgth*xgk[j];
        absc1 = centr-absc;
        absc2 = centr+absc;
        tabsc1 = boun+dinf*(0.1e1-absc1)/absc1;
        tabsc2 = boun+dinf*(0.1e1-absc2)/absc2;
        fval1 = f(tabsc1);
        fval2 = f(tabsc2);
        if(inf == 2) fval1 = fval1+f(-tabsc1);
        if(inf == 2) fval2 = fval2+f(-tabsc2);
        fval1 = (fval1/absc1)/absc1;
        fval2 = (fval2/absc2)/absc2;
        fv1[j] = fval1;
        fv2[j] = fval2;
        fsum = fval1+fval2;
        resg = resg+wg[j]*fsum;
        resk = resk+wgk[j]*fsum;
        resabs = resabs+wgk[j]*(fabs(fval1)+fabs(fval2));
 l10: ;}
      reskh = resk*0.5;
      resasc = wgk[8]*fabs(fc-reskh);
      for (j=1; j<=7; j++) { // end: 20
        resasc = resasc+wgk[j]*(fabs(fv1[j]-reskh)+fabs(fv2[j]-reskh));
 l20: ;}
      result = resk*hlgth;
      resasc = resasc*hlgth;
      resabs = resabs*hlgth;
      abserr = fabs((resk-resg)*hlgth);
      if(resasc != 0.0 && abserr != 0.0) abserr = resasc*
        min(0.1e1, ((cast(Real)0.2e3)*abserr/resasc)^^(cast(Real)1.5));
      if(resabs > uflow/(0.5e2*epmach)) abserr = max
       ((epmach*0.5e2)*resabs,abserr);
      return;
}


unittest
{
    alias qk15i!(float, float delegate(float)) fqk15i;
    alias qk15i!(double, double delegate(double)) dqk15i;
    alias qk15i!(double, double function(double)) dfqk15i;
    alias qk15i!(real, real delegate(real)) rqk15i;
}

