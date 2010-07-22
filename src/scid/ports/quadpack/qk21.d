/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qk21;


import std.algorithm: max, min;
import std.math;

import scid.core.fortran;




///
void qk21(Real, Func)(Func f, Real a, Real b, out Real result,
    out Real abserr, out Real resabs, out Real resasc)
{
//***begin prologue  dqk21
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a1a2
//***keywords  21-point gauss-kronrod rules
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
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
//                       declared e x t e r n a l in the driver program.
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
//                       result is computed by applying the 21-point
//                       kronrod rule (resk) obtained by optimal addition
//                       of abscissae to the 10-point gauss rule (resg).
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
//***end prologue  dqk21
//
      Real absc,centr,dhlgth,
        epmach,fc,fsum,fval1,fval2,hlgth,
        resg,resk,reskh,uflow;
      Real[10] fv1_, fv2_;
      int j,jtw,jtwm1;
//
//           the abscissae and weights are given for the interval (-1,1).
//           because of symmetry only the positive abscissae and their
//           corresponding weights are given.
//
//           xgk    - abscissae of the 21-point kronrod rule
//                    xgk(2), xgk(4), ...  abscissae of the 10-point
//                    gauss rule
//                    xgk(1), xgk(3), ...  abscissae which are optimally
//                    added to the 10-point gauss rule
//
//           wgk    - weights of the 21-point kronrod rule
//
//           wg     - weights of the 10-point gauss rule
//
//
// gauss quadrature weights and kronron quadrature abscissae and weights
// as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
// bell labs, nov. 1981.
//
      static immutable Real[5] wg_ = [
        0.0666713443_0868813759_3568809893_332,
        0.1494513491_5058059314_5776339657_697,
        0.2190863625_1598204399_5534934228_163,
        0.2692667193_0999635509_1226921569_469,
        0.2955242247_1475287017_3892994651_338];
//
      static immutable Real[11] xgk_ = [
        0.9956571630_2580808073_5527280689_003,
        0.9739065285_1717172007_7964012084_452,
        0.9301574913_5570822600_1207180059_508,
        0.8650633666_8898451073_2096688423_493,
        0.7808177265_8641689706_3717578345_042,
        0.6794095682_9902440623_4327365114_874,
        0.5627571346_6860468333_9000099272_694,
        0.4333953941_2924719079_9265943165_784,
        0.2943928627_0146019813_1126603103_866,
        0.1488743389_8163121088_4826001129_720,
        0.0000000000_0000000000_0000000000_000];
//
      static immutable Real[11] wgk_ = [
        0.0116946388_6737187427_8064396062_192,
        0.0325581623_0796472747_8818972459_390,
        0.0547558965_7435199603_1381300244_580,
        0.0750396748_1091995276_7043140916_190,
        0.0931254545_8369760553_5065465083_366,
        0.1093871588_0229764189_9210590325_805,
        0.1234919762_6206585107_7958109831_074,
        0.1347092173_1147332592_8054001771_707,
        0.1427759385_7706008079_7094273138_717,
        0.1477391049_0133849137_4841515972_068,
        0.1494455540_0291690566_4936468389_821];
//
      auto fv1 = dimension(fv1_.ptr, 10);
      auto fv2 = dimension(fv2_.ptr, 10);
      auto wg = dimension(wg_.ptr, 5);
      auto wgk = dimension(wgk_.ptr, 11);
      auto xgk = dimension(xgk_.ptr, 11);
//
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the interval
//           hlgth  - half-length of the interval
//           absc   - abscissa
//           fval*  - function value
//           resg   - result of the 10-point gauss formula
//           resk   - result of the 21-point kronrod formula
//           reskh  - approximation to the mean value of f over (a,b),
//                    i.e. to i/(b-a)
//
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqk21
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
      centr = 0.5*(a+b);
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
//
//           compute the 21-point kronrod approximation to
//           the integral, and estimate the absolute error.
//
      resg = 0.0;
      fc = f(centr);
      resk = wgk[11]*fc;
      resabs = fabs(resk);
      for (j=1; j<=5; j++) { //do 10 j = 1,5
        jtw = 2*j;
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
      for (j=1; j<=5; j++) { //do 15 j = 1,5
        jtwm1 = 2*j-1;
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
      resasc = wgk[11]*fabs(fc-reskh);
      for (j=1; j<=10; j++) { //do 20 j=1,10
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
    alias qk21!(float, float delegate(float)) fqk21;
    alias qk21!(double, double delegate(double)) dqk21;
    alias qk21!(double, double function(double)) dfqk21;
    alias qk21!(real, real delegate(real)) rqk21;

    double f(double x) { return x^^3; }
    double result, abserr, resabs, resasc;
    qk21(&f, 0.0, 1.0, result, abserr, resabs, resasc);
    check (isAccurate(result, abserr, 0.25, 1e-6));
}
