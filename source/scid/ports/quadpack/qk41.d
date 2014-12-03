// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.
// An idiomatic D port can be found in scid.internal.calculus.integrate_qk.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qk41;


import std.algorithm: max, min;
import std.math;

import scid.core.fortran;




///
void qk41(Real, Func)(Func f, Real a, Real b, out Real result, out Real abserr,
    out Real resabs, out Real resasc)
{
//***begin prologue  dqk41
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a1a2
//***keywords  41-point gauss-kronrod rules
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
//                       result is computed by applying the 41-point
//                       gauss-kronrod rule (resk) obtained by optimal
//                       addition of abscissae to the 20-point gauss
//                       rule (resg).
//
//              abserr - double precision
//                       estimate of the modulus of the absolute error,
//                       which should not exceed abs(i-result)
//
//              resabs - double precision
//                       approximation to the integral j
//
//              resasc - double precision
//                       approximation to the integal of abs(f-i/(b-a))
//                       over (a,b)
//
//***references  (none)
//***routines called  d1mach
//***end prologue  dqk41
//
      Real absc,centr,dhlgth,
        epmach,fc,fsum,fval1,fval2,hlgth,
        resg,resk,reskh,uflow;
      Real[20] fv1_, fv2_;
      int j,jtw,jtwm1;
//
//           the abscissae and weights are given for the interval (-1,1).
//           because of symmetry only the positive abscissae and their
//           corresponding weights are given.
//
//           xgk    - abscissae of the 41-point gauss-kronrod rule
//                    xgk(2), xgk(4), ...  abscissae of the 20-point
//                    gauss rule
//                    xgk(1), xgk(3), ...  abscissae which are optimally
//                    added to the 20-point gauss rule
//
//           wgk    - weights of the 41-point gauss-kronrod rule
//
//           wg     - weights of the 20-point gauss rule
//
//
// gauss quadrature weights and kronron quadrature abscissae and weights
// as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
// bell labs, nov. 1981.
//
      static immutable Real[10] wg_ = [
          0.0176140071_3915211831_1861962351_853,
          0.0406014298_0038694133_1039952274_932,
          0.0626720483_3410906356_9506535187_042,
          0.0832767415_7670474872_4758143222_046,
          0.1019301198_1724043503_6750135480_350,
          0.1181945319_6151841731_2377377711_382,
          0.1316886384_4917662689_8494499748_163,
          0.1420961093_1838205132_9298325067_165,
          0.1491729864_7260374678_7828737001_969,
          0.1527533871_3072585069_8084331955_098];
//
      static immutable Real[21] xgk_ = [
          0.9988590315_8827766383_8315576545_863,
          0.9931285991_8509492478_6122388471_320,
          0.9815078774_5025025919_3342994720_217,
          0.9639719272_7791379126_7666131197_277,
          0.9408226338_3175475351_9982722212_443,
          0.9122344282_5132590586_7752441203_298,
          0.8782768112_5228197607_7442995113_078,
          0.8391169718_2221882339_4529061701_521,
          0.7950414288_3755119835_0638833272_788,
          0.7463319064_6015079261_4305070355_642,
          0.6932376563_3475138480_5490711845_932,
          0.6360536807_2651502545_2836696226_286,
          0.5751404468_1971031534_2946036586_425,
          0.5108670019_5082709800_4364050955_251,
          0.4435931752_3872510319_9992213492_640,
          0.3737060887_1541956067_2548177024_927,
          0.3016278681_1491300432_0555356858_592,
          0.2277858511_4164507808_0496195368_575,
          0.1526054652_4092267550_5220241022_678,
          0.0765265211_3349733375_4640409398_838,
          0.0000000000_0000000000_0000000000_000];
//
      static immutable Real[21] wgk_ = [
          0.0030735837_1852053150_1218293246_031,
          0.0086002698_5564294219_8661787950_102,
          0.0146261692_5697125298_3787960308_868,
          0.0203883734_6126652359_8010231432_755,
          0.0258821336_0495115883_4505067096_153,
          0.0312873067_7703279895_8543119323_801,
          0.0366001697_5820079803_0557240707_211,
          0.0416688733_2797368626_3788305936_895,
          0.0464348218_6749767472_0231880926_108,
          0.0509445739_2372869193_2707670050_345,
          0.0551951053_4828599474_4832372419_777,
          0.0591114008_8063957237_4967220648_594,
          0.0626532375_5478116802_5870122174_255,
          0.0658345971_3361842211_1563556969_398,
          0.0686486729_2852161934_5623411885_368,
          0.0710544235_5344406830_5790361723_210,
          0.0730306903_3278666749_5189417658_913,
          0.0745828754_0049918898_6581418362_488,
          0.0757044976_8455667465_9542775376_617,
          0.0763778676_7208073670_5502835038_061,
          0.0766007119_1799965644_5049901530_102];
//
      auto fv1 = dimension(fv1_.ptr, 20);
      auto fv2 = dimension(fv2_.ptr, 20);
      auto xgk = dimension(xgk_.ptr, 21);
      auto wgk = dimension(wgk_.ptr, 21);
      auto wg = dimension(wg_.ptr, 10);
//
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the interval
//           hlgth  - half-length of the interval
//           absc   - abscissa
//           fval*  - function value
//           resg   - result of the 20-point gauss formula
//           resk   - result of the 41-point kronrod formula
//           reskh  - approximation to mean value of f over (a,b), i.e.
//                    to i/(b-a)
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqk41
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
      centr = 0.5*(a+b);
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
//
//           compute the 41-point gauss-kronrod approximation to
//           the integral, and estimate the absolute error.
//
      resg = 0.0;
      fc = f(centr);
      resk = wgk[21]*fc;
      resabs = fabs(resk);
      for (j=1; j<=10; j++) { //do 10 j=1,10
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
      for (j=1; j<=10; j++) { //do 15 j = 1,10
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
      resasc = wgk[21]*fabs(fc-reskh);
      for (j=1; j<=20; j++) { //do 20 j=1,20
        resasc = resasc+wgk[j]*(fabs(fv1[j]-reskh)+fabs(fv2[j]-reskh));
 l20: ;}
      result = resk*hlgth;
      resabs = resabs*dhlgth;
      resasc = resasc*dhlgth;
      abserr = fabs((resk-resg)*hlgth);
      if(resasc != 0.0 && abserr != 0.)
        abserr = resasc*min(0.1e1,((cast(Real)0.2e3)*abserr/resasc)^^(cast(Real)1.5));
      if(resabs > uflow/(0.5e2*epmach)) abserr = max
        ((epmach*0.5e2)*resabs,abserr);
      return;
}

version(unittest) import scid.core.testing;
unittest
{
    alias qk41!(float, float delegate(float)) fqk41;
    alias qk41!(double, double delegate(double)) dqk41;
    alias qk41!(double, double function(double)) dfqk41;
    alias qk41!(real, real delegate(real)) rqk41;

    double f(double x) { return x^^3; }
    double result, abserr, resabs, resasc;
    qk41(&f, 0.0, 1.0, result, abserr, resabs, resasc);
    check (isAccurate(result, abserr, 0.25, 1e-6));
}
