// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.
// An idiomatic D port can be found in scid.internal.calculus.integrate_qk.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qk31;


import std.algorithm: max, min;
import std.math;

import scid.core.fortran;




///
void qk31(Real, Func)(Func f, Real a, Real b, out Real result,
    out Real abserr, out Real resabs, out Real resasc)
{
//***begin prologue  dqk31
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a1a2
//***keywords  31-point gauss-kronrod rules
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  to compute i = integral of f over (a,b) with error
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
//                       result is computed by applying the 31-point
//                       gauss-kronrod rule (resk), obtained by optimal
//                       addition of abscissae to the 15-point gauss
//                       rule (resg).
//
//              abserr - double precison
//                       estimate of the modulus of the modulus,
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
//***end prologue  dqk31
      Real absc,centr,dhlgth,
        epmach,fc,fsum,fval1,fval2,hlgth,
        resg,resk,reskh,uflow;
      Real[15] fv1_, fv2_;
      int j=1,jtw=1,jtwm1=1;
//
//           the abscissae and weights are given for the interval (-1,1).
//           because of symmetry only the positive abscissae and their
//           corresponding weights are given.
//
//           xgk    - abscissae of the 31-point kronrod rule
//                    xgk(2), xgk(4), ...  abscissae of the 15-point
//                    gauss rule
//                    xgk(1), xgk(3), ...  abscissae which are optimally
//                    added to the 15-point gauss rule
//
//           wgk    - weights of the 31-point kronrod rule
//
//           wg     - weights of the 15-point gauss rule
//
//
// gauss quadrature weights and kronron quadrature abscissae and weights
// as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
// bell labs, nov. 1981.
//
      immutable static Real[8] wg_ = [
          0.0307532419_9611726835_4628393577_204,
          0.0703660474_8810812470_9267416450_667,
          0.1071592204_6717193501_1869546685_869,
          0.1395706779_2615431444_7804794511_028,
          0.1662692058_1699393355_3200860481_209,
          0.1861610000_1556221102_6800561866_423,
          0.1984314853_2711157645_6118326443_839,
          0.2025782419_2556127288_0620199967_519];
//
      immutable static Real[16] xgk_ = [
          0.9980022986_9339706028_5172840152_271,
          0.9879925180_2048542848_9565718586_613,
          0.9677390756_7913913425_7347978784_337,
          0.9372733924_0070590430_7758947710_209,
          0.8972645323_4408190088_2509656454_496,
          0.8482065834_1042721620_0648320774_217,
          0.7904185014_4246593296_7649294817_947,
          0.7244177313_6017004741_6186054613_938,
          0.6509967412_9741697053_3735895313_275,
          0.5709721726_0853884753_7226737253_911,
          0.4850818636_4023968069_3655740232_351,
          0.3941513470_7756336989_7207370981_045,
          0.2991800071_5316881216_6780024266_389,
          0.2011940939_9743452230_0628303394_596,
          0.1011420669_1871749902_7074231447_392,
          0.0000000000_0000000000_0000000000_000];
//
      immutable static Real[16] wgk_ = [
          0.0053774798_7292334898_7792051430_128,
          0.0150079473_2931612253_8374763075_807,
          0.0254608473_2671532018_6874001019_653,
          0.0353463607_9137584622_2037948478_360,
          0.0445897513_2476487660_8227299373_280,
          0.0534815246_9092808726_5343147239_430,
          0.0620095678_0067064028_5139230960_803,
          0.0698541213_1872825870_9520077099_147,
          0.0768496807_5772037889_4432777482_659,
          0.0830805028_2313302103_8289247286_104,
          0.0885644430_5621177064_7275443693_774,
          0.0931265981_7082532122_5486872747_346,
          0.0966427269_8362367850_5179907627_589,
          0.0991735987_2179195933_2393173484_603,
          0.1007698455_2387559504_4946662617_570,
          0.1013300070_1479154901_7374792767_493];
//
      auto fv1 = dimension(fv1_.ptr, 15);
      auto fv2 = dimension(fv2_.ptr, 15);
      auto wg = dimension(wg_.ptr, 8);
      auto wgk = dimension(wgk_.ptr, 16);
      auto xgk = dimension(xgk_.ptr, 16);
//
//
//           list of major variables
//           -----------------------
//           centr  - mid point of the interval
//           hlgth  - half-length of the interval
//           absc   - abscissa
//           fval*  - function value
//           resg   - result of the 15-point gauss formula
//           resk   - result of the 31-point kronrod formula
//           reskh  - approximation to the mean value of f over (a,b),
//                    i.e. to i/(b-a)
//
//           machine dependent constants
//           ---------------------------
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//***first executable statement  dqk31
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
      centr = 0.5*(a+b);
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
//
//           compute the 31-point kronrod approximation to
//           the integral, and estimate the absolute error.
//
      fc = f(centr);
      resg = wg[8]*fc;
      resk = wgk[16]*fc;
      resabs = fabs(resk);
      for (j=1; j<=7; j++) { //do 10 j=1,7
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
      for (j=1; j<=8; j++) { //do 15 j = 1,8
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
      resasc = wgk[16]*fabs(fc-reskh);
      for (j=1; j<=15; j++) { //do 20 j=1,15
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
    alias qk31!(float, float delegate(float)) fqk31;
    alias qk31!(double, double delegate(double)) dqk31;
    alias qk31!(double, double function(double)) dfqk31;
    alias qk31!(real, real delegate(real)) rqk31;

    double f(double x) { return x^^3; }
    double result, abserr, resabs, resasc;
    qk31(&f, 0.0, 1.0, result, abserr, resabs, resasc);
    assert (isAccurate(result, abserr, 0.25, 1e-6));
}
