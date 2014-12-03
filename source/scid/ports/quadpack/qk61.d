// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.
// An idiomatic D port can be found in scid.internal.calculus.integrate_qk.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qk61;


import std.algorithm: max, min;
import std.math;

import scid.core.fortran;




///
void qk61(Real, Func)(Func f, Real a, Real b, out Real result, out Real abserr,
    out Real resabs, out Real resasc)
{
//***begin prologue  dqk61
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a1a2
//***keywords  61-point gauss-kronrod rules
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  to compute i = integral of f over (a,b) with error
//                           estimate
//                       j = integral of dabs(f) over (a,b)
//***description
//
//        integration rule
//        standard fortran subroutine
//        double precision version
//
//
//        parameters
//         on entry
//           f      - double precision
//                    function subprogram defining the integrand
//                    function f(x). the actual name for f needs to be
//                    declared e x t e r n a l in the calling program.
//
//           a      - double precision
//                    lower limit of integration
//
//           b      - double precision
//                    upper limit of integration
//
//         on return
//           result - double precision
//                    approximation to the integral i
//                    result is computed by applying the 61-point
//                    kronrod rule (resk) obtained by optimal addition of
//                    abscissae to the 30-point gauss rule (resg).
//
//           abserr - double precision
//                    estimate of the modulus of the absolute error,
//                    which should equal or exceed dabs(i-result)
//
//           resabs - double precision
//                    approximation to the integral j
//
//           resasc - double precision
//                    approximation to the integral of dabs(f-i/(b-a))
//
//
//***references  (none)
//***routines called  d1mach
//***end prologue  dqk61
//
      Real dabsc,centr,dhlgth,
        epmach,fc,fsum,fval1,fval2,hlgth,
        resg,resk,reskh,uflow;
      Real[30] fv1_, fv2_;
      int j,jtw,jtwm1;
//
//           the abscissae and weights are given for the
//           interval (-1,1). because of symmetry only the positive
//           abscissae and their corresponding weights are given.
//
//           xgk   - abscissae of the 61-point kronrod rule
//                   xgk(2), xgk(4)  ... abscissae of the 30-point
//                   gauss rule
//                   xgk(1), xgk(3)  ... optimally added abscissae
//                   to the 30-point gauss rule
//
//           wgk   - weights of the 61-point kronrod rule
//
//           wg    - weigths of the 30-point gauss rule
//
//
// gauss quadrature weights and kronron quadrature abscissae and weights
// as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
// bell labs, nov. 1981.
//
      static immutable Real[15] wg_ = [
          0.0079681924_9616660561_5465883474_674,
          0.0184664683_1109095914_2302131912_047,
          0.0287847078_8332336934_9719179611_292,
          0.0387991925_6962704959_6801936446_348,
          0.0484026728_3059405290_2938140422_808,
          0.0574931562_1761906648_1721689402_056,
          0.0659742298_8218049512_8128515115_962,
          0.0737559747_3770520626_8243850022_191,
          0.0807558952_2942021535_4694938460_530,
          0.0868997872_0108297980_2387530715_126,
          0.0921225222_3778612871_7632707087_619,
          0.0963687371_7464425963_9468626351_810,
          0.0995934205_8679526706_2780282103_569,
          0.1017623897_4840550459_6428952168_554,
          0.1028526528_9355884034_1285636705_415];
//
      static immutable Real[31] xgk_ = [
          0.9994844100_5049063757_1325895705_811,
          0.9968934840_7464954027_1630050918_695,
          0.9916309968_7040459485_8628366109_486,
          0.9836681232_7974720997_0032581605_663,
          0.9731163225_0112626837_4693868423_707,
          0.9600218649_6830751221_6871025581_798,
          0.9443744447_4855997941_5831324037_439,
          0.9262000474_2927432587_9324277080_474,
          0.9055733076_9990779854_6522558925_958,
          0.8825605357_9205268154_3116462530_226,
          0.8572052335_4606109895_8658510658_944,
          0.8295657623_8276839744_2898119732_502,
          0.7997278358_2183908301_3668942322_683,
          0.7677774321_0482619491_7977340974_503,
          0.7337900624_5322680472_6171131369_528,
          0.6978504947_9331579693_2292388026_640,
          0.6600610641_2662696137_0053668149_271,
          0.6205261829_8924286114_0477556431_189,
          0.5793452358_2636169175_6024932172_540,
          0.5366241481_4201989926_4169793311_073,
          0.4924804678_6177857499_3693061207_709,
          0.4470337695_3808917678_0609900322_854,
          0.4004012548_3039439253_5476211542_661,
          0.3527047255_3087811347_1037207089_374,
          0.3040732022_7362507737_2677107199_257,
          0.2546369261_6788984643_9805129817_805,
          0.2045251166_8230989143_8957671002_025,
          0.1538699136_0858354696_3794672743_256,
          0.1028069379_6673703014_7096751318_001,
          0.0514718425_5531769583_3025213166_723,
          0.0000000000_0000000000_0000000000_000];
//
      static immutable Real[31] wgk_ = [
          0.0013890136_9867700762_4551591226_760,
          0.0038904611_2709988405_1267201844_516,
          0.0066307039_1593129217_3319826369_750,
          0.0092732796_5951776342_8441146892_024,
          0.0118230152_5349634174_2232898853_251,
          0.0143697295_0704580481_2451432443_580,
          0.0169208891_8905327262_7572289420_322,
          0.0194141411_9394238117_3408951050_128,
          0.0218280358_2160919229_7167485738_339,
          0.0241911620_7808060136_5686370725_232,
          0.0265099548_8233310161_0601709335_075,
          0.0287540487_6504129284_3978785354_334,
          0.0309072575_6238776247_2884252943_092,
          0.0329814470_5748372603_1814191016_854,
          0.0349793380_2806002413_7499670731_468,
          0.0368823646_5182122922_3911065617_136,
          0.0386789456_2472759295_0348651532_281,
          0.0403745389_5153595911_1995279752_468,
          0.0419698102_1516424614_7147541285_970,
          0.0434525397_0135606931_6831728117_073,
          0.0448148001_3316266319_2355551616_723,
          0.0460592382_7100698811_6271735559_374,
          0.0471855465_6929915394_5261478181_099,
          0.0481858617_5708712914_0779492298_305,
          0.0490554345_5502977888_7528165367_238,
          0.0497956834_2707420635_7811569379_942,
          0.0504059214_0278234684_0893085653_585,
          0.0508817958_9874960649_2297473049_805,
          0.0512215478_4925877217_0656282604_944,
          0.0514261285_3745902593_3862879215_781,
          0.0514947294_2945156755_8340433647_099];
//
      auto fv1 = dimension(fv1_.ptr, 30);
      auto fv2 = dimension(fv2_.ptr, 30);
      auto xgk = dimension(xgk_.ptr, 31);
      auto wgk = dimension(wgk_.ptr, 31);
      auto wg = dimension(wg_.ptr, 15);
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the interval
//           hlgth  - half-length of the interval
//           dabsc  - abscissa
//           fval*  - function value
//           resg   - result of the 30-point gauss rule
//           resk   - result of the 61-point kronrod rule
//           reskh  - approximation to the mean value of f
//                    over (a,b), i.e. to i/(b-a)
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
      centr = 0.5*(b+a);
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
//
//           compute the 61-point kronrod approximation to the
//           integral, and estimate the absolute error.
//
//***first executable statement  dqk61
      resg = 0.0;
      fc = f(centr);
      resk = wgk[31]*fc;
      resabs = fabs(resk);
      for (j=1; j<=15; j++) { //do 10 j=1,15
        jtw = j*2;
        dabsc = hlgth*xgk[jtw];
        fval1 = f(centr-dabsc);
        fval2 = f(centr+dabsc);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1+fval2;
        resg = resg+wg[j]*fsum;
        resk = resk+wgk[jtw]*fsum;
        resabs = resabs+wgk[jtw]*(fabs(fval1)+fabs(fval2));
 l10: ;}
      for (j=1; j<=15; j++) { //do 15 j=1,15
        jtwm1 = j*2-1;
        dabsc = hlgth*xgk[jtwm1];
        fval1 = f(centr-dabsc);
        fval2 = f(centr+dabsc);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1+fval2;
        resk = resk+wgk[jtwm1]*fsum;
        resabs = resabs+wgk[jtwm1]*(fabs(fval1)+fabs(fval2));
 l15: ;}
      reskh = resk*0.5;
      resasc = wgk[31]*fabs(fc-reskh);
      for (j=1; j<=30; j++) { //do 20 j=1,30
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
    alias qk61!(float, float delegate(float)) fqk61;
    alias qk61!(double, double delegate(double)) dqk61;
    alias qk61!(double, double function(double)) dfqk61;
    alias qk61!(real, real delegate(real)) rqk61;

    double f(double x) { return x^^3; }
    double result, abserr, resabs, resasc;
    qk61(&f, 0.0, 1.0, result, abserr, resabs, resasc);
    check (isAccurate(result, abserr, 0.25, 1e-6));
}
