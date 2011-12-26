// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.
// An idiomatic D port can be found in scid.internal.calculus.integrate_qng.

/** This module is deprecated in favour of scid.internal.calculus.integrate_qng.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qng;


deprecated:


import std.algorithm: max, min;
import std.conv;
import std.math;
import std.traits;

import scid.core.fortran;
import scid.types;

version(unittest)
{
    import scid.core.testing;
}




///
void qng(Real, Func)(Func f, Real a, Real b, Real epsabs, Real epsrel,
    out Real result, out Real abserr, out int neval, out int ier)
{
//***begin prologue  dqng
//***date written   800101   (yymmdd)
//***revision date  810101   (yymmdd)
//***category no.  h2a1a1
//***keywords  automatic integrator, smooth integrand,
//             non-adaptive, gauss-kronrod(patterson)
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl math & progr. div. - k.u.leuven
//           kahaner,david,nbs - modified (2/82)
//***purpose  the routine calculates an approximation result to a
//            given definite integral i = integral of f over (a,b),
//            hopefully satisfying following claim for accuracy
//            abs(i-result).le.max(epsabs,epsrel*abs(i)).
//***description
//
// non-adaptive integration
// standard fortran subroutine
// double precision version
//
//           f      - double precision
//                    function subprogram defining the integrand function
//                    f(x). the actual name for f needs to be declared
//                    e x t e r n a l in the driver program.
//
//           a      - double precision
//                    lower limit of integration
//
//           b      - double precision
//                    upper limit of integration
//
//           epsabs - double precision
//                    absolute accuracy requested
//           epsrel - double precision
//                    relative accuracy requested
//                    if  epsabs.le.0
//                    and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                    the routine will end with ier = 6.
//
//         on return
//           result - double precision
//                    approximation to the integral i
//                    result is obtained by applying the 21-point
//                    gauss-kronrod rule (res21) obtained by optimal
//                    addition of abscissae to the 10-point gauss rule
//                    (res10), or by applying the 43-point rule (res43)
//                    obtained by optimal addition of abscissae to the
//                    21-point gauss-kronrod rule, or by applying the
//                    87-point rule (res87) obtained by optimal addition
//                    of abscissae to the 43-point rule.
//
//           abserr - double precision
//                    estimate of the modulus of the absolute error,
//                    which should equal or exceed abs(i-result)
//
//           neval  - integer
//                    number of integrand evaluations
//
//           ier    - ier = 0 normal and reliable termination of the
//                            routine. it is assumed that the requested
//                            accuracy has been achieved.
//                    ier.gt.0 abnormal termination of the routine. it is
//                            assumed that the requested accuracy has
//                            not been achieved.
//           error messages
//                    ier = 1 the maximum number of steps has been
//                            executed. the integral is probably too
//                            difficult to be calculated by dqng.
//                        = 6 the input is invalid, because
//                            epsabs.le.0 and
//                            epsrel.lt.max(50*rel.mach.acc.,0.5d-28).
//                            result, abserr and neval are set to zero.
//
//***references  (none)
//***routines called  d1mach,xerror
//***end prologue  dqng
//
      Real absc,centr,dhlgth,
        epmach,fcentr,fval,fval1,fval2,
        hlgth,res10,res21,res43,res87,resabs,resasc,
        reskh,uflow;
      Real[5] fv1_, fv2_, fv3_, fv4_;
      Real[21] savfun_;
      int ipx=1,k=1,l=1;
//
//           the following data statements contain the
//           abscissae and weights of the integration rules used.
//
//           x1      abscissae common to the 10-, 21-, 43- and 87-
//                   point rule
//           x2      abscissae common to the 21-, 43- and 87-point rule
//           x3      abscissae common to the 43- and 87-point rule
//           x4      abscissae of the 87-point rule
//           w10     weights of the 10-point formula
//           w21a    weights of the 21-point formula for abscissae x1
//           w21b    weights of the 21-point formula for abscissae x2
//           w43a    weights of the 43-point formula for abscissae x1, x3
//           w43b    weights of the 43-point formula for abscissae x3
//           w87a    weights of the 87-point formula for abscissae x1,
//                   x2, x3
//           w87b    weights of the 87-point formula for abscissae x4
//
//
// gauss-kronrod-patterson quadrature coefficients for use in
// quadpack routine qng.  these coefficients were calculated with
// 101 decimal digit arithmetic by l. w. fullerton, bell labs, nov 1981.
//
      static immutable Real[5] x1_ = [
          0.9739065285_1717172007_7964012084_452,
          0.8650633666_8898451073_2096688423_493,
          0.6794095682_9902440623_4327365114_874,
          0.4333953941_2924719079_9265943165_784,
          0.1488743389_8163121088_4826001129_720];
      static immutable Real[5] w10_ = [
          0.0666713443_0868813759_3568809893_332,
          0.1494513491_5058059314_5776339657_697,
          0.2190863625_1598204399_5534934228_163,
          0.2692667193_0999635509_1226921569_469,
          0.2955242247_1475287017_3892994651_338];
//
      static immutable Real[5] x2_ = [
          0.9956571630_2580808073_5527280689_003,
          0.9301574913_5570822600_1207180059_508,
          0.7808177265_8641689706_3717578345_042,
          0.5627571346_6860468333_9000099272_694,
          0.2943928627_0146019813_1126603103_866];
      static immutable Real[5] w21a_ = [
          0.0325581623_0796472747_8818972459_390,
          0.0750396748_1091995276_7043140916_190,
          0.1093871588_0229764189_9210590325_805,
          0.1347092173_1147332592_8054001771_707,
          0.1477391049_0133849137_4841515972_068];
      static immutable Real[6] w21b_ = [
          0.0116946388_6737187427_8064396062_192,
          0.0547558965_7435199603_1381300244_580,
          0.0931254545_8369760553_5065465083_366,
          0.1234919762_6206585107_7958109831_074,
          0.1427759385_7706008079_7094273138_717,
          0.1494455540_0291690566_4936468389_821];
//
      static immutable Real[11] x3_ = [
          0.9993333609_0193208139_4099323919_911,
          0.9874334029_0808886979_5961478381_209,
          0.9548079348_1426629925_7919200290_473,
          0.9001486957_4832829362_5099494069_092,
          0.8251983149_8311415084_7066732588_520,
          0.7321483889_8930498261_2354848755_461,
          0.6228479705_3772523864_1159120344_323,
          0.4994795740_7105649995_2214885499_755,
          0.3649016613_4658076804_3989548502_644,
          0.2222549197_7660129649_8260928066_212,
          0.0746506174_6138332204_3914435796_506];
      static immutable Real[10] w43a_ = [
          0.0162967342_8966656492_4281974617_663,
          0.0375228761_2086950146_1613795898_115,
          0.0546949020_5825544214_7212685465_005,
          0.0673554146_0947808607_5553166302_174,
          0.0738701996_3239395343_2140695251_367,
          0.0057685560_5976979618_4184327908_655,
          0.0273718905_9324884208_1276069289_151,
          0.0465608269_1042883074_3339154433_824,
          0.0617449952_0144256449_6240336030_883,
          0.0713872672_6869339776_8559114425_516];
      static immutable Real[12] w43b_ = [
          0.0018444776_4021241410_0389106552_965,
          0.0107986895_8589165174_0465406741_293,
          0.0218953638_6779542810_2523123075_149,
          0.0325974639_7534568944_3882222526_137,
          0.0421631379_3519181184_7627924327_955,
          0.0507419396_0018457778_0189020092_084,
          0.0583793955_4261924837_5475369330_206,
          0.0647464049_5144588554_4689259517_511,
          0.0695661979_1235648452_8633315038_405,
          0.0728244414_7183320815_0939535192_842,
          0.0745077510_1417511827_3571813842_889,
          0.0747221475_1740300559_4425168280_423];
//
      static immutable Real[22] x4_ = [
          0.9999029772_6272923449_0529830591_582,
          0.9979898959_8667874542_7496322365_960,
          0.9921754978_6068722280_8523352251_425,
          0.9813581635_7271277357_1916941623_894,
          0.9650576238_5838461912_8284110607_926,
          0.9431676131_3367059681_6416634507_426,
          0.9158064146_8550720959_1826430720_050,
          0.8832216577_7131650137_2117548744_163,
          0.8457107484_6241566660_5902011504_855,
          0.8035576580_3523098278_8739474980_964,
          0.7570057306_8549555832_8942793432_020,
          0.7062732097_8732181982_4094274740_840,
          0.6515894665_0117792253_4422205016_736,
          0.5932233740_5796108887_5273770349_144,
          0.5314936059_7083193228_5268948562_671,
          0.4667636230_4202284487_1966781659_270,
          0.3994248478_5921880473_2101665817_923,
          0.3298748771_0618828826_5053371824_597,
          0.2585035592_0216155180_2280975429_025,
          0.1856953965_6834665201_5917141167_606,
          0.1118422131_7990746817_2398359241_362,
          0.0373521233_9461987081_4998165437_704];
      static immutable Real[21] w87a_ = [
          0.0081483773_8414917290_0002878448_190,
          0.0187614382_0156282224_3935059003_794,
          0.0273474510_5005228616_1582829741_283,
          0.0336777073_1163793004_6581056957_588,
          0.0369350998_2042790761_4589586742_499,
          0.0028848724_3021153050_1334156248_695,
          0.0136859460_2271270188_8950035273_128,
          0.0232804135_0288831112_3409291030_404,
          0.0308724976_1171335867_5466394126_442,
          0.0356936336_3941877071_9351355457_044,
          0.0009152833_4520224136_0843392549_948,
          0.0053992802_1930047136_7738743391_053,
          0.0109476796_0111893113_4327826856_808,
          0.0162987316_9678733526_2665703223_280,
          0.0210815688_8920383511_2433060188_190,
          0.0253709697_6925382724_3467999831_710,
          0.0291896977_5647575250_1446154084_920,
          0.0323732024_6720278968_5788194889_595,
          0.0347830989_5036514275_0781997949_596,
          0.0364122207_3135178756_2801163687_577,
          0.0372538755_0304770853_9592001191_226];
      static immutable Real[23] w87b_ = [
          0.0002741455_6376207235_0016527092_881,
          0.0018071241_5505794294_8341311753_254,
          0.0040968692_8275916486_4458070683_480,
          0.0067582900_5184737869_9816577897_424,
          0.0095499576_7220164653_6053581325_377,
          0.0123294476_5224485369_4626639963_780,
          0.0150104473_4638895237_6697286041_943,
          0.0175489679_8624319109_9665352925_900,
          0.0199380377_8644088820_2278192730_714,
          0.0221949359_6101228679_6332102959_499,
          0.0243391471_2600080547_0360647041_454,
          0.0263745054_1483920724_1503786552_615,
          0.0282869107_8877120065_9968002987_960,
          0.0300525811_2809269532_2521110347_341,
          0.0316467513_7143992940_4586051078_883,
          0.0330504134_1997850329_0785944862_689,
          0.0342550997_0422606178_7082821046_821,
          0.0352624126_6015668103_3782717998_428,
          0.0360769896_2288870118_5500318003_895,
          0.0366986044_9845609449_8018047441_094,
          0.0371205492_6983257611_4119958413_599,
          0.0373342287_5193504032_1235449094_698,
          0.0373610737_6267902341_0321241766_599];
//
      auto fv1 = dimension(fv1_.ptr, 5);
      auto fv2 = dimension(fv2_.ptr, 5);
      auto fv3 = dimension(fv3_.ptr, 5);
      auto fv4 = dimension(fv4_.ptr, 5);
      auto x1 = dimension(x1_.ptr, 5);
      auto x2 = dimension(x2_.ptr, 5);
      auto x3 = dimension(x3_.ptr, 11);
      auto x4 = dimension(x4_.ptr, 22);
      auto w10 = dimension(w10_.ptr, 5);
      auto w21a = dimension(w21a_.ptr, 5);
      auto w21b = dimension(w21b_.ptr, 6);
      auto w43a = dimension(w43a_.ptr, 10);
      auto w43b = dimension(w43b_.ptr, 12);
      auto w87a = dimension(w87a_.ptr, 21);
      auto w87b = dimension(w87b_.ptr, 23);
      auto savfun = dimension(savfun_.ptr, 21);
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the integration interval
//           hlgth  - half-length of the integration interval
//           fcentr - function value at mid point
//           absc   - abscissa
//           fval   - function value
//           savfun - array of function values which have already been
//                    computed
//           res10  - 10-point gauss result
//           res21  - 21-point kronrod result
//           res43  - 43-point result
//           res87  - 87-point result
//           resabs - approximation to the integral of abs(f)
//           resasc - approximation to the integral of abs(f-i/(b-a))
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqng
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
//           test on validity of parameters
//           ------------------------------
//
      result = 0.0;
      abserr = 0.0;
      neval = 0;
      ier = 6;
      if(epsabs <= 0.0 && epsrel < max(0.5e2*epmach,0.5e-28))
        goto l80;
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
      centr = 0.5*(b+a);
      fcentr = f(centr);
      neval = 21;
      ier = 1;
//
//          compute the integral using the 10- and 21-point formula.
//
      for (l=1; l<=3; l++) { //do 70 l = 1,3
      if (l == 1) goto l5;
      else if (l == 2) goto l25;
      else if (l == 3) goto l45;
      //goto l(5,25,45),l
  l5: res10 = 0.0;
      res21 = w21b[6]*fcentr;
      resabs = w21b[6]*fabs(fcentr);
      for (k=1; k<=5; k++) { //do 10 k=1,5
        absc = hlgth*x1[k];
        fval1 = f(centr+absc);
        fval2 = f(centr-absc);
        fval = fval1+fval2;
        res10 = res10+w10[k]*fval;
        res21 = res21+w21a[k]*fval;
        resabs = resabs+w21a[k]*(fabs(fval1)+fabs(fval2));
        savfun[k] = fval;
        fv1[k] = fval1;
        fv2[k] = fval2;
 l10: ;}
      ipx = 5;
      for (k=1; k<=5; k++) { //do 15 k=1,5
        ipx = ipx+1;
        absc = hlgth*x2[k];
        fval1 = f(centr+absc);
        fval2 = f(centr-absc);
        fval = fval1+fval2;
        res21 = res21+w21b[k]*fval;
        resabs = resabs+w21b[k]*(fabs(fval1)+fabs(fval2));
        savfun[ipx] = fval;
        fv3[k] = fval1;
        fv4[k] = fval2;
 l15: ;}
//
//          test for convergence.
//
      result = res21*hlgth;
      resabs = resabs*dhlgth;
      reskh = 0.5*res21;
      resasc = w21b[6]*fabs(fcentr-reskh);
      for (k=1; k<=5; k++) { //do 20 k = 1,5
        resasc = resasc+w21a[k]*(fabs(fv1[k]-reskh)+fabs(fv2[k]-reskh))
                        +w21b[k]*(fabs(fv3[k]-reskh)+fabs(fv4[k]-reskh));
 l20: ;}
      abserr = fabs((res21-res10)*hlgth);
      resasc = resasc*dhlgth;
      goto l65;
//
//          compute the integral using the 43-point formula.
//
 l25: res43 = w43b[12]*fcentr;
      neval = 43;
      for (k=1; k<=10; k++) { //do 30 k=1,10
        res43 = res43+savfun[k]*w43a[k];
 l30: ;}
      for (k=1; k<=11; k++) { //do 40 k=1,11
        ipx = ipx+1;
        absc = hlgth*x3[k];
        fval = f(absc+centr)+f(centr-absc);
        res43 = res43+fval*w43b[k];
        savfun[ipx] = fval;
 l40: ;}
//
//          test for convergence.
//
      result = res43*hlgth;
      abserr = fabs((res43-res21)*hlgth);
      goto l65;
//
//          compute the integral using the 87-point formula.
//
 l45: res87 = w87b[23]*fcentr;
      neval = 87;
      for (k=1; k<=21; k++) { //do 50 k=1,21
        res87 = res87+savfun[k]*w87a[k];
 l50: ;}
      for (k=1; k<=22; k++) { // do 60 k=1,22
        absc = hlgth*x4[k];
        res87 = res87+w87b[k]*(f(absc+centr)+f(centr-absc));
 l60: ;}
      result = res87*hlgth;
      abserr = fabs((res87-res43)*hlgth);
 l65: if(resasc != 0.0 && abserr != 0.0)
        abserr = resasc*min(0.1e1,((cast(Real)0.2e3)*abserr/resasc)^^(cast(Real)1.5));
      if (resabs > uflow/(0.5e2*epmach)) abserr = max
        ((epmach*0.5e2)*resabs,abserr);
      if (abserr <= max(epsabs,epsrel*fabs(result))) ier = 0;
// ***jump out of do-loop
      if (ier == 0) goto l999;
 l70: ;}
 l80: throw new Exception("abnormal return from qng: "~to!string(ier));
l999: return;
}

unittest
{
    alias qng!(float, float delegate(float)) fqng;
    alias qng!(double, double delegate(double)) dqng;
    alias qng!(double, double function(double)) dfqng;
    alias qng!(real, real delegate(real)) rqng;

    float f(float x) { return x <= 0.0 ? 0.0 : sqrt(x)*log(x); }

    enum : float
    {
        a = 0.0,
        b = 1.0,
        epsabs = 0.0,
        epsrel = 0.001
    }
    float result, abserr;
    int neval, ier;

    qng(&f, a, b, epsabs, epsrel, result, abserr, neval, ier);
    
    check(isAccurate(result, abserr, -4.0f/9.0f, epsrel, epsabs));
}
