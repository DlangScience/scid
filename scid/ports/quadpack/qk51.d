/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qk51;


import std.algorithm: max, min;
import std.math;

import scid.common.fortran;




///
void qk51(Real, Func)(Func f, Real a, Real b, out Real result, out Real abserr,
    out Real resabs, out Real resasc)
{
//***begin prologue  dqk51
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a1a2
//***keywords  51-point gauss-kronrod rules
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math & progr. div. - k.u.leuven
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
//                       function subroutine defining the integrand
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
//                       result is computed by applying the 51-point
//                       kronrod rule (resk) obtained by optimal addition
//                       of abscissae to the 25-point gauss rule (resg).
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
//***end prologue  dqk51
//
      Real absc,centr,dhlgth,
        epmach,fc,fsum,fval1,fval2,hlgth,
        resg,resk,reskh,uflow;
      Real[25] fv1_, fv2_;
      int j,jtw,jtwm1;
//
//           the abscissae and weights are given for the interval (-1,1).
//           because of symmetry only the positive abscissae and their
//           corresponding weights are given.
//
//           xgk    - abscissae of the 51-point kronrod rule
//                    xgk(2), xgk(4), ...  abscissae of the 25-point
//                    gauss rule
//                    xgk(1), xgk(3), ...  abscissae which are optimally
//                    added to the 25-point gauss rule
//
//           wgk    - weights of the 51-point kronrod rule
//
//           wg     - weights of the 25-point gauss rule
//
//
// gauss quadrature weights and kronron quadrature abscissae and weights
// as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
// bell labs, nov. 1981.
//
      static immutable Real[13] wg_ = [
          0.0113937985_0102628794_7902964113_235,
          0.0263549866_1503213726_1901815295_299,
          0.0409391567_0130631265_5623487711_646,
          0.0549046959_7583519192_5936891540_473,
          0.0680383338_1235691720_7187185656_708,
          0.0801407003_3500101801_3234959669_111,
          0.0910282619_8296364981_1497220702_892,
          0.1005359490_6705064420_2206890392_686,
          0.1085196244_7426365311_6093957050_117,
          0.1148582591_4571164833_9325545869_556,
          0.1194557635_3578477222_8178126512_901,
          0.1222424429_9031004168_8959518945_852,
          0.1231760537_2671545120_3902873079_050];
//
      static immutable Real[26] xgk_ = [
          0.9992621049_9260983419_3457486540_341,
          0.9955569697_9049809790_8784946893_902,
          0.9880357945_3407724763_7331014577_406,
          0.9766639214_5951751149_8315386479_594,
          0.9616149864_2584251241_8130033660_167,
          0.9429745712_2897433941_4011169658_471,
          0.9207471152_8170156174_6346084546_331,
          0.8949919978_7827536885_1042006782_805,
          0.8658470652_9327559544_8996969588_340,
          0.8334426287_6083400142_1021108693_570,
          0.7978737979_9850005941_0410904994_307,
          0.7592592630_3735763057_7282865204_361,
          0.7177664068_1308438818_6654079773_298,
          0.6735663684_7346836448_5120633247_622,
          0.6268100990_1031741278_8122681624_518,
          0.5776629302_4122296772_3689841612_654,
          0.5263252843_3471918259_9623778158_010,
          0.4730027314_4571496052_2182115009_192,
          0.4178853821_9303774885_1814394594_572,
          0.3611723058_0938783773_5821730127_641,
          0.3030895389_3110783016_7478909980_339,
          0.2438668837_2098843204_5190362797_452,
          0.1837189394_2104889201_5969888759_528,
          0.1228646926_1071039638_7359818808_037,
          0.0615444830_0568507888_6546392366_797,
          0.0000000000_0000000000_0000000000_000];
//
      static immutable Real[26] wgk_ = [
          0.0019873838_9233031592_6507851882_843,
          0.0055619321_3535671375_8040236901_066,
          0.0094739733_8617415160_7207710523_655,
          0.0132362291_9557167481_3656405846_976,
          0.0168478177_0912829823_1516667536_336,
          0.0204353711_4588283545_6568292235_939,
          0.0240099456_0695321622_0092489164_881,
          0.0274753175_8785173780_2948455517_811,
          0.0307923001_6738748889_1109020215_229,
          0.0340021302_7432933783_6748795229_551,
          0.0371162714_8341554356_0330625367_620,
          0.0400838255_0403238207_4839284467_076,
          0.0428728450_2017004947_6895792439_495,
          0.0455029130_4992178890_9870584752_660,
          0.0479825371_3883671390_6392255756_915,
          0.0502776790_8071567196_3325259433_440,
          0.0523628858_0640747586_4366712137_873,
          0.0542511298_8854549014_4543370459_876,
          0.0559508112_2041231730_8240686382_747,
          0.0574371163_6156783285_3582693939_506,
          0.0586896800_2239420796_1974175856_788,
          0.0597203403_2417405997_9099291932_562,
          0.0605394553_7604586294_5360267517_565,
          0.0611285097_1705304830_5859030416_293,
          0.0614711898_7142531666_1544131965_264,
//       note: wgk (26) was calculated from the values of wgk(1..25)
          0.0615808180_6783293507_8759824240_066];
//
      auto fv1 = dimension(fv1_.ptr, 25);
      auto fv2 = dimension(fv2_.ptr, 25);
      auto xgk = dimension(xgk_.ptr, 26);
      auto wgk = dimension(wgk_.ptr, 26);
      auto wg = dimension(wg_.ptr, 13);
//
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the interval
//           hlgth  - half-length of the interval
//           absc   - abscissa
//           fval*  - function value
//           resg   - result of the 25-point gauss formula
//           resk   - result of the 51-point kronrod formula
//           reskh  - approximation to the mean value of f over (a,b),
//                    i.e. to i/(b-a)
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqk51
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
      centr = 0.5*(a+b);
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
//
//           compute the 51-point kronrod approximation to
//           the integral, and estimate the absolute error.
//
      fc = f(centr);
      resg = wg[13]*fc;
      resk = wgk[26]*fc;
      resabs = fabs(resk);
      for (j=1; j<=12; j++) { //do 10 j=1,12
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
      for (j=1; j<=13; j++) { //do 15 j = 1,13
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
      resasc = wgk[26]*fabs(fc-reskh);
      for (j=1; j<=25; j++) { //do 20 j=1,25
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

version(unittest) import scid.common.testing;
unittest
{
    alias qk51!(float, float delegate(float)) fqk51;
    alias qk51!(double, double delegate(double)) dqk51;
    alias qk51!(double, double function(double)) dfqk51;
    alias qk51!(real, real delegate(real)) rqk51;

    double f(double x) { return x^^3; }
    double result, abserr, resabs, resasc;
    qk51(&f, 0.0, 1.0, result, abserr, resabs, resasc);
    check (isAccurate(result, abserr, 0.25, 1e-6));
}
