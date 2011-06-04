/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qk15w;


import std.algorithm: min, max;
import std.math;

import scid.common.fortran;




///
void qk15w(Real, Func)(Func f, Real function(Real,Real,Real,Real,Real,int) w,
    Real p1, Real p2, Real p3, Real p4, int kp, Real a, Real b,
    out Real result, out Real abserr, out Real resabs, out Real resasc)
{
//***begin prologue  dqk15w
//***date written   810101   (yymmdd)
//***revision date  830518   (mmddyy)
//***category no.  h2a2a2
//***keywords  15-point gauss-kronrod rules
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  to compute i = integral of f*w over (a,b), with error
//                           estimate
//                       j = integral of abs(f*w) over (a,b)
//***description
//
//           integration rules
//           standard fortran subroutine
//           double precision version
//
//           parameters
//             on entry
//              f      - double precision
//                       function subprogram defining the integrand
//                       function f(x). the actual name for f needs to be
//                       declared e x t e r n a l in the driver program.
//
//              w      - double precision
//                       function subprogram defining the integrand
//                       weight function w(x). the actual name for w
//                       needs to be declared e x t e r n a l in the
//                       calling program.
//
//              p1, p2, p3, p4 - double precision
//                       parameters in the weight function
//
//              kp     - integer
//                       key for indicating the type of weight function
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
//                       of abscissae to the 7-point gauss rule (resg).
//
//              abserr - double precision
//                       estimate of the modulus of the absolute error,
//                       which should equal or exceed abs(i-result)
//
//              resabs - double precision
//                       approximation to the integral of abs(f)
//
//              resasc - double precision
//                       approximation to the integral of abs(f-i/(b-a))
//
//
//***references  (none)
//***routines called  d1mach
//***end prologue  dqk15w
//
      Real absc,absc1,absc2,centr,dhlgth,
        epmach,fc,fsum,fval1,fval2,hlgth,
        resg,resk,reskh,uflow;
      Real[7] fv1_, fv2_;
      int j=1,jtw=1,jtwm1=1;
//
//           the abscissae and weights are given for the interval (-1,1).
//           because of symmetry only the positive abscissae and their
//           corresponding weights are given.
//
//           xgk    - abscissae of the 15-point gauss-kronrod rule
//                    xgk(2), xgk(4), ... abscissae of the 7-point
//                    gauss rule
//                    xgk(1), xgk(3), ... abscissae which are optimally
//                    added to the 7-point gauss rule
//
//           wgk    - weights of the 15-point gauss-kronrod rule
//
//           wg     - weights of the 7-point gauss rule
//
      static immutable Real[8] xgk_ = [
           0.9914553711208126,
           0.9491079123427585,
           0.8648644233597691,
           0.7415311855993944,
           0.5860872354676911,
           0.4058451513773972,
           0.2077849550078985,
           0.0000000000000000
        ];
//
     static immutable Real[8] wgk_ = [
           0.2293532201052922e-1,
           0.6309209262997855e-1,
           0.1047900103222502,
           0.1406532597155259,
           0.1690047266392679,
           0.1903505780647854,
           0.2044329400752989,
           0.2094821410847278
        ];
//
      static immutable Real[4] wg_ = [
           0.1294849661688697,
           0.2797053914892767,
           0.3818300505051889,
           0.4179591836734694
        ];
//
      auto fv1 = dimension(fv1_.ptr, 7);
      auto fv2 = dimension(fv2_.ptr, 7);
      auto xgk = dimension(xgk_.ptr, 8);
      auto wgk = dimension(wgk_.ptr, 8);
      auto wg = dimension(wg_.ptr, 4);
//
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the interval
//           hlgth  - half-length of the interval
//           absc*  - abscissa
//           fval*  - function value
//           resg   - result of the 7-point gauss formula
//           resk   - result of the 15-point kronrod formula
//           reskh  - approximation to the mean value of f*w over (a,b),
//                    i.e. to i/(b-a)
//
//           machine dependent constants
//           ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqk15w
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
      centr = 0.5*(a+b);
      hlgth = 0.5*(b-a);
      dhlgth = fabs(hlgth);
//
//           compute the 15-point kronrod approximation to the
//           integral, and estimate the error.
//
      fc = f(centr)*w(centr,p1,p2,p3,p4,kp);
      resg = wg[4]*fc;
      resk = wgk[8]*fc;
      resabs = fabs(resk);
      for (j=1; j<=3; j++) { // end: 10
        jtw = j*2;
        absc = hlgth*xgk[jtw];
        absc1 = centr-absc;
        absc2 = centr+absc;
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp);
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp);
        fv1[jtw] = fval1;
        fv2[jtw] = fval2;
        fsum = fval1+fval2;
        resg = resg+wg[j]*fsum;
        resk = resk+wgk[jtw]*fsum;
        resabs = resabs+wgk[jtw]*(fabs(fval1)+fabs(fval2));
 l10:;}
      for(j=1; j<=4; j++) { // end: 15
        jtwm1 = j*2-1;
        absc = hlgth*xgk[jtwm1];
        absc1 = centr-absc;
        absc2 = centr+absc;
        fval1 = f(absc1)*w(absc1,p1,p2,p3,p4,kp);
        fval2 = f(absc2)*w(absc2,p1,p2,p3,p4,kp);
        fv1[jtwm1] = fval1;
        fv2[jtwm1] = fval2;
        fsum = fval1+fval2;
        resk = resk+wgk[jtwm1]*fsum;
        resabs = resabs+wgk[jtwm1]*(fabs(fval1)+fabs(fval2));
 l15:;}
      reskh = resk*0.5;
      resasc = wgk[8]*fabs(fc-reskh);
      for (j=1; j<=7; j++) { // end: 20
        resasc = resasc+wgk[j]*(fabs(fv1[j]-reskh)+fabs(fv2[j]-reskh));
 l20:;}
      result = resk*hlgth;
      resabs = resabs*dhlgth;
      resasc = resasc*dhlgth;
      abserr = fabs((resk-resg)*hlgth);
      if(resasc != 0.0 && abserr != 0.0)
        abserr = resasc*min(0.1e1, ((cast(Real)0.2e3)*abserr/resasc)^^(cast(Real)1.5));
      if(resabs > uflow/(0.5e2*epmach)) abserr = max((epmach*
        0.5e2)*resabs,abserr);
      return;
}


unittest
{
    alias qk15w!(float, float delegate(float)) fqk15w;
    alias qk15w!(double, double delegate(double)) dqk15w;
    alias qk15w!(double, double function(double)) dfqk15w;
    alias qk15w!(real, real delegate(real)) rqk15w;
}
