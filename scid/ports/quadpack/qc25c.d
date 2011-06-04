/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qc25c;


import std.math: fabs, log;

import scid.common.fortran;
import scid.ports.quadpack.qcheb;
import scid.ports.quadpack.qk15w;
import scid.ports.quadpack.qwgtc;




///
void qc25c(Real, Func)(Func f, Real a, Real b, Real c, out Real result,
    out Real abserr, ref int krul, out int neval)
{
//***begin prologue  dqc25c
//***date written   810101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a2,j4
//***keywords  25-point clenshaw-curtis integration
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  to compute i = integral of f*w over (a,b) with
//            error estimate, where w(x) = 1/(x-c)
//***description
//
//        integration rules for the computation of cauchy
//        principal value integrals
//        standard fortran subroutine
//        double precision version
//
//        parameters
//           f      - double precision
//                    function subprogram defining the integrand function
//                    f(x). the actual name for f needs to be declared
//                    e x t e r n a l  in the driver program.
//
//           a      - double precision
//                    left end point of the integration interval
//
//           b      - double precision
//                    right end point of the integration interval, b.gt.a
//
//           c      - double precision
//                    parameter in the weight function
//
//           result - double precision
//                    approximation to the integral
//                    result is computed by using a generalized
//                    clenshaw-curtis method if c lies within ten percent
//                    of the integration interval. in the other case the
//                    15-point kronrod rule obtained by optimal addition
//                    of abscissae to the 7-point gauss rule, is applied.
//
//           abserr - double precision
//                    estimate of the modulus of the absolute error,
//                    which should equal or exceed abs(i-result)
//
//           krul   - integer
//                    key which is decreased by 1 if the 15-point
//                    gauss-kronrod scheme has been used
//
//           neval  - integer
//                    number of integrand evaluations
//
//.......................................................................
//***references  (none)
//***routines called  dqcheb,dqk15w,dqwgtc
//***end prologue  dqc25c
//
      Real ak22,amom0,amom1,amom2,cc,centr,
        hlgth,p2,p3,p4,resabs,
        resasc,res12,res24,u;
      Real[13] cheb12_;
      Real[25] fval_, cheb24_;
      int i=1,isym=1,k=1,kp=1;
//
//           the vector x contains the values cos(k*pi/24),
//           k = 1, ..., 11, to be used for the chebyshev series
//           expansion of f
//
      static immutable Real[11] x_ = [
          0.9914448613_7381041114_4557526928_563,
          0.9659258262_8906828674_9743199728_897,
          0.9238795325_1128675612_8183189396_788,
          0.8660254037_8443864676_3723170752_936,
          0.7933533402_9123516457_9776961501_299,
          0.7071067811_8654752440_0844362104_849,
          0.6087614290_0872063941_6097542898_164,
          0.5000000000_0000000000_0000000000_000,
          0.3826834323_6508977172_8459984030_399,
          0.2588190451_0252076234_8898837624_048,
          0.1305261922_2005159154_8406227895_489];
//
    auto x = dimension(x_.ptr, 11);
    auto fval = dimension(fval_.ptr, 25);
    auto cheb12 = dimension(cheb12_.ptr, 13);
    auto cheb24 = dimension(cheb24_.ptr, 25);
//
//           list of major variables
//           ----------------------
//           fval   - value of the function f at the points
//                    cos(k*pi/24),  k = 0, ..., 24
//           cheb12 - chebyshev series expansion coefficients,
//                    for the function f, of degree 12
//           cheb24 - chebyshev series expansion coefficients,
//                    for the function f, of degree 24
//           res12  - approximation to the integral corresponding
//                    to the use of cheb12
//           res24  - approximation to the integral corresponding
//                    to the use of cheb24
//           dqwgtc - external function subprogram defining
//                    the weight function
//           hlgth  - half-length of the interval
//           centr  - mid point of the interval
//
//
//           check the position of c.
//
//***first executable statement  dqc25c
      cc = (0.2e1*c-b-a)/(b-a);
      if(fabs(cc) < 0.11e1) goto l10;
//
//           apply the 15-point gauss-kronrod scheme.
//
      krul = krul-1;
      qk15w!(Real,Func)(f,&qwgtc!Real,c,p2,p3,p4,kp,a,b,result,abserr,
        resabs,resasc);
      neval = 15;
      if (resasc == abserr) krul = krul+1;
      goto l50;
//
//           use the generalized clenshaw-curtis method.
//
 l10: hlgth = 0.5*(b-a);
      centr = 0.5*(b+a);
      neval = 25;
      fval[1] = 0.5*f(hlgth+centr);
      fval[13] = f(centr);
      fval[25] = 0.5*f(centr-hlgth);
      for (i=2; i<=12; i++) { //do 20 i=2,12
        u = hlgth*x[i-1];
        isym = 26-i;
        fval[i] = f(u+centr);
        fval[isym] = f(centr-u);
 l20: ;}
//
//           compute the chebyshev series expansion.
//
      qcheb!Real(x.ptr,fval.ptr,cheb12.ptr,cheb24.ptr);
//
//           the modified chebyshev moments are computed by forward
//           recursion, using amom0 and amom1 as starting values.
//
      amom0 = log(fabs((0.1e1-cc)/(0.1e1+cc)));
      amom1 = 0.2e1+cc*amom0;
      res12 = cheb12[1]*amom0+cheb12[2]*amom1;
      res24 = cheb24[1]*amom0+cheb24[2]*amom1;
      for (k=3; k<=13; k++) { //do 30 k=3,13
        amom2 = 0.2e1*cc*amom1-amom0;
        ak22 = (k-2)*(k-2);
        if((k/2)*2 == k) amom2 = amom2-0.4e1/(ak22-0.1e1);
        res12 = res12+cheb12[k]*amom2;
        res24 = res24+cheb24[k]*amom2;
        amom0 = amom1;
        amom1 = amom2;
 l30: ;}
      for (k=14; k<=25; k++) { //do 40 k=14,25
        amom2 = 0.2e1*cc*amom1-amom0;
        ak22 = (k-2)*(k-2);
        if((k/2)*2 == k) amom2 = amom2-0.4e1/(ak22-0.1e1);
        res24 = res24+cheb24[k]*amom2;
        amom0 = amom1;
        amom1 = amom2;
 l40: ;}
      result = res24;
      abserr = fabs(res24-res12);
 l50: return;
}


unittest
{
    alias qc25c!(float, float delegate(float)) fqc25c;
    alias qc25c!(double, double delegate(double)) dqc25c;
    alias qc25c!(double, double function(double)) dfqc25c;
    alias qc25c!(real, real delegate(real)) rqc25c;
}
