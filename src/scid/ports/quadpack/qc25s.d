/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qc25s;


import std.math: fabs, log;

import scid.core.fortran;
import scid.ports.quadpack.qcheb;
import scid.ports.quadpack.qk15w;
import scid.ports.quadpack.qwgts;




///
void qc25s(Real, Func)(Func f, Real a, Real b, Real bl, Real br, Real alfa,
    Real beta, Real* ri_, Real* rj_, Real* rg_, Real* rh_, out Real result,
    out Real abserr, out Real resasc, int integr, out int nev)
{
//***begin prologue  dqc25s
//***date written   810101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a2
//***keywords  25-point clenshaw-curtis integration
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  to compute i = integral of f*w over (bl,br), with error
//            estimate, where the weight function w has a singular
//            behaviour of algebraico-logarithmic type at the points
//            a and/or b. (bl,br) is a part of (a,b).
//***description
//
//        integration rules for integrands having algebraico-logarithmic
//        end point singularities
//        standard fortran subroutine
//        double precision version
//
//        parameters
//           f      - double precision
//                    function subprogram defining the integrand
//                    f(x). the actual name for f needs to be declared
//                    e x t e r n a l  in the driver program.
//
//           a      - double precision
//                    left end point of the original interval
//
//           b      - double precision
//                    right end point of the original interval, b.gt.a
//
//           bl     - double precision
//                    lower limit of integration, bl.ge.a
//
//           br     - double precision
//                    upper limit of integration, br.le.b
//
//           alfa   - double precision
//                    parameter in the weight function
//
//           beta   - double precision
//                    parameter in the weight function
//
//           ri,rj,rg,rh - double precision
//                    modified chebyshev moments for the application
//                    of the generalized clenshaw-curtis
//                    method (computed in subroutine dqmomo)
//
//           result - double precision
//                    approximation to the integral
//                    result is computed by using a generalized
//                    clenshaw-curtis method if b1 = a or br = b.
//                    in all other cases the 15-point kronrod
//                    rule is applied, obtained by optimal addition of
//                    abscissae to the 7-point gauss rule.
//
//           abserr - double precision
//                    estimate of the modulus of the absolute error,
//                    which should equal or exceed abs(i-result)
//
//           resasc - double precision
//                    approximation to the integral of abs(f*w-i/(b-a))
//
//           integr - integer
//                    which determines the weight function
//                    = 1   w(x) = (x-a)**alfa*(b-x)**beta
//                    = 2   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)
//                    = 3   w(x) = (x-a)**alfa*(b-x)**beta*log(b-x)
//                    = 4   w(x) = (x-a)**alfa*(b-x)**beta*log(x-a)*
//                                 log(b-x)
//
//           nev    - integer
//                    number of integrand evaluations
//***references  (none)
//***routines called  dqcheb,dqk15w
//***end prologue  dqc25s
//
      Real centr,dc,factor,fix,hlgth,resabs,res12,
        res24,u;
      Real[13] cheb12_;
      Real[25] cheb24_, fval_;
      int i=1,isym=1;
//
//           the vector x contains the values cos(k*pi/24)
//           k = 1, ..., 11, to be used for the computation of the
//           chebyshev series expansion of f.
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
      auto cheb12 = dimension(cheb12_.ptr, 13);
      auto cheb24 = dimension(cheb24_.ptr, 25);
      auto fval = dimension(fval_.ptr, 25);
      auto rg = dimension(rg_, 25);
      auto rh = dimension(rh_, 25);
      auto ri = dimension(ri_, 25);
      auto rj = dimension(rj_, 25);
      auto x = dimension(x_.ptr, 11);
//
//           list of major variables
//           -----------------------
//
//           fval   - value of the function f at the points
//                    (br-bl)*0.5*cos(k*pi/24)+(br+bl)*0.5
//                    k = 0, ..., 24
//           cheb12 - coefficients of the chebyshev series expansion
//                    of degree 12, for the function f, in the
//                    interval (bl,br)
//           cheb24 - coefficients of the chebyshev series expansion
//                    of degree 24, for the function f, in the
//                    interval (bl,br)
//           res12  - approximation to the integral obtained from cheb12
//           res24  - approximation to the integral obtained from cheb24
//           dqwgts - external function subprogram defining
//                    the four possible weight functions
//           hlgth  - half-length of the interval (bl,br)
//           centr  - mid point of the interval (bl,br)
//
//***first executable statement  dqc25s
      nev = 25;
      if(bl == a && (alfa != 0.0 || integr == 2 || integr == 4))
        goto l10;
      if(br == b && (beta != 0.0 || integr == 3 || integr == 4))
        goto l140;
//
//           if a > bl and b < br, apply the 15-point gauss-kronrod
//           scheme.
//
//
      qk15w!(Real,Func)(f,&qwgts!Real,a,b,alfa,beta,integr,bl,br,
          result,abserr,resabs,resasc);
      nev = 15;
      goto l270;
//
//           this part of the program is executed only if a = bl.
//           ----------------------------------------------------
//
//           compute the chebyshev series expansion of the
//           following function
//           f1 = (0.5*(b+b-br-a)-0.5*(br-a)*x)**beta
//                  *f(0.5*(br-a)*x+0.5*(br+a))
//
 l10: hlgth = 0.5*(br-bl);
      centr = 0.5*(br+bl);
      fix = b-centr;
      fval[1] = 0.5*f(hlgth+centr)*((fix-hlgth)^^beta);
      fval[13] = f(centr)*(fix^^beta);
      fval[25] = 0.5*f(centr-hlgth)*((fix+hlgth)^^beta);
      for (i=2; i<=12; i++) { //do 20 i=2,12
        u = hlgth*x[i-1];
        isym = 26-i;
        fval[i] = f(u+centr)*((fix-u)^^beta);
        fval[isym] = f(centr-u)*((fix+u)^^beta);
 l20: ;}
      factor = hlgth^^(alfa+(cast(Real)0.1e1));
      result = 0.0;
      abserr = 0.0;
      res12 = 0.0;
      res24 = 0.0;
      if(integr > 2) goto l70;
      qcheb!Real(x.ptr,fval.ptr,cheb12.ptr,cheb24.ptr);
//
//           integr = 1  (or 2)
//
      for (i=1; i<=13; i++) { //do 30 i=1,13
        res12 = res12+cheb12[i]*ri[i];
        res24 = res24+cheb24[i]*ri[i];
 l30: ;}
      for (i=14; i<=25; i++) { //do 40 i=14,25
        res24 = res24+cheb24[i]*ri[i];
 l40: ;}
      if(integr == 1) goto l130;
//
//           integr = 2
//
      dc = log(br-bl);
      result = res24*dc;
      abserr = fabs((res24-res12)*dc);
      res12 = 0.0;
      res24 = 0.0;
      for (i=1; i<=13; i++) { //do 50 i=1,13
        res12 = res12+cheb12[i]*rg[i];
        res24 = res12+cheb24[i]*rg[i];
 l50: ;}
      for (i=14; i<=25; i++) { //do 60 i=14,25
        res24 = res24+cheb24[i]*rg[i];
 l60: ;}
      goto l130;
//
//           compute the chebyshev series expansion of the
//           following function
//           f4 = f1*log(0.5*(b+b-br-a)-0.5*(br-a)*x)
//
 l70: fval[1] = fval[1]*log(fix-hlgth);
      fval[13] = fval[13]*log(fix);
      fval[25] = fval[25]*log(fix+hlgth);
      for (i=2; i<=12; i++) { //do 80 i=2,12
        u = hlgth*x[i-1];
        isym = 26-i;
        fval[i] = fval[i]*log(fix-u);
        fval[isym] = fval[isym]*log(fix+u);
 l80: ;}
      qcheb!Real(x.ptr,fval.ptr,cheb12.ptr,cheb24.ptr);
//
//           integr = 3  (or 4)
//
      for (i=1; i<=13; i++) { //do 90 i=1,13
        res12 = res12+cheb12[i]*ri[i];
        res24 = res24+cheb24[i]*ri[i];
 l90: ;}
      for (i=14; i<=25; i++) { //do 100 i=14,25
        res24 = res24+cheb24[i]*ri[i];
l100: ;}
      if(integr == 3) goto l130;
//
//           integr = 4
//
      dc = log(br-bl);
      result = res24*dc;
      abserr = fabs((res24-res12)*dc);
      res12 = 0.0;
      res24 = 0.0;
      for (i=1; i<=13; i++) { //do 110 i=1,13
        res12 = res12+cheb12[i]*rg[i];
        res24 = res24+cheb24[i]*rg[i];
l110: ;}
      for (i=14; i<=25; i++) { //do 120 i=14,25
        res24 = res24+cheb24[i]*rg[i];
l120: ;}
l130: result = (result+res24)*factor;
      abserr = (abserr+fabs(res24-res12))*factor;
      goto l270;
//
//           this part of the program is executed only if b = br.
//           ----------------------------------------------------
//
//           compute the chebyshev series expansion of the
//           following function
//           f2 = (0.5*(b+bl-a-a)+0.5*(b-bl)*x)**alfa
//                *f(0.5*(b-bl)*x+0.5*(b+bl))
//
l140: hlgth = 0.5*(br-bl);
      centr = 0.5*(br+bl);
      fix = centr-a;
      fval[1] = 0.5*f(hlgth+centr)*((fix+hlgth)^^alfa);
      fval[13] = f(centr)*(fix^^alfa);
      fval[25] = 0.5*f(centr-hlgth)*((fix-hlgth)^^alfa);
      for (i=2; i<=12; i++) { //do 150 i=2,12
        u = hlgth*x[i-1];
        isym = 26-i;
        fval[i] = f(u+centr)*((fix+u)^^alfa);
        fval[isym] = f(centr-u)*((fix-u)^^alfa);
l150: ;}
      factor = hlgth^^(beta+(cast(Real)0.1e1));
      result = 0.0;
      abserr = 0.0;
      res12 = 0.0;
      res24 = 0.0;
      if(integr == 2 || integr == 4) goto l200;
//
//           integr = 1  (or 3)
//
      qcheb!Real(x.ptr,fval.ptr,cheb12.ptr,cheb24.ptr);
      for (i=1; i<=13; i++) { //do 160 i=1,13
        res12 = res12+cheb12[i]*rj[i];
        res24 = res24+cheb24[i]*rj[i];
l160: ;}
      for (i=14; i<=25; i++) { //do 170 i=14,25
        res24 = res24+cheb24[i]*rj[i];
l170: ;}
      if(integr == 1) goto l260;
//
//           integr = 3
//
      dc = log(br-bl);
      result = res24*dc;
      abserr = fabs((res24-res12)*dc);
      res12 = 0.0;
      res24 = 0.0;
      for (i=1; i<=13; i++) { //do 180 i=1,13
        res12 = res12+cheb12[i]*rh[i];
        res24 = res24+cheb24[i]*rh[i];
l180: ;}
      for (i=14; i<=25; i++) { //do 190 i=14,25
        res24 = res24+cheb24[i]*rh[i];
l190: ;}
      goto l260;
//
//           compute the chebyshev series expansion of the
//           following function
//           f3 = f2*log(0.5*(b-bl)*x+0.5*(b+bl-a-a))
//
l200: fval[1] = fval[1]*log(hlgth+fix);
      fval[13] = fval[13]*log(fix);
      fval[25] = fval[25]*log(fix-hlgth);
      for (i=2; i<=12; i++) { //do 210 i=2,12
        u = hlgth*x[i-1];
        isym = 26-i;
        fval[i] = fval[i]*log(u+fix);
        fval[isym] = fval[isym]*log(fix-u);
l210: ;}
      qcheb!Real(x.ptr,fval.ptr,cheb12.ptr,cheb24.ptr);
//
//           integr = 2  (or 4)
//
      for (i=1; i<=13; i++) { //do 220 i=1,13
        res12 = res12+cheb12[i]*rj[i];
        res24 = res24+cheb24[i]*rj[i];
l220: ;}
      for (i=14; i<=25; i++) { //do 230 i=14,25
        res24 = res24+cheb24[i]*rj[i];
l230: ;}
      if(integr == 2) goto l260;
      dc = log(br-bl);
      result = res24*dc;
      abserr = fabs((res24-res12)*dc);
      res12 = 0.0;
      res24 = 0.0;
//
//           integr = 4
//
      for (i=1; i<=13; i++) { //do 240 i=1,13
        res12 = res12+cheb12[i]*rh[i];
        res24 = res24+cheb24[i]*rh[i];
l240: ;}
      for (i=14; i<=25; i++) { //do 250 i=14,25
        res24 = res24+cheb24[i]*rh[i];
l250: ;}
l260: result = (result+res24)*factor;
      abserr = (abserr+fabs(res24-res12))*factor;
l270: return;
}


unittest
{
    alias qc25s!(float, float delegate(float)) fqc25s;
    alias qc25s!(double, double delegate(double)) dqc25s;
    alias qc25s!(double, double function(double)) dfqc25s;
    alias qc25s!(real, real delegate(real)) rqc25s;
}
