// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qc25f;


import std.math;

import scid.core.fortran;
import scid.ports.linpack.gtsl;
import scid.ports.quadpack.qk15w;
import scid.ports.quadpack.qwgtf;
import scid.ports.quadpack.qcheb;

version(unittest) import scid.core.testing;




///
void qc25f(Real, Func)(Func f, Real a, Real b, Real omega, int integr,
    int nrmom, int maxp1, int ksave, out Real result,
    out Real abserr, out int neval, out Real resabs, out Real resasc,
    out int momcom, Real* chebmo_)
{
//***begin prologue  dqc25f
//***date written   810101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a2
//***keywords  integration rules for functions with cos or sin
//             factor, clenshaw-curtis, gauss-kronrod
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  to compute the integral i=integral of f(x) over (a,b)
//            where w(x) = cos(omega*x) or w(x)=sin(omega*x) and to
//            compute j = integral of abs(f) over (a,b). for small value
//            of omega or small intervals (a,b) the 15-point gauss-kronro
//            rule is used. otherwise a generalized clenshaw-curtis
//            method is used.
//***description
//
//        integration rules for functions with cos or sin factor
//        standard fortran subroutine
//        double precision version
//
//        parameters
//         on entry
//           f      - double precision
//                    function subprogram defining the integrand
//                    function f(x). the actual name for f needs to
//                    be declared e x t e r n a l in the calling program.
//
//           a      - double precision
//                    lower limit of integration
//
//           b      - double precision
//                    upper limit of integration
//
//           omega  - double precision
//                    parameter in the weight function
//
//           integr - integer
//                    indicates which weight function is to be used
//                       integr = 1   w(x) = cos(omega*x)
//                       integr = 2   w(x) = sin(omega*x)
//
//           nrmom  - integer
//                    the length of interval (a,b) is equal to the length
//                    of the original integration interval divided by
//                    2**nrmom (we suppose that the routine is used in an
//                    adaptive integration process, otherwise set
//                    nrmom = 0). nrmom must be zero at the first call.
//
//           maxp1  - integer
//                    gives an upper bound on the number of chebyshev
//                    moments which can be stored, i.e. for the
//                    intervals of lengths abs(bb-aa)*2**(-l),
//                    l = 0,1,2, ..., maxp1-2.
//
//           ksave  - integer
//                    key which is one when the moments for the
//                    current interval have been computed
//
//         on return
//           result - double precision
//                    approximation to the integral i
//
//           abserr - double precision
//                    estimate of the modulus of the absolute
//                    error, which should equal or exceed abs(i-result)
//
//           neval  - integer
//                    number of integrand evaluations
//
//           resabs - double precision
//                    approximation to the integral j
//
//           resasc - double precision
//                    approximation to the integral of abs(f-i/(b-a))
//
//         on entry and return
//           momcom - integer
//                    for each interval length we need to compute the
//                    chebyshev moments. momcom counts the number of
//                    intervals for which these moments have already been
//                    computed. if nrmom.lt.momcom or ksave = 1, the
//                    chebyshev moments for the interval (a,b) have
//                    already been computed and stored, otherwise we
//                    compute them and we increase momcom.
//
//           chebmo - double precision
//                    array of dimension at least (maxp1,25) containing
//                    the modified chebyshev moments for the first momcom
//                    momcom interval lengths
//
// ......................................................................
//***references  (none)
//***routines called  d1mach,dgtsl,dqcheb,dqk15w,dqwgtf
//***end prologue  dqc25f
//
      Real ac,an,an2,as,asap,ass,centr,
       conc,cons,cospar,
       estc,ests,hlgth,oflow,parint,par2,par22,
       p2,p3,p4,resc12,resc24,ress12,ress24,
       sinpar;
      int i=1,iers=1,isym=1,j=1,k=1,m=1,
       noequ=1,noeq1=1;
//
      Real[13] cheb12_;
      Real[25] cheb24_, d_, d1_, d2_, fval_;
      Real[28] v_;
      auto chebmo = dimension(chebmo_, maxp1, 25);
      auto cheb12 = dimension(cheb12_.ptr, 13);
      auto cheb24 = dimension(cheb24_.ptr, 25);
      auto d = dimension(d_.ptr, 25);
      auto d1 = dimension(d1_.ptr, 25);
      auto d2 = dimension(d2_.ptr, 25);
      auto fval = dimension(fval_.ptr, 25);
      auto v = dimension(v_.ptr, 28);
//
//           the vector x contains the values cos(k*pi/24)
//           k = 1, ...,11, to be used for the chebyshev expansion of f
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
        0.1305261922_2005159154_8406227895_489
      ];
      auto x = dimension(x_.ptr, 11);
//
//           list of major variables
//           -----------------------
//
//           centr  - mid point of the integration interval
//           hlgth  - half-length of the integration interval
//           fval   - value of the function f at the points
//                    (b-a)*0.5*cos(k*pi/12) + (b+a)*0.5, k = 0, ..., 24
//           cheb12 - coefficients of the chebyshev series expansion
//                    of degree 12, for the function f, in the
//                    interval (a,b)
//           cheb24 - coefficients of the chebyshev series expansion
//                    of degree 24, for the function f, in the
//                    interval (a,b)
//           resc12 - approximation to the integral of
//                    cos(0.5*(b-a)*omega*x)*f(0.5*(b-a)*x+0.5*(b+a))
//                    over (-1,+1), using the chebyshev series
//                    expansion of degree 12
//           resc24 - approximation to the same integral, using the
//                    chebyshev series expansion of degree 24
//           ress12 - the analogue of resc12 for the sine
//           ress24 - the analogue of resc24 for the sine
//
//
//           machine dependent constant
//           --------------------------
//
//           oflow is the largest positive magnitude.
//
//***first executable statement  dqc25f
      oflow = Real.max;
//
      centr = 0.5*(b+a);
      hlgth = 0.5*(b-a);
      parint = omega*hlgth;
//
//           compute the integral using the 15-point gauss-kronrod
//           formula if the value of the parameter in the integrand
//           is small.
//
      if (fabs(parint) > 0.2e1) goto l10;
      qk15w!(Real, Func)(f,&qwgtf!Real,omega,p2,p3,p4,integr,a,b,result,
        abserr,resabs,resasc);
      neval = 15;
      goto l170;
//
//           compute the integral using the generalized clenshaw-
//           curtis method.
//
 l10: conc = hlgth*cos(centr*omega);
      cons = hlgth*sin(centr*omega);
      resasc = oflow;
      neval = 25;
//
//           check whether the chebyshev moments for this interval
//           have already been computed.
//
      if (nrmom < momcom || ksave == 1) goto l120;
//
//           compute a new set of chebyshev moments.
//
      m = momcom+1;
      par2 = parint*parint;
      par22 = par2+0.2e1;
      sinpar = sin(parint);
      cospar = cos(parint);
//
//           compute the chebyshev moments with respect to cosine.
//
      v[1] = 0.2e1*sinpar/parint;
      v[2] = (0.8e1*cospar+(par2+par2-0.8e1)*sinpar/parint)/par2;
      v[3] = (0.32e2*(par2-0.12e2)*cospar+(0.2e1*
       ((par2-0.80e2)*par2+0.192e3)*sinpar)/parint)/(par2*par2);
      ac = 0.8e1*cospar;
      as = 0.24e2*parint*sinpar;
      if (fabs(parint) > 0.24e2) goto l30;
//
//           compute the chebyshev moments as the solutions of a
//           boundary value problem with 1 initial value (v(3)) and 1
//           end value (computed using an asymptotic formula).
//
      noequ = 25;
      noeq1 = noequ-1;
      an = 0.6e1;
      for (k=1; k<=noeq1; k++) { // end: l20
        an2 = an*an;
        d[k] = -0.2e1*(an2-0.4e1)*(par22-an2-an2);
        d2[k] = (an-0.1e1)*(an-0.2e1)*par2;
        d1[k+1] = (an+0.3e1)*(an+0.4e1)*par2;
        v[k+3] = as-(an2-0.4e1)*ac;
        an = an+0.2e1;
 l20:;}
      an2 = an*an;
      d[noequ] = -0.2e1*(an2-0.4e1)*(par22-an2-an2);
      v[noequ+3] = as-(an2-0.4e1)*ac;
      v[4] = v[4]-0.56e2*par2*v[3];
      ass = parint*sinpar;
      asap = (((((0.210e3*par2-0.1e1)*cospar-(0.105e3*par2
        -0.63e2)*ass)/an2-(0.1e1-0.15e2*par2)*cospar
        +0.15e2*ass)/an2-cospar+0.3e1*ass)/an2-cospar)/an2;
      v[noequ+3] = v[noequ+3]-0.2e1*asap*par2*(an-0.1e1)*
         (an-0.2e1);
//
//           solve the tridiagonal system by means of gaussian
//           elimination with partial pivoting.
//
//***        call to dgtsl must be replaced by call to
//***        double precision version of linpack routine sgtsl
//
      gtsl!Real(noequ,d1.ptr,d.ptr,d2.ptr,v.ptr+3,iers);
      goto l50;
//
//           compute the chebyshev moments by means of forward
//           recursion.
//
 l30: an = 0.4e1;
      for (i=4; i<=13; i++) { // end: l40
        an2 = an*an;
        v[i] = ((an2-0.4e1)*(0.2e1*(par22-an2-an2)*v[i-1]-ac)
          +as-par2*(an+0.1e1)*(an+0.2e1)*v[i-2])/
          (par2*(an-0.1e1)*(an-0.2e1));
        an = an+0.2e1;
 l40:;}
 l50: for (j=1; j<=13; j++) { // end: l60
        chebmo[m,2*j-1] = v[j];
 l60:;}
//
//           compute the chebyshev moments with respect to sine.
//
      v[1] = 0.2e1*(sinpar-parint*cospar)/par2;
      v[2] = (0.18e2-0.48e2/par2)*sinpar/par2
        +(-0.2e1+0.48e2/par2)*cospar/parint;
      ac = -0.24e2*parint*cospar;
      as = -0.8e1*sinpar;
      if (fabs(parint) > 0.24e2) goto l80;
//
//           compute the chebyshev moments as the solutions of a boundary
//           value problem with 1 initial value (v(2)) and 1 end value
//           (computed using an asymptotic formula).
//
      an = 0.5e1;
      for (k=1 ; k<=noeq1; k++) { // end: l70
        an2 = an*an;
        d[k] = -0.2e1*(an2-0.4e1)*(par22-an2-an2);
        d2[k] = (an-0.1e1)*(an-0.2e1)*par2;
        d1[k+1] = (an+0.3e1)*(an+0.4e1)*par2;
        v[k+2] = ac+(an2-0.4e1)*as;
        an = an+0.2e1;
 l70:;}
      an2 = an*an;
      d[noequ] = -0.2e1*(an2-0.4e1)*(par22-an2-an2);
      v[noequ+2] = ac+(an2-0.4e1)*as;
      v[3] = v[3]-0.42e2*par2*v[2];
      ass = parint*cospar;
      asap = (((((0.105e3*par2-0.63e2)*ass+(0.210e3*par2
        -0.1e1)*sinpar)/an2+(0.15e2*par2-0.1e1)*sinpar-
        0.15e2*ass)/an2-0.3e1*ass-sinpar)/an2-sinpar)/an2;
      v[noequ+2] = v[noequ+2]-0.2e1*asap*par2*(an-0.1e1)
        *(an-0.2e1);
//
//           solve the tridiagonal system by means of gaussian
//           elimination with partial pivoting.
//
//***        call to dgtsl must be replaced by call to
//***        double precision version of linpack routine sgtsl
//
      gtsl!Real(noequ,d1.ptr,d.ptr,d2.ptr,v.ptr+2,iers);
      goto l100;
//
//           compute the chebyshev moments by means of forward recursion.
//
 l80: an = 0.3e1;
      for (i=3; i<=12; i++) { // end: l90
        an2 = an*an;
        v[i] = ((an2-0.4e1)*(0.2e1*(par22-an2-an2)*v[i-1]+as)
          +ac-par2*(an+0.1e1)*(an+0.2e1)*v[i-2])
          /(par2*(an-0.1e1)*(an-0.2e1));
        an = an+0.2e1;
 l90:;}  
l100: for (j=1; j<=12; j++) { // end: l110
        chebmo[m,2*j] = v[j];
l110:;}
l120: if (nrmom < momcom) m = nrmom+1;
       if (momcom < (maxp1-1) && nrmom >= momcom) momcom = momcom+1;
//
//           compute the coefficients of the chebyshev expansions
//           of degrees 12 and 24 of the function f.
//
      fval[1] = 0.5*f(centr+hlgth);
      fval[13] = f(centr);
      fval[25] = 0.5*f(centr-hlgth);
      for(i=2; i<=12; i++) { // end: l130
        isym = 26-i;
        fval[i] = f(hlgth*x[i-1]+centr);
        fval[isym] = f(centr-hlgth*x[i-1]);
l130:;}
      qcheb!Real(x.ptr,fval.ptr,cheb12.ptr,cheb24.ptr);
//
//           compute the integral and error estimates.
//
      resc12 = cheb12[13]*chebmo[m,13];
      ress12 = 0.0;
      k = 11;
      for (j=1; j<=6; j++) { // end: l140
        resc12 = resc12+cheb12[k]*chebmo[m,k];
        ress12 = ress12+cheb12[k+1]*chebmo[m,k+1];
        k = k-2;
l140:;}
      resc24 = cheb24[25]*chebmo[m,25];
      ress24 = 0.0;
      resabs = fabs(cheb24[25]);
      k = 23;
      for (j=1; j<=12; j++) { // end: l150
        resc24 = resc24+cheb24[k]*chebmo[m,k];
        ress24 = ress24+cheb24[k+1]*chebmo[m,k+1];
        resabs = fabs(cheb24[k])+fabs(cheb24[k+1]);
        k = k-2;
l150:;}
      estc = fabs(resc24-resc12);
      ests = fabs(ress24-ress12);
      resabs = resabs*fabs(hlgth);
      if(integr == 2) goto l160;
      result = conc*resc24-cons*ress24;
      abserr = fabs(conc*estc)+fabs(cons*ests);
      goto l170;
l160: result = conc*ress24+cons*resc24;
      abserr = fabs(conc*ests)+fabs(cons*estc);
l170: return;
}


unittest
{
    alias qc25f!(float, float delegate(float)) fqc25f;
    alias qc25f!(double, double delegate(double)) dqc25f;
    alias qc25f!(double, double function(double)) dfqc25f;
    alias qc25f!(real, real delegate(real)) rqc25f;
}


unittest
{
    double f(double x) { return x*x; }

    double a = 0.0, b = 2.0;
    int integr = 1, nrmom = 0, maxp1 = 21, ksave = 0, momcom = 0;
    double result, abserr, resabs, resasc;
    int neval;
    double[] chebmo = new double[maxp1*25];

    // Compute using 15-point Gauss-Kronrod formula.
    double omega = 1.0;
    double ans = 2*(2*cos(2.0) + sin(2.0));
    qc25f(&f, a, b, omega, integr, nrmom, maxp1, ksave, result, abserr,
        neval, resabs, resasc, momcom, chebmo.ptr);
    assert (isAccurate(result, abserr, ans, 1e-8));


    // Compute using Clenshaw-Curtis method.
    omega = 3.0;
    ans = 2*(6*cos(6.0) + 17*sin(6.0))/27;
    qc25f(&f, a, b, omega, integr, nrmom, maxp1, ksave, result, abserr,
        neval, resabs, resasc, momcom, chebmo.ptr);
    assert (isAccurate(result, 1e-10, ans, 1e-8));

    // Compute again, using the stored moments.
    ksave = 1;
    qc25f(&f, a, b, omega, integr, nrmom, maxp1, ksave, result, abserr,
        neval, resabs, resasc, momcom, chebmo.ptr);
    assert (isAccurate(result, 1e-10, ans, 1e-8));
}


unittest
{
    double alpha = 3.0;
    double f(double x) { return x <= 0.0 ? 0.0 : exp(-(2.0^^(-alpha))*x)/sqrt(x); }

    double a = 1.0, omega = 1.0;
    int integr = 1, nrmom = 0, maxp1 = 21, ksave = 0, momcom = 0;
    double result, abserr, resabs, resasc;
    int neval;
    double[] chebmo = new double[maxp1*25];

    // Compute using 15-point Gauss-Kronrod formula.
    double b = 2.0;
    double ans = 0.0729155596656107239492740504503422453;
    qc25f(&f, a, b, omega, integr, nrmom, maxp1, ksave, result, abserr,
        neval, resabs, resasc, momcom, chebmo.ptr);
    assert (approxEqual(result, ans, 1e-8));


    // Compute using Clenshaw-Curtis method.
    b = 6.0;
    ans = -0.50789465569969289260403404655063469;
    qc25f(&f, a, b, omega, integr, nrmom, maxp1, ksave, result, abserr,
        neval, resabs, resasc, momcom, chebmo.ptr);
    assert (approxEqual(result, ans, 1e-8));

    // Compute again, using the stored moments.
    ksave = 1;
    qc25f(&f, a, b, omega, integr, nrmom, maxp1, ksave, result, abserr,
        neval, resabs, resasc, momcom, chebmo.ptr);
    assert (approxEqual(result, ans, 1e-8));
}

