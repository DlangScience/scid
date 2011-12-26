// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qawo;


import std.conv;
import scid.ports.quadpack.qawoe;

version(unittest)
{
    import std.math;
    import scid.core.testing;
}




///
void qawo(Real, Func)(Func f, Real a, Real b, Real omega, int integr,
    Real epsabs, Real epsrel, out Real result, out Real abserr,
    out int neval, out int ier, int leniw, int maxp1, int lenw,
    out int last, int* iwork, Real* work)
{
//***begin prologue  qawo
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a1
//***keywords  automatic integrator, special-purpose,
//             integrand with oscillatory cos or sin factor,
//             clenshaw-curtis method, (end point) singularities,
//             extrapolation, globally adaptive
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            definite integral
//            i = integral of f(x)*w(x) over (a,b)
//            where w(x) = cos(omega*x) or w(x) = sin(omega*x),
//            hopefully satisfying following claim for accuracy
//            abs(i-result).le.max(epsabs,epsrel*abs(i)).
//***description
//
//        computation of oscillatory integrals
//        standard fortran subroutine
//        real version
//
//        parameters
//         on entry
//            f      - real
//                     function subprogram defining the function
//                     f(x).  the actual name for f needs to be
//                     declared e x t e r n a l in the driver program.
//
//            a      - real
//                     lower limit of integration
//
//            b      - real
//                     upper limit of integration
//
//            omega  - real
//                     parameter in the integrand weight function
//
//            integr - integer
//                     indicates which of the weight functions is used
//                     integr = 1      w(x) = cos(omega*x)
//                     integr = 2      w(x) = sin(omega*x)
//                     if integr.ne.1.and.integr.ne.2, the routine will
//                     end with ier = 6.
//
//            epsabs - real
//                     absolute accuracy requested
//            epsrel - real
//                     relative accuracy requested
//                     if epsabs.le.0 and
//                     epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                     the routine will end with ier = 6.
//
//         on return
//            result - real
//                     approximation to the integral
//
//            abserr - real
//                     estimate of the modulus of the absolute error,
//                     which should equal or exceed abs(i-result)
//
//            neval  - integer
//                     number of  integrand evaluations
//
//            ier    - integer
//                     ier = 0 normal and reliable termination of the
//                             routine. it is assumed that the requested
//                             accuracy has been achieved.
//                   - ier.gt.0 abnormal termination of the routine.
//                             the estimates for integral and error are
//                             less reliable. it is assumed that the
//                             requested accuracy has not been achieved.
//            error messages
//                     ier = 1 maximum number of subdivisions allowed
//                             (= leniw/2) has been achieved. one can
//                             allow more subdivisions by increasing the
//                             value of leniw (and taking the according
//                             dimension adjustments into account).
//                             however, if this yields no improvement it
//                             is advised to analyze the integrand in
//                             order to determine the integration
//                             difficulties. if the position of a local
//                             difficulty can be determined (e.g.
//                             singularity, discontinuity within the
//                             interval) one will probably gain from
//                             splitting up the interval at this point
//                             and calling the integrator on the
//                             subranges. if possible, an appropriate
//                             special-purpose integrator should be used
//                             which is designed for handling the type of
//                             difficulty involved.
//                         = 2 the occurrence of roundoff error is
//                             detected, which prevents the requested
//                             tolerance from being achieved.
//                             the error may be under-estimated.
//                         = 3 extremely bad integrand behaviour occurs
//                             at some interior points of the
//                             integration interval.
//                         = 4 the algorithm does not converge.
//                             roundoff error is detected in the
//                             extrapolation table. it is presumed that
//                             the requested tolerance cannot be achieved
//                             due to roundoff in the extrapolation
//                             table, and that the returned result is
//                             the best which can be obtained.
//                         = 5 the integral is probably divergent, or
//                             slowly convergent. it must be noted that
//                             divergence can occur with any other value
//                             of ier.
//                         = 6 the input is invalid, because
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
//                             or (integr.ne.1 and integr.ne.2),
//                             or leniw.lt.2 or maxp1.lt.1 or
//                             lenw.lt.leniw*2+maxp1*25.
//                             result, abserr, neval, last are set to
//                             zero. except when leniw, maxp1 or lenw are
//                             invalid, work(limit*2+1), work(limit*3+1),
//                             iwork(1), iwork(limit+1) are set to zero,
//                             work(1) is set to a and work(limit+1) to
//                             b.
//
//         dimensioning parameters
//            leniw  - integer
//                     dimensioning parameter for iwork.
//                     leniw/2 equals the maximum number of subintervals
//                     allowed in the partition of the given integration
//                     interval (a,b), leniw.ge.2.
//                     if leniw.lt.2, the routine will end with ier = 6.
//
//            maxp1  - integer
//                     gives an upper bound on the number of chebyshev
//                     moments which can be stored, i.e. for the
//                     intervals of lengths abs(b-a)*2**(-l),
//                     l=0,1, ..., maxp1-2, maxp1.ge.1
//                     if maxp1.lt.1, the routine will end with ier = 6.
//
//            lenw   - integer
//                     dimensioning parameter for work
//                     lenw must be at least leniw*2+maxp1*25.
//                     if lenw.lt.(leniw*2+maxp1*25), the routine will
//                     end with ier = 6.
//
//            last   - integer
//                     on return, last equals the number of subintervals
//                     produced in the subdivision process, which
//                     determines the number of significant elements
//                     actually in the work arrays.
//
//         work arrays
//            iwork  - integer
//                     vector of dimension at least leniw
//                     on return, the first k elements of which contain
//                     pointers to the error estimates over the
//                     subintervals, such that work(limit*3+iwork(1)), ..
//                     work(limit*3+iwork(k)) form a decreasing
//                     sequence, with limit = lenw/2 , and k = last
//                     if last.le.(limit/2+2), and k = limit+1-last
//                     otherwise.
//                     furthermore, iwork(limit+1), ..., iwork(limit+
//                     last) indicate the subdivision levels of the
//                     subintervals, such that iwork(limit+i) = l means
//                     that the subinterval numbered i is of length
//                     abs(b-a)*2**(1-l).
//
//            work   - real
//                     vector of dimension at least lenw
//                     on return
//                     work(1), ..., work(last) contain the left
//                      end points of the subintervals in the
//                      partition of (a,b),
//                     work(limit+1), ..., work(limit+last) contain
//                      the right end points,
//                     work(limit*2+1), ..., work(limit*2+last) contain
//                      the integral approximations over the
//                      subintervals,
//                     work(limit*3+1), ..., work(limit*3+last)
//                      contain the error estimates.
//                     work(limit*4+1), ..., work(limit*4+maxp1*25)
//                      provide space for storing the chebyshev moments.
//                     note that limit = lenw/2.
//
//***references  (none)
//***routines called  qawoe,xerror
//***end prologue  qawo
//
       int lvl=1,l1=1,l2=1,l3=1,l4=1,limit=1,momcom=1;
//
       //dimension iwork(leniw),work(lenw)
//
//         check validity of leniw, maxp1 and lenw.
//
//***first executable statement  qawo
      ier = 6;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      if(leniw < 2 || maxp1 < 1 || lenw < (leniw*2+maxp1*25))
        goto l10;
//
//         prepare call for qawoe
//
      limit = leniw/2;
      l1 = limit;
      l2 = limit+l1;
      l3 = limit+l2;
      l4 = limit+l3;
      qawoe!(Real,Func)(f,a,b,omega,integr,epsabs,epsrel,limit,1,maxp1,result,
        abserr,neval,ier,last,work,work+l1,work+l2,work+l3,
        iwork,iwork+l1,momcom,work+l4);
//
//         call error handler if necessary
//
      lvl = 0;
l10:  if(ier == 6) lvl = 1;
      if(ier != 0)
        throw new Exception("abnormal return from  qawo: "~to!string(ier));
      return;
}


unittest
{
    alias qawo!(float, float delegate(float)) fqawo;
    alias qawo!(double, double delegate(double)) dqawo;
    alias qawo!(double, double function(double)) dfqawo;
    alias qawo!(real, real delegate(real)) rqawo;
}


unittest
{
    // From the QUADPACK book
    double dlog(double x) { return x > 0.0 ? log(x) : 0.0; }

    double a = 0.0, b = 1.0, omega = 10*PI;
    int integr = 2;
    double epsabs = 0.0, epsrel = 1e-8;
    double result, abserr;
    int neval, ier, last;
    const int leniw = 1000, maxp1=21, lenw=leniw*2+maxp1*25;

    int[leniw] iwork;
    double[lenw] work;

    qawo(&dlog, a, b, omega, integr, epsabs, epsrel, result, abserr,
        neval, ier, leniw, maxp1, lenw, last, iwork.ptr, work.ptr);

    double ans = -0.128136848399167;
    check (isAccurate!double(result, abserr, ans, epsrel, epsabs));
}


unittest
{
    // This is integral 14 in the QUADPACK book.
    real alpha = 3.0;
    real f(real x) { return x <= 0.0 ? 0.0 : exp(-(2.0L^^(-alpha))*x)/sqrt(x); }
    
    real a = 0.0;
    real b = 20*(2.0L^^alpha);
    real omega = 1.0;
    int integr = 1; // cos(x)
    real epsabs = 0.0;
    real epsrel = 1e-8;

    real result, abserr;
    int neval, ier, last;

    enum
    {
        leniw = 1000,
        maxp1 = 21,
        lenw = leniw*2 + maxp1*25
    }

    int[leniw] iwork;
    real[lenw] work;

    qawo(&f, a, b, omega, integr, epsabs, epsrel, result, abserr,
        neval, ier, leniw, maxp1, lenw, last, iwork.ptr, work.ptr);
    
    real eps = 1e-20;
    real ans = (1 + eps)*sqrt(PI) * ((1+(4.0L^^(-alpha)))^^(-0.25L))
        * cos(atan(2.0L^^alpha)/2);
    check (isAccurate(result, abserr, ans, epsrel, epsabs));
}

