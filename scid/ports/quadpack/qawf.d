// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qawf;


import std.conv;
import scid.ports.quadpack.qawfe;

version(unittest)
{
    import std.math;
    import scid.core.testing;
}



///
void qawf(Real, Func)(Func f, Real a, Real omega, int integr, Real epsabs,
    out Real result, out Real abserr, out int neval, out int ier,
    int limlst, out int lst, int leniw, int maxp1, int lenw,
    int* iwork, Real* work)
{
//***begin prologue  dqawf
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a3a1
//***keywords  automatic integrator, special-purpose,fourier
//             integral, integration between zeros with dqawoe,
//             convergence acceleration with dqelg
//***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            fourier integral i=integral of f(x)*w(x) over (a,infinity)
//            where w(x) = cos(omega*x) or w(x) = sin(omega*x).
//            hopefully satisfying following claim for accuracy
//            abs(i-result).le.epsabs.
//***description
//
//        computation of fourier integrals
//        standard fortran subroutine
//        double precision version
//
//
//        parameters
//         on entry
//            f      - double precision
//                     function subprogram defining the integrand
//                     function f(x). the actual name for f needs to be
//                     declared e x t e r n a l in the driver program.
//
//            a      - double precision
//                     lower limit of integration
//
//            omega  - double precision
//                     parameter in the integrand weight function
//
//            integr - integer
//                     indicates which of the weight functions is used
//                     integr = 1      w(x) = cos(omega*x)
//                     integr = 2      w(x) = sin(omega*x)
//                     if integr.ne.1.and.integr.ne.2, the routine
//                     will end with ier = 6.
//
//            epsabs - double precision
//                     absolute accuracy requested, epsabs.gt.0.
//                     if epsabs.le.0, the routine will end with ier = 6.
//
//         on return
//            result - double precision
//                     approximation to the integral
//
//            abserr - double precision
//                     estimate of the modulus of the absolute error,
//                     which should equal or exceed abs(i-result)
//
//            neval  - integer
//                     number of integrand evaluations
//
//            ier    - integer
//                     ier = 0 normal and reliable termination of the
//                             routine. it is assumed that the requested
//                             accuracy has been achieved.
//                     ier.gt.0 abnormal termination of the routine.
//                             the estimates for integral and error are
//                             less reliable. it is assumed that the
//                             requested accuracy has not been achieved.
//            error messages
//                    if omega.ne.0
//                     ier = 1 maximum number of cycles allowed
//                             has been achieved, i.e. of subintervals
//                             (a+(k-1)c,a+kc) where
//                             c = (2*int(abs(omega))+1)*pi/abs(omega),
//                             for k = 1, 2, ..., lst.
//                             one can allow more cycles by increasing
//                             the value of limlst (and taking the
//                             according dimension adjustments into
//                             account). examine the array iwork which
//                             contains the error flags on the cycles, in
//                             order to look for eventual local
//                             integration difficulties.
//                             if the position of a local difficulty
//                             can be determined (e.g. singularity,
//                             discontinuity within the interval) one
//                             will probably gain from splitting up the
//                             interval at this point and calling
//                             appropriate integrators on the subranges.
//                         = 4 the extrapolation table constructed for
//                             convergence accelaration of the series
//                             formed by the integral contributions over
//                             the cycles, does not converge to within
//                             the requested accuracy.
//                             as in the case of ier = 1, it is advised
//                             to examine the array iwork which contains
//                             the error flags on the cycles.
//                         = 6 the input is invalid because
//                             (integr.ne.1 and integr.ne.2) or
//                              epsabs.le.0 or limlst.lt.1 or
//                              leniw.lt.(limlst+2) or maxp1.lt.1 or
//                              lenw.lt.(leniw*2+maxp1*25).
//                              result, abserr, neval, lst are set to
//                              zero.
//                         = 7 bad integrand behaviour occurs within
//                             one or more of the cycles. location and
//                             type of the difficulty involved can be
//                             determined from the first lst elements of
//                             vector iwork.  here lst is the number of
//                             cycles actually needed (see below).
//                             iwork(k) = 1 the maximum number of
//                                          subdivisions (=(leniw-limlst)
//                                          /2) has been achieved on the
//                                          k th cycle.
//                                      = 2 occurrence of roundoff error
//                                          is detected and prevents the
//                                          tolerance imposed on the k th
//                                          cycle, from being achieved
//                                          on this cycle.
//                                      = 3 extremely bad integrand
//                                          behaviour occurs at some
//                                          points of the k th cycle.
//                                      = 4 the integration procedure
//                                          over the k th cycle does
//                                          not converge (to within the
//                                          required accuracy) due to
//                                          roundoff in the extrapolation
//                                          procedure invoked on this
//                                          cycle. it is assumed that the
//                                          result on this interval is
//                                          the best which can be
//                                          obtained.
//                                      = 5 the integral over the k th
//                                          cycle is probably divergent
//                                          or slowly convergent. it must
//                                          be noted that divergence can
//                                          occur with any other value of
//                                          iwork(k).
//                    if omega = 0 and integr = 1,
//                    the integral is calculated by means of dqagie,
//                    and ier = iwork(1) (with meaning as described
//                    for iwork(k),k = 1).
//
//         dimensioning parameters
//            limlst - integer
//                     limlst gives an upper bound on the number of
//                     cycles, limlst.ge.3.
//                     if limlst.lt.3, the routine will end with ier = 6.
//
//            lst    - integer
//                     on return, lst indicates the number of cycles
//                     actually needed for the integration.
//                     if omega = 0, then lst is set to 1.
//
//            leniw  - integer
//                     dimensioning parameter for iwork. on entry,
//                     (leniw-limlst)/2 equals the maximum number of
//                     subintervals allowed in the partition of each
//                     cycle, leniw.ge.(limlst+2).
//                     if leniw.lt.(limlst+2), the routine will end with
//                     ier = 6.
//
//            maxp1  - integer
//                     maxp1 gives an upper bound on the number of
//                     chebyshev moments which can be stored, i.e. for
//                     the intervals of lengths abs(b-a)*2**(-l),
//                     l = 0,1, ..., maxp1-2, maxp1.ge.1.
//                     if maxp1.lt.1, the routine will end with ier = 6.
//            lenw   - integer
//                     dimensioning parameter for work
//                     lenw must be at least leniw*2+maxp1*25.
//                     if lenw.lt.(leniw*2+maxp1*25), the routine will
//                     end with ier = 6.
//
//         work arrays
//            iwork  - integer
//                     vector of dimension at least leniw
//                     on return, iwork(k) for k = 1, 2, ..., lst
//                     contain the error flags on the cycles.
//
//            work   - double precision
//                     vector of dimension at least
//                     on return,
//                     work(1), ..., work(lst) contain the integral
//                      approximations over the cycles,
//                     work(limlst+1), ..., work(limlst+lst) contain
//                      the error extimates over the cycles.
//                     further elements of work have no specific
//                     meaning for the user.
//
//***references  (none)
//***routines called  dqawfe,xerror
//***end prologue  dqawf
//
       int last, limit, ll2,l1,l2,l3,l4,l5,l6;
//
//         check validity of limlst, leniw, maxp1 and lenw.
//
//***first executable statement  dqawf
      ier = 6;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      if(limlst < 3 || leniw < (limlst+2) || maxp1 < 1 || lenw <
         (leniw*2+maxp1*25)) goto l10;
//
//         prepare call for dqawfe
//
      limit = (leniw-limlst)/2;
      l1 = limlst;
      l2 = limlst+l1;
      l3 = limit+l2;
      l4 = limit+l3;
      l5 = limit+l4;
      l6 = limit+l5;
      ll2 = limit+l1;
      qawfe!(Real,Func)(f,a,omega,integr,epsabs,limlst,limit,maxp1,result,
        abserr,neval,ier,work,work+l1,iwork,lst,work+l2,
        work+l3,work+l4,work+l5,iwork+l1,iwork+ll2,work+l6);
//
//         call error handler if necessary
//
l10:  if(ier != 0)
        throw new Exception("abnormal return from qawf: "~to!string(ier));
      return;
}


unittest
{
    alias qawf!(float, float delegate(float)) fqawf;
    alias qawf!(double, double delegate(double)) dqawf;
    alias qawf!(double, double function(double)) dfqawf;
    alias qawf!(real, real delegate(real)) rqawf;
}


unittest
{
    double f(double x) { return x > 0.0 ? 1/sqrt(x) : 0.0; }

    enum : double
    {
        a = 0.0,
        omega = PI_2,
        epsabs = 1e-8,
    }
    double result, abserr;
    int neval, ier, lst;
    enum
    {
        integr = 1, // cos(pi x/2)
        limlst = 50,
        leniw = 1050,
        maxp1 = 21,
        lenw = leniw*2 + maxp1*25
    }

    int[leniw] iwork;
    double[lenw] work;

    qawf(&f, a, omega, integr, epsabs, result, abserr, neval, ier, limlst,
        lst, leniw, maxp1, lenw, iwork.ptr, work.ptr);
}


unittest
{
    // This is integral 14 in the QUADPACK book.
    real alpha;
    real f(real x) { return x <= 0.0 ? 0.0 : exp(-(2.0L^^(-alpha))*x)/sqrt(x); }
    enum real omega = 1.0;
    enum integr = 1; // cos(x)
    
    enum : real
    {
        a = 0.0,
        epsabs = 1e-8
    }
    real result, abserr;
    int neval, ier, lst;
    enum
    {
        limlst = 50,
        leniw = 1050,
        maxp1 = 21,
        lenw = leniw*2 + maxp1*25
    }

    int[leniw] iwork;
    real[lenw] work;

    alpha = 0.0;
    qawf(&f, a, omega, integr, epsabs, result, abserr, neval, ier, limlst,
        lst, leniw, maxp1, lenw, iwork.ptr, work.ptr);
    real ans = sqrt(PI) * ((1+(4.0L^^(-alpha)))^^(-0.25L))
        * cos(atan(2.0L^^alpha)/2);
    check (isAccurate(result, abserr, ans, 0.0L, epsabs));
}
