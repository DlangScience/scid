// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qagi;


import std.conv;
import scid.ports.quadpack.qagie;

version(unittest)
{
    import std.math;
    import scid.core.testing;
}




///
void qagi(Real, Func)(Func f, Real bound, int inf, Real epsabs,
    Real epsrel, out Real result, out Real abserr, out int neval,
    out int ier, int limit, int lenw, out int last,
    int* iwork, Real* work)
{
//***begin prologue  dqagi
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a3a1,h2a4a1
//***keywords  automatic integrator, infinite intervals,
//             general-purpose, transformation, extrapolation,
//             globally adaptive
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. -k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            integral   i = integral of f over (bound,+infinity)
//            or i = integral of f over (-infinity,bound)
//            or i = integral of f over (-infinity,+infinity)
//            hopefully satisfying following claim for accuracy
//            abs(i-result).le.max(epsabs,epsrel*abs(i)).
//***description
//
//        integration over infinite intervals
//        standard fortran subroutine
//
//        parameters
//         on entry
//            f      - double precision
//                     function subprogram defining the integrand
//                     function f(x). the actual name for f needs to be
//                     declared e x t e r n a l in the driver program.
//
//            bound  - double precision
//                     finite bound of integration range
//                     (has no meaning if interval is doubly-infinite)
//
//            inf    - integer
//                     indicating the kind of integration range involved
//                     inf = 1 corresponds to  (bound,+infinity),
//                     inf = -1            to  (-infinity,bound),
//                     inf = 2             to (-infinity,+infinity).
//
//            epsabs - double precision
//                     absolute accuracy requested
//            epsrel - double precision
//                     relative accuracy requested
//                     if  epsabs.le.0
//                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                     the routine will end with ier = 6.
//
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
//                   - ier.gt.0 abnormal termination of the routine. the
//                             estimates for result and error are less
//                             reliable. it is assumed that the requested
//                             accuracy has not been achieved.
//            error messages
//                     ier = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more
//                             subdivisions by increasing the value of
//                             limit (and taking the according dimension
//                             adjustments into account). however, if
//                             this yields no improvement it is advised
//                             to analyze the integrand in order to
//                             determine the integration difficulties. if
//                             the position of a local difficulty can be
//                             determined (e.g. singularity,
//                             discontinuity within the interval) one
//                             will probably gain from splitting up the
//                             interval at this point and calling the
//                             integrator on the subranges. if possible,
//                             an appropriate special-purpose integrator
//                             should be used, which is designed for
//                             handling the type of difficulty involved.
//                         = 2 the occurrence of roundoff error is
//                             detected, which prevents the requested
//                             tolerance from being achieved.
//                             the error may be under-estimated.
//                         = 3 extremely bad integrand behaviour occurs
//                             at some points of the integration
//                             interval.
//                         = 4 the algorithm does not converge.
//                             roundoff error is detected in the
//                             extrapolation table.
//                             it is assumed that the requested tolerance
//                             cannot be achieved, and that the returned
//                             result is the best which can be obtained.
//                         = 5 the integral is probably divergent, or
//                             slowly convergent. it must be noted that
//                             divergence can occur with any other value
//                             of ier.
//                         = 6 the input is invalid, because
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
//                              or limit.lt.1 or leniw.lt.limit*4.
//                             result, abserr, neval, last are set to
//                             zero. exept when limit or leniw is
//                             invalid, iwork(1), work(limit*2+1) and
//                             work(limit*3+1) are set to zero, work(1)
//                             is set to a and work(limit+1) to b.
//
//         dimensioning parameters
//            limit - integer
//                    dimensioning parameter for iwork
//                    limit determines the maximum number of subintervals
//                    in the partition of the given integration interval
//                    (a,b), limit.ge.1.
//                    if limit.lt.1, the routine will end with ier = 6.
//
//            lenw  - integer
//                    dimensioning parameter for work
//                    lenw must be at least limit*4.
//                    if lenw.lt.limit*4, the routine will end
//                    with ier = 6.
//
//            last  - integer
//                    on return, last equals the number of subintervals
//                    produced in the subdivision process, which
//                    determines the number of significant elements
//                    actually in the work arrays.
//
//         work arrays
//            iwork - integer
//                    vector of dimension at least limit, the first
//                    k elements of which contain pointers
//                    to the error estimates over the subintervals,
//                    such that work(limit*3+iwork(1)),... ,
//                    work(limit*3+iwork(k)) form a decreasing
//                    sequence, with k = last if last.le.(limit/2+2), and
//                    k = limit+1-last otherwise
//
//            work  - double precision
//                    vector of dimension at least lenw
//                    on return
//                    work(1), ..., work(last) contain the left
//                     end points of the subintervals in the
//                     partition of (a,b),
//                    work(limit+1), ..., work(limit+last) contain
//                     the right end points,
//                    work(limit*2+1), ...,work(limit*2+last) contain the
//                     integral approximations over the subintervals,
//                    work(limit*3+1), ..., work(limit*3)
//                     contain the error estimates.
//***references  (none)
//***routines called  dqagie,xerror
//***end prologue  dqagi
//
      int lvl=1,l1=1,l2=1,l3=1;
//
//         check validity of limit and lenw.
//
//***first executable statement  dqagi
      ier = 6;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      if(limit < 1 || lenw < limit*4) goto l10;
//
//         prepare call for dqagie.
//
      l1 = limit;
      l2 = limit+l1;
      l3 = limit+l2;
//
      qagie!(Real,Func)(f,bound,inf,epsabs,epsrel,limit,result,abserr,
       neval,ier,work,work+l1,work+l2,work+l3,iwork,last);
//
//         call error handler if necessary.
//
       lvl = 0;
l10:  if(ier == 6) lvl = 1;
      if(ier != 0)
        throw new Exception("abnormal return from  qagi: "~to!string(ier));
      return;
}


unittest
{
    alias qagi!(float, float delegate(float)) fqagi;
    alias qagi!(double, double delegate(double)) dqagi;
    alias qagi!(double, double function(double)) dfqagi;
    alias qagi!(real, real delegate(real)) rqagi;
}

unittest
{
    real f(real x) { return exp(-x*x); }

    enum : real
    {
        bound = 0.0,
        epsabs = 0.0,
        epsrel = 1e-10
    }
    real result, abserr;
    int neval, ier;
    enum
    {
        limit = 500,
        lenw = 4*limit
    }
    int last;
    
    int[limit] iwork;
    real[lenw] work;

    enum real sqrtPi = sqrt(PI);

    int inf = 1; // (bound,+inf)
    qagi(&f, bound, inf, epsabs, epsrel, result, abserr, neval, ier,
        limit, lenw, last, iwork.ptr, work.ptr);
    assert (isAccurate(result, abserr, sqrtPi/2, epsrel, epsabs));

    inf = -1; // (-inf,bound)
    qagi(&f, bound, inf, epsabs, epsrel, result, abserr, neval, ier,
        limit, lenw, last, iwork.ptr, work.ptr);
    assert (isAccurate(result, abserr, sqrtPi/2, epsrel, epsabs));
    
    inf = 2; // (-inf,+inf)
    qagi(&f, bound, inf, epsabs, epsrel, result, abserr, neval, ier,
        limit, lenw, last, iwork.ptr, work.ptr);

    assert (isAccurate(result, abserr, sqrtPi, epsrel, epsabs));
}

unittest
{
    // Integral 15 in the QUADPACK book.
    double alpha;
    double f(double x) { return (x^^2.0) * exp(-(2.0^^(-alpha))*x); }
    
    double bound = 0.0, epsabs = 0.0, epsrel = 1e-8, result, abserr;
    int inf = 1, neval, ier, last;
    enum { limit = 500, lenw = 4*limit }

    int[limit] iwork;
    double[lenw] work;

    for (alpha=0.0; alpha<5.001; alpha+=1.0)
    {
        qagi (&f, bound, inf, epsabs, epsrel, result, abserr, neval, ier,
            limit, lenw, last, iwork.ptr, work.ptr);

        double ans = 2.0^^(3*alpha+1.0);
        assert (isAccurate(result, abserr, ans, epsrel, epsabs));
    }
}


unittest
{
    // Integral 16 in the QUADPACK book.
    double alpha;
    double f(double x) { return (x^^(alpha-1))/((1+10*x)^^2); }
    
    double bound = 0.0, epsabs = 0.0, epsrel = 1e-8, result, abserr;
    int inf = 1, neval, ier, last;
    enum { limit = 500, lenw = 4*limit }

    int[limit] iwork;
    double[lenw] work;

    for (int i=-9; i<=9; i++)
    {
        alpha = 1.0 + i*0.1;
        qagi (&f, bound, inf, epsabs, epsrel, result, abserr, neval, ier,
            limit, lenw, last, iwork.ptr, work.ptr);

        double ans;
        if (alpha == 1.0)  ans = 0.1;
        else ans = (10.0^^(-alpha))*(1-alpha)*PI/sin(PI*alpha);

        assert (isAccurate(result, abserr, ans, epsrel, epsabs));
    }
}
