/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qagp;


import std.conv;

import scid.core.testing;
import scid.ports.quadpack.qagpe;

version(unittest)
{
    import std.math;
    import scid.core.testing;
}




///
void qagp(Real, Func)(Func f, Real a, Real b, int npts2, const Real* points,
    Real epsabs, Real epsrel, out Real result, out Real abserr,
    out int neval, out int ier, int leniw, int lenw, out int last,
    int* iwork, Real* work)
{
//***begin prologue  dqagp
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a1
//***keywords  automatic integrator, general-purpose,
//             singularities at user specified points,
//             extrapolation, globally adaptive
//***author  piessens,robert,appl. math. & progr. div - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            definite integral i = integral of f over (a,b),
//            hopefully satisfying following claim for accuracy
//            break points of the integration interval, where local
//            difficulties of the integrand may occur (e.g.
//            singularities, discontinuities), are provided by the user.
//***description
//
//        computation of a definite integral
//        standard fortran subroutine
//        double precision version
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
//            b      - double precision
//                     upper limit of integration
//
//            npts2  - integer
//                     number equal to two more than the number of
//                     user-supplied break points within the integration
//                     range, npts.ge.2.
//                     if npts2.lt.2, the routine will end with ier = 6.
//
//            points - double precision
//                     vector of dimension npts2, the first (npts2-2)
//                     elements of which are the user provided break
//                     points. if these points do not constitute an
//                     ascending sequence there will be an automatic
//                     sorting.
//
//            epsabs - double precision
//                     absolute accuracy requested
//            epsrel - double precision
//                     relative accuracy requested
//                     if  epsabs.le.0
//                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                     the routine will end with ier = 6.
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
//                     ier = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more
//                             subdivisions by increasing the value of
//                             limit (and taking the according dimension
//                             adjustments into account). however, if
//                             this yields no improvement it is advised
//                             to analyze the integrand in order to
//                             determine the integration difficulties. if
//                             the position of a local difficulty can be
//                             determined (i.e. singularity,
//                             discontinuity within the interval), it
//                             should be supplied to the routine as an
//                             element of the vector points. if necessary
//                             an appropriate special-purpose integrator
//                             must be used, which is designed for
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
//                             it is presumed that the requested
//                             tolerance cannot be achieved, and that
//                             the returned result is the best which
//                             can be obtained.
//                         = 5 the integral is probably divergent, or
//                             slowly convergent. it must be noted that
//                             divergence can occur with any other value
//                             of ier.gt.0.
//                         = 6 the input is invalid because
//                             npts2.lt.2 or
//                             break points are specified outside
//                             the integration range or
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
//                             result, abserr, neval, last are set to
//                             zero. exept when leniw or lenw or npts2 is
//                             invalid, iwork(1), iwork(limit+1),
//                             work(limit*2+1) and work(limit*3+1)
//                             are set to zero.
//                             work(1) is set to a and work(limit+1)
//                             to b (where limit = (leniw-npts2)/2).
//
//         dimensioning parameters
//            leniw - integer
//                    dimensioning parameter for iwork
//                    leniw determines limit = (leniw-npts2)/2,
//                    which is the maximum number of subintervals in the
//                    partition of the given integration interval (a,b),
//                    leniw.ge.(3*npts2-2).
//                    if leniw.lt.(3*npts2-2), the routine will end with
//                    ier = 6.
//
//            lenw  - integer
//                    dimensioning parameter for work
//                    lenw must be at least leniw*2-npts2.
//                    if lenw.lt.leniw*2-npts2, the routine will end
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
//                    vector of dimension at least leniw. on return,
//                    the first k elements of which contain
//                    pointers to the error estimates over the
//                    subintervals, such that work(limit*3+iwork(1)),...,
//                    work(limit*3+iwork(k)) form a decreasing
//                    sequence, with k = last if last.le.(limit/2+2), and
//                    k = limit+1-last otherwise
//                    iwork(limit+1), ...,iwork(limit+last) contain the
//                     subdivision levels of the subintervals, i.e.
//                     if (aa,bb) is a subinterval of (p1,p2)
//                     where p1 as well as p2 is a user-provided
//                     break point or integration limit, then (aa,bb) has
//                     level l if abs(bb-aa) = abs(p2-p1)*2**(-l),
//                    iwork(limit*2+1), ..., iwork(limit*2+npts2) have
//                     no significance for the user,
//                    note that limit = (leniw-npts2)/2.
//
//            work  - double precision
//                    vector of dimension at least lenw
//                    on return
//                    work(1), ..., work(last) contain the left
//                     end points of the subintervals in the
//                     partition of (a,b),
//                    work(limit+1), ..., work(limit+last) contain
//                     the right end points,
//                    work(limit*2+1), ..., work(limit*2+last) contain
//                     the integral approximations over the subintervals,
//                    work(limit*3+1), ..., work(limit*3+last)
//                     contain the corresponding error estimates,
//                    work(limit*4+1), ..., work(limit*4+npts2)
//                     contain the integration limits and the
//                     break points sorted in an ascending sequence.
//                    note that limit = (leniw-npts2)/2.
//
//***references  (none)
//***routines called  dqagpe,xerror
//***end prologue  dqagp
//
      int limit,l1,l2,l3,l4;
//
//         check validity of limit and lenw.
//
//***first executable statement  dqagp
      ier = 6;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      if(leniw < (3*npts2-2) || lenw < (leniw*2-npts2) || npts2 < 2)
        goto l10;
//
//         prepare call for dqagpe.
//
      limit = (leniw-npts2)/2;
      l1 = limit;
      l2 = limit+l1;
      l3 = limit+l2;
      l4 = limit+l3;
//
      qagpe!(Real, Func)(f,a,b,npts2,points,epsabs,epsrel,limit,result,abserr,
        neval,ier,work,work+l1,work+l2,work+l3,work+l4,
        iwork,iwork+l1,iwork+l2,last);
//
//         call error handler if necessary.
//
l10:  if(ier != 0)
        throw new Exception("abnormal return from qagp: "~to!string(ier));
      return;
}

unittest
{
    alias qagp!(float, float delegate(float)) fqagp;
    alias qagp!(double, double delegate(double)) dqagp;
    alias qagp!(double, double function(double)) dfqagp;
    alias qagp!(real, real delegate(real)) rqagp;
}

unittest
{

    double f(double x)
    {
        return (x^^3)*log(fabs(((x^^2)-1) * ((x^^2)-2)));
    }

    enum : double
    {
        a = 0.0,
        b = 3.0,
        epsabs = 0.0,
        epsrel = 1e-10,
    }
    double[] points = [1.0, sqrt(2.0)]; // Singularities
    double result, abserr;
    int neval, ier, last;
    enum
    {
        leniw = 1000,
        lenw = leniw*2,
    }
    int[leniw] iwork;
    double[lenw] work;

    qagp(&f, a, b, points.length, points.ptr, epsabs, epsrel,
        result, abserr, neval, ier, leniw, lenw, last, iwork.ptr, work.ptr);

    double ans = 61*log(2.0) + 77*log(7.0)/4 - 27;
    check (isAccurate(result, abserr, ans, epsrel, epsabs));
}
