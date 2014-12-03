// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qawc;


import std.conv;
import scid.ports.quadpack.qawce;

version(unittest)
{
    import std.math;
    import scid.core.testing;
}




///
void qawc(Real, Func)(Func f, Real a, Real b, Real c, Real epsabs, Real epsrel,
    out Real result, out Real abserr, out int neval, out int ier, int limit,
    int lenw, out int last, int* iwork, Real* work)
{
//***begin prologue  dqawc
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a1,j4
//***keywords  automatic integrator, special-purpose,
//             cauchy principal value,
//             clenshaw-curtis, globally adaptive
//***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a
//            cauchy principal value i = integral of f*w over (a,b)
//            (w(x) = 1/((x-c), c.ne.a, c.ne.b), hopefully satisfying
//            following claim for accuracy
//            abs(i-result).le.max(epsabe,epsrel*abs(i)).
//***description
//
//        computation of a cauchy principal value
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
//                     under limit of integration
//
//            b      - double precision
//                     upper limit of integration
//
//            c      - parameter in the weight function, c.ne.a, c.ne.b.
//                     if c = a or c = b, the routine will end with
//                     ier = 6 .
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
//                     estimate or the modulus of the absolute error,
//                     which should equal or exceed abs(i-result)
//
//            neval  - integer
//                     number of integrand evaluations
//
//            ier    - integer
//                     ier = 0 normal and reliable termination of the
//                             routine. it is assumed that the requested
//                             accuracy has been achieved.
//                     ier.gt.0 abnormal termination of the routine
//                             the estimates for integral and error are
//                             less reliable. it is assumed that the
//                             requested accuracy has not been achieved.
//            error messages
//                     ier = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more sub-
//                             divisions by increasing the value of limit
//                             (and taking the according dimension
//                             adjustments into account). however, if
//                             this yields no improvement it is advised
//                             to analyze the integrand in order to
//                             determine the integration difficulties.
//                             if the position of a local difficulty
//                             can be determined (e.g. singularity,
//                             discontinuity within the interval) one
//                             will probably gain from splitting up the
//                             interval at this point and calling
//                             appropriate integrators on the subranges.
//                         = 2 the occurrence of roundoff error is detec-
//                             ted, which prevents the requested
//                             tolerance from being achieved.
//                         = 3 extremely bad integrand behaviour occurs
//                             at some points of the integration
//                             interval.
//                         = 6 the input is invalid, because
//                             c = a or c = b or
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
//                             or limit.lt.1 or lenw.lt.limit*4.
//                             result, abserr, neval, last are set to
//                             zero. exept when lenw or limit is invalid,
//                             iwork(1), work(limit*2+1) and
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
//           lenw   - integer
//                    dimensioning parameter for work
//                    lenw must be at least limit*4.
//                    if lenw.lt.limit*4, the routine will end with
//                    ier = 6.
//
//            last  - integer
//                    on return, last equals the number of subintervals
//                    produced in the subdivision process, which
//                    determines the number of significant elements
//                    actually in the work arrays.
//
//         work arrays
//            iwork - integer
//                    vector of dimension at least limit, the first k
//                    elements of which contain pointers
//                    to the error estimates over the subintervals,
//                    such that work(limit*3+iwork(1)), ... ,
//                    work(limit*3+iwork(k)) form a decreasing
//                    sequence, with k = last if last.le.(limit/2+2),
//                    and k = limit+1-last otherwise
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
//                     contain the error estimates.
//
//***references  (none)
//***routines called  dqawce,xerror
//***end prologue  dqawc
//
      int l1,l2,l3;
//
//         check validity of limit and lenw.
//
//***first executable statement  dqawc
      ier = 6;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      if(limit < 1 || lenw < limit*4) goto l10;
//
//         prepare call for dqawce.
//
      l1 = limit;
      l2 = limit+l1;
      l3 = limit+l2;
      qawce!(Real, Func)(f,a,b,c,epsabs,epsrel,limit,result,abserr,neval,ier,
        work,work+l1,work+l2,work+l3,iwork,last);
//
//         call error handler if necessary.
//
l10:  if(ier != 0)
        throw new Exception("abnormal return from qawc: "~to!string(ier));
      return;
}


unittest
{
    alias qawc!(float, float delegate(float)) fqawc;
    alias qawc!(double, double delegate(double)) dqawc;
    alias qawc!(double, double function(double)) dfqawc;
    alias qawc!(real, real delegate(real)) rqawc;
}


unittest
{
    // Integral 17 in the QUADPACK book.
    enum : double
    {
        a = 0,
        b = 5,
        c = 2,
        epsabs = 0,
        epsrel = 1e-8,
    }
    double result, abserr;
    int neval, ier, last;
    enum limit = 500, lenw = limit*4;
    int[limit] iwork;
    double[lenw] work;

    foreach (alpha; 0 .. 10)
    {
        double f(double x)
        {
            return 2.0^^(-alpha) / ((x-1)^^2 + 4.0^^(-alpha));
        }

        qawc(&f, a, b, c, epsabs, epsrel,
            result, abserr, neval, ier, limit,
            lenw, last, iwork.ptr, work.ptr);

        double exact =
            (
                2.0^^(-alpha) * log(3.0/2)
              - 2.0^^(-alpha-1) * log((16+4.0^^(-alpha))/(1+4.0^^(-alpha)))
              - atan(2.0^^(alpha+2))
              - atan(2.0^^alpha)
            )
          / (1 + 4.0^^(-alpha));

        check(isAccurate(result, abserr, exact, epsrel, epsabs));
    }
}
