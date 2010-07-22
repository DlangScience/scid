/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qag;


import std.conv;
import scid.ports.quadpack.qage;

version(unittest)
{
    import std.math;
    import scid.core.testing;
}




///
void qag(Real, Func)(Func f, Real a, Real b, Real epsabs, Real epsrel, int key,
    out Real result, out Real abserr, out int neval, out int ier, 
    int limit, int lenw, out int last, int* iwork, Real* work)
{
//***begin prologue  dqag
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a1a1
//***keywords  automatic integrator, general-purpose,
//             integrand examinator, globally adaptive,
//             gauss-kronrod
//***author  piessens,robert,appl. math. & progr. div - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            definite integral i = integral of f over (a,b),
//            hopefully satisfying following claim for accuracy
//            abs(i-result)le.max(epsabs,epsrel*abs(i)).
//***description
//
//        computation of a definite integral
//        standard fortran subroutine
//        double precision version
//
//            f      - double precision
//                     function subprogam defining the integrand
//                     function f(x). the actual name for f needs to be
//                     declared e x t e r n a l in the driver program.
//
//            a      - double precision
//                     lower limit of integration
//
//            b      - double precision
//                     upper limit of integration
//
//            epsabs - double precision
//                     absolute accoracy requested
//            epsrel - double precision
//                     relative accuracy requested
//                     if  epsabs.le.0
//                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                     the routine will end with ier = 6.
//
//            key    - integer
//                     key for choice of local integration rule
//                     a gauss-kronrod pair is used with
//                       7 - 15 points if key.lt.2,
//                      10 - 21 points if key = 2,
//                      15 - 31 points if key = 3,
//                      20 - 41 points if key = 4,
//                      25 - 51 points if key = 5,
//                      30 - 61 points if key.gt.5.
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
//                     ier.gt.0 abnormal termination of the routine
//                             the estimates for result and error are
//                             less reliable. it is assumed that the
//                             requested accuracy has not been achieved.
//                      error messages
//                     ier = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more
//                             subdivisions by increasing the value of
//                             limit (and taking the according dimension
//                             adjustments into account). however, if
//                             this yield no improvement it is advised
//                             to analyze the integrand in order to
//                             determine the integration difficulaties.
//                             if the position of a local difficulty can
//                             be determined (i.e.singularity,
//                             discontinuity within the interval) one
//                             will probably gain from splitting up the
//                             interval at this point and calling the
//                             integrator on the subranges. if possible,
//                             an appropriate special-purpose integrator
//                             should be used which is designed for
//                             handling the type of difficulty involved.
//                         = 2 the occurrence of roundoff error is
//                             detected, which prevents the requested
//                             tolerance from being achieved.
//                         = 3 extremely bad integrand behaviour occurs
//                             at some points of the integration
//                             interval.
//                         = 6 the input is invalid, because
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
//                             or limit.lt.1 or lenw.lt.limit*4.
//                             result, abserr, neval, last are set
//                             to zero.
//                             except when lenw is invalid, iwork(1),
//                             work(limit*2+1) and work(limit*3+1) are
//                             set to zero, work(1) is set to a and
//                             work(limit+1) to b.
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
//                    if lenw.lt.limit*4, the routine will end with
//                    ier = 6.
//
//            last  - integer
//                    on return, last equals the number of subintervals
//                    produced in the subdiviosion process, which
//                    determines the number of significant elements
//                    actually in the work arrays.
//
//         work arrays
//            iwork - integer
//                    vector of dimension at least limit, the first k
//                    elements of which contain pointers to the error
//                    estimates over the subintervals, such that
//                    work(limit*3+iwork(1)),... , work(limit*3+iwork(k))
//                    form a decreasing sequence with k = last if
//                    last.le.(limit/2+2), and k = limit+1-last otherwise
//
//            work  - double precision
//                    vector of dimension at least lenw
//                    on return
//                    work(1), ..., work(last) contain the left end
//                    points of the subintervals in the partition of
//                     (a,b),
//                    work(limit+1), ..., work(limit+last) contain the
//                     right end points,
//                    work(limit*2+1), ..., work(limit*2+last) contain
//                     the integral approximations over the subintervals,
//                    work(limit*3+1), ..., work(limit*3+last) contain
//                     the error estimates.
//
//***references  (none)
//***routines called  dqage,xerror
//***end prologue  dqag
      int l1,l2,l3;
//
//         check validity of lenw.
//
//***first executable statement  dqag
      ier = 6;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      if(limit < 1 || lenw < limit*4) goto l10;
//
//         prepare call for dqage.
//
      l1 = limit;
      l2 = limit+l1;
      l3 = limit+l2;
//
      qage!(Real, Func)(f,a,b,epsabs,epsrel,key,limit,result,abserr,neval,
        ier,work,work+l1,work+l2,work+l3,iwork,last);
//
//         call error handler if necessary.
//
l10:  
      if(ier != 0)
        throw new Exception("abnormal return from qag: "~to!string(ier));
      return;
}

unittest
{
    alias qag!(float, float delegate(float)) fqag;
    alias qag!(double, double delegate(double)) dqag;
    alias qag!(double, double function(double)) dfqag;
    alias qag!(real, real delegate(real)) rqag;
}

unittest
{
    double f(double x) { return cos(100*sin(x)); }

    enum : double
    {
        a = 0.0,
        b = PI,
        epsabs = 0.0
    }
    double epsrel=1e-10, result, abserr;
    int key, neval, ier, last;
    enum
    {
        limit = 500,
        lenw = 4*limit
    }

    int[limit] iwork;
    double[lenw] work;

    // The answer is pi*J_0(100)
    enum double ans = 0.062787400491492695655;

    // Check that the routine works with all keys.
    for (key=1; key<=6; key++)
    {
        qag(&f, a, b, epsabs, epsrel, key, result, abserr, neval, ier,
            limit, lenw, last, iwork.ptr, work.ptr);
        check (isAccurate(result, abserr, ans, epsrel, epsabs));
    }
}

unittest
{
    // Integral 1 in the QUADPACK book.
    double alpha;
    double f(double x) { return (x^^alpha)*log(1/x); }

    double a=0.0, b=1.0, epsabs=0.0, epsrel=1e-8;
    double result, abserr;
    int neval, ier, last;
    enum { limit=500, lenw=4*limit }
    int[limit] iwork;
    double[lenw] work;

    for (alpha=0.0; alpha<=2.6+10*double.epsilon;
        alpha+=0.2)
    {
        double ans = (alpha+1)^^(-2.0);
        foreach (key; [1,3,6])
        {
            qag(&f, a, b, epsabs, epsrel, key, result, abserr, neval, ier,
                limit, lenw, last, iwork.ptr, work.ptr);

            check (isAccurate(result, abserr, ans, epsrel, epsabs));
        }
    }
}

unittest
{
    // Integral 15 in the QUADPACK book.
    double alpha;
    double f(double x) { return (x^^2.0) * exp(-(2.0^^(-alpha))*x); }
    
    double a = 0.0, epsabs = 0.0, epsrel = 1e-8, result, abserr;
    int key = 6, neval, ier, last;
    enum { limit = 500, lenw = 4*limit }

    int[limit] iwork;
    double[lenw] work;

    for (alpha=0.0; alpha<5.001; alpha+=1.0)
    {
        double b = 40*(2.0^^alpha);
        qag (&f, a, b, epsabs, epsrel, key, result, abserr, neval, ier,
            limit, lenw, last, iwork.ptr, work.ptr);
        
        double eps = 1e-15;
        double ans = (1+eps)*(2.0^^(3*alpha+1.0));
        check (isAccurate(result, abserr, ans, epsrel, epsabs));
    }
}
