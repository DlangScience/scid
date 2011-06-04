/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qawce;


import std.algorithm: min, max;
import std.math: fabs;

import scid.common.fortran;
import scid.ports.quadpack.qc25c;
import scid.ports.quadpack.qpsrt;




///
void qawce(Real, Func)(Func f, Real a, Real b, Real c, Real epsabs,
    Real epsrel, int limit, out Real result,out Real abserr, out int neval,
    out int ier, Real* alist_, Real* blist_, Real* rlist_, Real* elist_,
    int* iord_, out int last)
{
//***begin prologue  dqawce
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a1,j4
//***keywords  automatic integrator, special-purpose,
//             cauchy principal value, clenshaw-curtis method
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***  purpose  the routine calculates an approximation result to a
//              cauchy principal value i = integral of f*w over (a,b)
//              (w(x) = 1/(x-c), (c.ne.a, c.ne.b), hopefully satisfying
//              following claim for accuracy
//              abs(i-result).le.max(epsabs,epsrel*abs(i))
//***description
//
//        computation of a cauchy principal value
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
//            c      - double precision
//                     parameter in the weight function, c.ne.a, c.ne.b
//                     if c = a or c = b, the routine will end with
//                     ier = 6.
//
//            epsabs - double precision
//                     absolute accuracy requested
//            epsrel - double precision
//                     relative accuracy requested
//                     if  epsabs.le.0
//                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                     the routine will end with ier = 6.
//
//            limit  - integer
//                     gives an upper bound on the number of subintervals
//                     in the partition of (a,b), limit.ge.1
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
//                             the estimates for integral and error are
//                             less reliable. it is assumed that the
//                             requested accuracy has not been achieved.
//            error messages
//                     ier = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more sub-
//                             divisions by increasing the value of
//                             limit. however, if this yields no
//                             improvement it is advised to analyze the
//                             the integrand, in order to determine the
//                             the integration difficulties. if the
//                             position of a local difficulty can be
//                             determined (e.g. singularity,
//                             discontinuity within the interval) one
//                             will probably gain from splitting up the
//                             interval at this point and calling
//                             appropriate integrators on the subranges.
//                         = 2 the occurrence of roundoff error is detec-
//                             ted, which prevents the requested
//                             tolerance from being achieved.
//                         = 3 extremely bad integrand behaviour
//                             occurs at some interior points of
//                             the integration interval.
//                         = 6 the input is invalid, because
//                             c = a or c = b or
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
//                             or limit.lt.1.
//                             result, abserr, neval, rlist(1), elist(1),
//                             iord(1) and last are set to zero. alist(1)
//                             and blist(1) are set to a and b
//                             respectively.
//
//            alist   - double precision
//                      vector of dimension at least limit, the first
//                       last  elements of which are the left
//                      end points of the subintervals in the partition
//                      of the given integration range (a,b)
//
//            blist   - double precision
//                      vector of dimension at least limit, the first
//                       last  elements of which are the right
//                      end points of the subintervals in the partition
//                      of the given integration range (a,b)
//
//            rlist   - double precision
//                      vector of dimension at least limit, the first
//                       last  elements of which are the integral
//                      approximations on the subintervals
//
//            elist   - double precision
//                      vector of dimension limit, the first  last
//                      elements of which are the moduli of the absolute
//                      error estimates on the subintervals
//
//            iord    - integer
//                      vector of dimension at least limit, the first k
//                      elements of which are pointers to the error
//                      estimates over the subintervals, so that
//                      elist(iord(1)), ..., elist(iord(k)) with k = last
//                      if last.le.(limit/2+2), and k = limit+1-last
//                      otherwise, form a decreasing sequence
//
//            last    - integer
//                      number of subintervals actually produced in
//                      the subdivision process
//
//***references  (none)
//***routines called  d1mach,dqc25c,dqpsrt
//***end prologue  dqawce
//
      Real aa,area,area1,area12,area2,a1,a2,
        bb,b1,b2,epmach,
        errbnd,errmax,error1,erro12,error2,errsum,uflow;
      int iroff1,iroff2,k,krule,maxerr,nev,
        nrmax;
//
      auto alist = dimension(alist_, limit);
      auto blist = dimension(blist_, limit);
      auto rlist = dimension(rlist_, limit);
      auto elist = dimension(elist_, limit);
      auto iord  = dimension(iord_, limit);
//
//            list of major variables
//            -----------------------
//
//           alist     - list of left end points of all subintervals
//                       considered up to now
//           blist     - list of right end points of all subintervals
//                       considered up to now
//           rlist(i)  - approximation to the integral over
//                       (alist(i),blist(i))
//           elist(i)  - error estimate applying to rlist(i)
//           maxerr    - pointer to the interval with largest
//                       error estimate
//           errmax    - elist(maxerr)
//           area      - sum of the integrals over the subintervals
//           errsum    - sum of the errors over the subintervals
//           errbnd    - requested accuracy max(epsabs,epsrel*
//                       abs(result))
//           *****1    - variable for the left subinterval
//           *****2    - variable for the right subinterval
//           last      - index for subdivision
//
//
//            machine dependent constants
//            ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//
//***first executable statement  dqawce
      epmach = Real.epsilon;
      uflow = Real.min_normal;
//
//
//           test on validity of parameters
//           ------------------------------
//
      ier = 6;
      neval = 0;
      last = 0;
      alist[1] = a;
      blist[1] = b;
      rlist[1] = 0.0;
      elist[1] = 0.0;
      iord[1] = 0;
      result = 0.0;
      abserr = 0.0;
      if(c == a || c == b || (epsabs <= 0.0 &&
        epsrel < max(0.5e2*epmach,0.5e-28))) goto l999;
//
//           first approximation to the integral
//           -----------------------------------
//
      aa=a;
      bb=b;
      if (a <= b) goto l10;
      aa=b;
      bb=a;
 l10: ier=0;
      krule = 1;
      qc25c!(Real,Func)(f,aa,bb,c,result,abserr,krule,neval);
      last = 1;
      rlist[1] = result;
      elist[1] = abserr;
      iord[1] = 1;
      alist[1] = a;
      blist[1] = b;
//
//           test on accuracy
//
      errbnd = max(epsabs,epsrel*fabs(result));
      if(limit == 1) ier = 1;
      if(abserr < min(0.1e-1*fabs(result),errbnd)
         || ier == 1) goto l70;
//
//           initialization
//           --------------
//
      alist[1] = aa;
      blist[1] = bb;
      rlist[1] = result;
      errmax = abserr;
      maxerr = 1;
      area = result;
      errsum = abserr;
      nrmax = 1;
      iroff1 = 0;
      iroff2 = 0;
//
//           main do-loop
//           ------------
//
      for (last=2; last<=limit; last++) { //do 40 last = 2,limit
//
//           bisect the subinterval with nrmax-th largest
//           error estimate.
//
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr]+blist[maxerr]);
        b2 = blist[maxerr];
        if(c <= b1 && c > a1) b1 = 0.5*(c+b2);
        if(c > b1 && c < b2) b1 = 0.5*(a1+c);
        a2 = b1;
        krule = 2;
        qc25c!(Real, Func)(f,a1,b1,c,area1,error1,krule,nev);
        neval = neval+nev;
        qc25c!(Real, Func)(f,a2,b2,c,area2,error2,krule,nev);
        neval = neval+nev;
//
//           improve previous approximations to integral
//           and error and test for accuracy.
//
        area12 = area1+area2;
        erro12 = error1+error2;
        errsum = errsum+erro12-errmax;
        area = area+area12-rlist[maxerr];
        if(fabs(rlist[maxerr]-area12) < 0.1e-4*fabs(area12)
           && erro12 >= 0.99*errmax && krule == 0)
          iroff1 = iroff1+1;
        if(last > 10 && erro12 > errmax && krule == 0)
          iroff2 = iroff2+1;
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel*fabs(area));
        if(errsum <= errbnd) goto l15;
//
//           test for roundoff error and eventually set error flag.
//
        if(iroff1 >= 6 && iroff2 > 20) ier = 2;
//
//           set error flag in the case that number of interval
//           bisections exceeds limit.
//
        if(last == limit) ier = 1;
//
//           set error flag in the case of bad integrand behaviour
//           at a point of the integration range.
//
        if(max(fabs(a1),fabs(b2)) <= (0.1e1+0.1e3*epmach)
          *(fabs(a2)+0.1e4*uflow)) ier = 3;
//
//           append the newly-created intervals to the list.
//
 l15:   if(error2 > error1) goto l20;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto l30;
 l20:   alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;
//
//           call subroutine dqpsrt to maintain the descending ordering
//           in the list of error estimates and select the subinterval
//           with nrmax-th largest error estimate (to be bisected next).
//
 l30:    qpsrt!Real(limit,last,maxerr,errmax,elist.ptr,iord.ptr,nrmax);
// ***jump out of do-loop
        if(ier != 0 || errsum <= errbnd) goto l50;
 l40: ;}
//
//           compute final result.
//           ---------------------
//
 l50: result = 0.0;
      for (k=1; k<=last; k++) { //do 60 k=1,last
        result = result+rlist[k];
 l60: ;}
      abserr = errsum;
 l70: if (aa == b) result=-result;
l999: return;
}


unittest
{
    alias qawce!(float, float delegate(float)) fqawce;
    alias qawce!(double, double delegate(double)) dqawce;
    alias qawce!(double, double function(double)) dfqawce;
    alias qawce!(real, real delegate(real)) rqawce;
}
