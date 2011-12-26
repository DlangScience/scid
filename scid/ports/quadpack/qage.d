// This code has been mechanically translated from the original FORTRAN
// code at http://netlib.org/quadpack.

/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qage;


import std.algorithm: max;
import std.math: fabs;

import scid.ports.quadpack.qk15;
import scid.ports.quadpack.qk21;
import scid.ports.quadpack.qk31;
import scid.ports.quadpack.qk41;
import scid.ports.quadpack.qk51;
import scid.ports.quadpack.qk61;
import scid.ports.quadpack.qpsrt;

import scid.core.fortran;




///
void qage(Real, Func)(Func f, Real a, Real b, Real epsabs, Real epsrel,
    int key, int limit, out Real result, out Real abserr, out int neval,
    out int ier, Real* alist_, Real* blist_, Real* rlist_, Real* elist_,
    int* iord_, out int last)
{
//***begin prologue  dqage
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a1a1
//***keywords  automatic integrator, general-purpose,
//             integrand examinator, globally adaptive,
//             gauss-kronrod
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            definite integral   i = integral of f over (a,b),
//            hopefully satisfying following claim for accuracy
//            abs(i-reslt).le.max(epsabs,epsrel*abs(i)).
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
//            epsabs - double precision
//                     absolute accuracy requested
//            epsrel - double precision
//                     relative accuracy requested
//                     if  epsabs.le.0
//                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                     the routine will end with ier = 6.
//
//            key    - integer
//                     key for choice of local integration rule
//                     a gauss-kronrod pair is used with
//                          7 - 15 points if key.lt.2,
//                         10 - 21 points if key = 2,
//                         15 - 31 points if key = 3,
//                         20 - 41 points if key = 4,
//                         25 - 51 points if key = 5,
//                         30 - 61 points if key.gt.5.
//
//            limit  - integer
//                     gives an upperbound on the number of subintervals
//                     in the partition of (a,b), limit.ge.1.
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
//            error messages
//                     ier = 1 maximum number of subdivisions allowed
//                             has been achieved. one can allow more
//                             subdivisions by increasing the value
//                             of limit.
//                             however, if this yields no improvement it
//                             is rather advised to analyze the integrand
//                             in order to determine the integration
//                             difficulties. if the position of a local
//                             difficulty can be determined(e.g.
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
//                         = 3 extremely bad integrand behaviour occurs
//                             at some points of the integration
//                             interval.
//                         = 6 the input is invalid, because
//                             (epsabs.le.0 and
//                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28),
//                             result, abserr, neval, last, rlist(1) ,
//                             elist(1) and iord(1) are set to zero.
//                             alist(1) and blist(1) are set to a and b
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
//                       last  elements of which are the
//                      integral approximations on the subintervals
//
//            elist   - double precision
//                      vector of dimension at least limit, the first
//                       last  elements of which are the moduli of the
//                      absolute error estimates on the subintervals
//
//            iord    - integer
//                      vector of dimension at least limit, the first k
//                      elements of which are pointers to the
//                      error estimates over the subintervals,
//                      such that elist(iord(1)), ...,
//                      elist(iord(k)) form a decreasing sequence,
//                      with k = last if last.le.(limit/2+2), and
//                      k = limit+1-last otherwise
//
//            last    - integer
//                      number of subintervals actually produced in the
//                      subdivision process
//
//***references  (none)
//***routines called  d1mach,dqk15,dqk21,dqk31,
//                    dqk41,dqk51,dqk61,dqpsrt
//***end prologue  dqage
//
      Real area,area1,area12,area2,a1,a2,
        b1,b2,defabs,defab1,defab2,epmach,
        errbnd,errmax,error1,error2,erro12,errsum,
        resabs,uflow;
      int iroff1=1,iroff2=1,k=1,keyf=1,maxerr=1, nrmax;
//
      auto alist = dimension(alist_, limit);
      auto blist = dimension(blist_, limit);
      auto elist = dimension(elist_, limit);
      auto iord = dimension(iord_, limit);
      auto rlist = dimension(rlist_, limit);
//
//            list of major variables
//            -----------------------
//
//           alist     - list of left end points of all subintervals
//                       considered up to now
//           blist     - list of right end points of all subintervals
//                       considered up to now
//           rlist(i)  - approximation to the integral over
//                      (alist(i),blist(i))
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
//           machine dependent constants
//           ---------------------------
//
//           epmach  is the largest relative spacing.
//           uflow  is the smallest positive magnitude.
//
//***first executable statement  dqage
      epmach = Real.epsilon;
      uflow = Real.epsilon;
//
//           test on validity of parameters
//           ------------------------------
//
      ier = 0;
      neval = 0;
      last = 0;
      result = 0.0;
      abserr = 0.0;
      alist[1] = a;
      blist[1] = b;
      rlist[1] = 0.0;
      elist[1] = 0.0;
      iord[1] = 0;
      if(epsabs <= 0.0 && 
        epsrel < max(0.5e2*epmach,0.5e-28)) ier = 6;
      if(ier == 6) goto l999;
//
//           first approximation to the integral
//           -----------------------------------
//
      keyf = key;
      if(key <= 0) keyf = 1;
      if(key >= 7) keyf = 6;
      neval = 0;
      typeof(&qk15!(Real, Func)) qkXX;
      if(keyf == 1)      qkXX = &qk15!(Real, Func);
      else if(keyf == 2) qkXX = &qk21!(Real, Func);
      else if(keyf == 3) qkXX = &qk31!(Real, Func);
      else if(keyf == 4) qkXX = &qk41!(Real, Func);
      else if(keyf == 5) qkXX = &qk51!(Real, Func);
      else if(keyf == 6) qkXX = &qk61!(Real, Func);
      qkXX(f,a,b,result,abserr,defabs,resabs);
      last = 1;
      rlist[1] = result;
      elist[1] = abserr;
      iord[1] = 1;
//
//           test on accuracy.
//
      errbnd = max(epsabs,epsrel*fabs(result));
      if(abserr <= 0.5e2*epmach*defabs && abserr > errbnd) ier = 2;
      if(limit == 1) ier = 1;
      if(ier != 0 || (abserr <= errbnd && abserr != resabs)
         || abserr == 0.0) goto l60;
//
//           initialization
//           --------------
//
//
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
      for (last=2; last<=limit; last++) { //do 30 last = 2,limit
//
//           bisect the subinterval with the largest error estimate.
//
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr]+blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        qkXX(f,a1,b1,area1,error1,resabs,defab1);
        qkXX(f,a2,b2,area2,error2,resabs,defab2);
//
//           improve previous approximations to integral
//           and error and test for accuracy.
//
        neval = neval+1;
        area12 = area1+area2;
        erro12 = error1+error2;
        errsum = errsum+erro12-errmax;
        area = area+area12-rlist[maxerr];
        if(defab1 == error1 || defab2 == error2) goto l5;
        if(fabs(rlist[maxerr]-area12) <= 0.1e-4*fabs(area12)
         && erro12 >= 0.99*errmax) iroff1 = iroff1+1;
        if(last > 10 && erro12 > errmax) iroff2 = iroff2+1;
  l5:   rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel*fabs(area));
        if(errsum <= errbnd) goto l8;
//
//           test for roundoff error and eventually set error flag.
//
        if(iroff1 >= 6 || iroff2 >= 20) ier = 2;
//
//           set error flag in the case that the number of subintervals
//           equals limit.
//
        if(last == limit) ier = 1;
//
//           set error flag in the case of bad integrand behaviour
//           at a point of the integration range.
//
        if(max(fabs(a1),fabs(b2)) <= (0.1e1+0.1e3*
          epmach)*(fabs(a2)+0.1e4*uflow)) ier = 3;
//
//           append the newly-created intervals to the list.
//
  l8:   if(error2 > error1) goto l10;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto l20;
 l10:   alist[maxerr] = a2;
        alist[last] = a1;
        blist[last] = b1;
        rlist[maxerr] = area2;
        rlist[last] = area1;
        elist[maxerr] = error2;
        elist[last] = error1;
//
//           call subroutine dqpsrt to maintain the descending ordering
//           in the list of error estimates and select the subinterval
//           with the largest error estimate (to be bisected next).
//
 l20:   qpsrt!Real(limit,last,maxerr,errmax,elist.ptr,iord.ptr,nrmax);
// ***jump out of do-loop
        if(ier != 0 || errsum <= errbnd) goto l40;
 l30: ;}
//
//           compute final result.
//           ---------------------
//
 l40: result = 0.0;
      for (k=1; k<=last; k++) { //do 50 k=1,last
        result = result+rlist[k];
 l50: ;}
      abserr = errsum;
 l60: if(keyf != 1) neval = (10*keyf+1)*(2*neval+1);
      else neval = 30*neval+15;
l999: return;
}

unittest
{
    alias qage!(float, float delegate(float)) fqage;
    alias qage!(double, double delegate(double)) dqage;
    alias qage!(double, double function(double)) dfqage;
    alias qage!(real, real delegate(real)) rqage;
}
