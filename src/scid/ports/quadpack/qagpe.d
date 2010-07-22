/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qagpe;


import std.algorithm: min, max;
import std.math: fabs;

import scid.ports.quadpack.qelg;
import scid.ports.quadpack.qk21;
import scid.ports.quadpack.qpsrt;

import scid.core.fortran;




///
void qagpe(Real, Func) (Func f, Real a, Real b, int npts2, const Real* points_,
    Real epsabs, Real epsrel, int limit, out Real result, out Real abserr,
    out int neval, out int ier, Real* alist_, Real* blist_, Real* rlist_,
    Real* elist_, Real* pts_, int* iord_, int* level_, int* ndin_,
    out int last)
{
//***begin prologue  dqagpe
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a1
//***keywords  automatic integrator, general-purpose,
//             singularities at user specified points,
//             extrapolation, globally adaptive.
//***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a given
//            definite integral i = integral of f over (a,b), hopefully
//            satisfying following claim for accuracy abs(i-result).le.
//            max(epsabs,epsrel*abs(i)). break points of the integration
//            interval, where local difficulties of the integrand may
//            occur(e.g. singularities,discontinuities),provided by user.
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
//                     range, npts2.ge.2.
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
//            limit  - integer
//                     gives an upper bound on the number of subintervals
//                     in the partition of (a,b), limit.ge.npts2
//                     if limit.lt.npts2, the routine will end with
//                     ier = 6.
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
//                             extrapolation table. it is presumed that
//                             the requested tolerance cannot be
//                             achieved, and that the returned result is
//                             the best which can be obtained.
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
//                             or limit.lt.npts2.
//                             result, abserr, neval, last, rlist(1),
//                             and elist(1) are set to zero. alist(1) and
//                             blist(1) are set to a and b respectively.
//
//            alist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the left end points
//                     of the subintervals in the partition of the given
//                     integration range (a,b)
//
//            blist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the right end points
//                     of the subintervals in the partition of the given
//                     integration range (a,b)
//
//            rlist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the integral
//                     approximations on the subintervals
//
//            elist  - double precision
//                     vector of dimension at least limit, the first
//                      last  elements of which are the moduli of the
//                     absolute error estimates on the subintervals
//
//            pts    - double precision
//                     vector of dimension at least npts2, containing the
//                     integration limits and the break points of the
//                     interval in ascending sequence.
//
//            level  - integer
//                     vector of dimension at least limit, containing the
//                     subdivision levels of the subinterval, i.e. if
//                     (aa,bb) is a subinterval of (p1,p2) where p1 as
//                     well as p2 is a user-provided break point or
//                     integration limit, then (aa,bb) has level l if
//                     abs(bb-aa) = abs(p2-p1)*2**(-l).
//
//            ndin   - integer
//                     vector of dimension at least npts2, after first
//                     integration over the intervals (pts(i)),pts(i+1),
//                     i = 0,1, ..., npts2-2, the error estimates over
//                     some of the intervals may have been increased
//                     artificially, in order to put their subdivision
//                     forward. if this happens for the subinterval
//                     numbered k, ndin(k) is put to 1, otherwise
//                     ndin(k) = 0.
//
//            iord   - integer
//                     vector of dimension at least limit, the first k
//                     elements of which are pointers to the
//                     error estimates over the subintervals,
//                     such that elist(iord(1)), ..., elist(iord(k))
//                     form a decreasing sequence, with k = last
//                     if last.le.(limit/2+2), and k = limit+1-last
//                     otherwise
//
//            last   - integer
//                     number of subintervals actually produced in the
//                     subdivisions process
//
//***references  (none)
//***routines called  d1mach,dqelg,dqk21,dqpsrt
//***end prologue  dqagpe
      Real abseps,area,area1,area12,area2,a1,
        a2,b1,b2,correc,defabs,defab1,defab2,
        dres,epmach,erlarg,erlast,errbnd,
        errmax,error1,erro12,error2,errsum,ertest,oflow,
        resa,resabs,reseps,sign,temp,uflow;
      Real[3] res3la_;
      Real[52] rlist2_;
      int i,id,ierro,ind1,ind2,ip1,iroff1,iroff2,iroff3,j,
        jlow,jupbnd,k,ksgn,ktmin,levcur,levmax,maxerr,
        nint,nintp1,npts,nres,nrmax,numrl2;
      bool extrap,noext;
//
//
      auto alist = dimension(alist_, limit);
      auto blist = dimension(blist_, limit);
      auto elist = dimension(elist_, limit);
      auto iord = dimension(iord_, limit);
      auto level = dimension(level_, limit);
      auto ndin = dimension(ndin_, npts2);
      auto points = dimension(points_, npts2);
      auto pts = dimension(pts_, npts2);
      auto res3la = dimension(res3la_.ptr, 3);
      auto rlist = dimension(rlist_, limit);
      auto rlist2 = dimension(rlist2_.ptr, 52);
//
//            the dimension of rlist2 is determined by the value of
//            limexp in subroutine epsalg (rlist2 should be of dimension
//            (limexp+2) at least).
//
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
//           rlist2    - array of dimension at least limexp+2
//                       containing the part of the epsilon table which
//                       is still needed for further computations
//           elist(i)  - error estimate applying to rlist(i)
//           maxerr    - pointer to the interval with largest error
//                       estimate
//           errmax    - elist(maxerr)
//           erlast    - error on the interval currently subdivided
//                       (before that subdivision has taken place)
//           area      - sum of the integrals over the subintervals
//           errsum    - sum of the errors over the subintervals
//           errbnd    - requested accuracy max(epsabs,epsrel*
//                       abs(result))
//           *****1    - variable for the left subinterval
//           *****2    - variable for the right subinterval
//           last      - index for subdivision
//           nres      - number of calls to the extrapolation routine
//           numrl2    - number of elements in rlist2. if an appropriate
//                       approximation to the compounded integral has
//                       been obtained, it is put in rlist2(numrl2) after
//                       numrl2 has been increased by one.
//           erlarg    - sum of the errors over the intervals larger
//                       than the smallest interval considered up to now
//           extrap    - logical variable denoting that the routine
//                       is attempting to perform extrapolation. i.e.
//                       before subdividing the smallest interval we
//                       try to decrease the value of erlarg.
//           noext     - logical variable denoting that extrapolation is
//                       no longer allowed (true-value)
//
//            machine dependent constants
//            ---------------------------
//
//           epmach is the largest relative spacing.
//           uflow is the smallest positive magnitude.
//           oflow is the largest positive magnitude.
//
//***first executable statement  dqagpe
      epmach = Real.epsilon;
//
//            test on validity of parameters
//            -----------------------------
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
      level[1] = 0;
      npts = npts2-2;
      if(npts2 < 2 || limit <= npts || (epsabs <= 0.0 && 
        epsrel < max(0.5e2*epmach,0.5e-28))) ier = 6;
      if(ier == 6) goto l999;
//
//            if any break points are provided, sort them into an
//            ascending sequence.
//
      sign = 1.0;
      if(a > b) sign = -1.0;
      pts[1] = min(a,b);
      if(npts == 0) goto l15;
      for (i=1; i<=npts; i++) { //do 10 i = 1,npts
        pts[i+1] = points[i];
 l10: ;}
 l15: pts[npts+2] = max(a,b);
      nint = npts+1;
      a1 = pts[1];
      if(npts == 0) goto l40;
      nintp1 = nint+1;
      for (i=1; i<=nint; i++) { //do 20 i = 1,nint
        ip1 = i+1;
        for (j=ip1; j<=nintp1; j++) { //do 20 j = ip1,nintp1
          if(pts[i] <= pts[j]) goto l20;
          temp = pts[i];
          pts[i] = pts[j];
          pts[j] = temp;
 l20: ;}}
      if(pts[1] != min(a,b) || pts[nintp1] != max(a,b)) ier = 6;
      if(ier == 6) goto l999;
//
//            compute first integral and error approximations.
//            ------------------------------------------------
//
 l40: resabs = 0.0;
      for (i=1; i<=nint; i++) { //do 50 i = 1,nint
        b1 = pts[i+1];
        qk21!(Real,Func)(f,a1,b1,area1,error1,defabs,resa);
        abserr = abserr+error1;
        result = result+area1;
        ndin[i] = 0;
        if(error1 == resa && error1 != 0.0) ndin[i] = 1;
        resabs = resabs+defabs;
        level[i] = 0;
        elist[i] = error1;
        alist[i] = a1;
        blist[i] = b1;
        rlist[i] = area1;
        iord[i] = i;
        a1 = b1;
 l50: ;}
      errsum = 0.0;
      for (i=1; i<=nint; i++) { //do 55 i = 1,nint
        if(ndin[i] == 1) elist[i] = abserr;
        errsum = errsum+elist[i];
 l55: ;}
//
//           test on accuracy.
//
      last = nint;
      neval = 21*nint;
      dres = fabs(result);
      errbnd = max(epsabs,epsrel*dres);
      if(abserr <= 0.1e3*epmach*resabs && abserr > errbnd) ier = 2;
      if(nint == 1) goto l80;
      for (i=1; i<=npts; i++) { //do 70 i = 1,npts
        jlow = i+1;
        ind1 = iord[i];
        for (j=jlow; j<=nint; j++) { //do 60 j = jlow,nint
          ind2 = iord[j];
          if(elist[ind1] > elist[ind2]) goto l60;
          ind1 = ind2;
          k = j;
 l60:   ;}
        if(ind1 == iord[i]) goto l70;
        iord[k] = iord[i];
        iord[i] = ind1;
 l70: ;}
      if(limit < npts2) ier = 1;
 l80: if(ier != 0 || abserr <= errbnd) goto l210;
//
//           initialization
//           --------------
//
      rlist2[1] = result;
      maxerr = iord[1];
      errmax = elist[maxerr];
      area = result;
      nrmax = 1;
      nres = 0;
      numrl2 = 1;
      ktmin = 0;
      extrap = false;
      noext = false;
      erlarg = errsum;
      ertest = errbnd;
      levmax = 1;
      iroff1 = 0;
      iroff2 = 0;
      iroff3 = 0;
      ierro = 0;
      uflow = Real.min_normal;
      oflow = Real.max;
      abserr = oflow;
      ksgn = -1;
      if(dres >= (0.1e1-0.5e2*epmach)*resabs) ksgn = 1;
//
//           main do-loop
//           ------------
//
      for (last=npts2; last<=limit; last++) { //do 160 last = npts2,limit
//
//           bisect the subinterval with the nrmax-th largest error
//           estimate.
//
        levcur = level[maxerr]+1;
        a1 = alist[maxerr];
        b1 = 0.5*(alist[maxerr]+blist[maxerr]);
        a2 = b1;
        b2 = blist[maxerr];
        erlast = errmax;
        qk21!(Real,Func)(f,a1,b1,area1,error1,resa,defab1);
        qk21!(Real,Func)(f,a2,b2,area2,error2,resa,defab2);
//
//           improve previous approximations to integral
//           and error and test for accuracy.
//
        neval = neval+42;
        area12 = area1+area2;
        erro12 = error1+error2;
        errsum = errsum+erro12-errmax;
        area = area+area12-rlist[maxerr];
        if(defab1 == error1 || defab2 == error2) goto l95;
        if(fabs(rlist[maxerr]-area12) > 0.1e-4*fabs(area12)
         || erro12 < 0.99*errmax) goto l90;
        if(extrap) iroff2 = iroff2+1;
        if(!extrap) iroff1 = iroff1+1;
 l90:   if(last > 10 && erro12 > errmax) iroff3 = iroff3+1;
 l95:   level[maxerr] = levcur;
        level[last] = levcur;
        rlist[maxerr] = area1;
        rlist[last] = area2;
        errbnd = max(epsabs,epsrel*fabs(area));
//
//           test for roundoff error and eventually set error flag.
//
        if(iroff1+iroff2 >= 10 || iroff3 >= 20) ier = 2;
        if(iroff2 >= 5) ierro = 3;
//
//           set error flag in the case that the number of
//           subintervals equals limit.
//
        if(last == limit) ier = 1;
//
//           set error flag in the case of bad integrand behaviour
//           at a point of the integration range
//
        if(max(fabs(a1),fabs(b2)) <= (0.1e1+0.1e3*epmach)*
          (fabs(a2)+0.1e4*uflow)) ier = 4;
//
//           append the newly-created intervals to the list.
//
        if(error2 > error1) goto l100;
        alist[last] = a2;
        blist[maxerr] = b1;
        blist[last] = b2;
        elist[maxerr] = error1;
        elist[last] = error2;
        goto l110;
l100:   alist[maxerr] = a2;
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
l110:   qpsrt!Real(limit,last,maxerr,errmax,elist.ptr,iord.ptr,nrmax);
// ***jump out of do-loop
        if(errsum <= errbnd) goto l190;
// ***jump out of do-loop
        if(ier != 0) goto l170;
        if(noext) goto l160;
        erlarg = erlarg-erlast;
        if(levcur+1 <= levmax) erlarg = erlarg+erro12;
        if(extrap) goto l120;
//
//           test whether the interval to be bisected next is the
//           smallest interval.
//
        if(level[maxerr]+1 <= levmax) goto l160;
        extrap = true;
        nrmax = 2;
l120:   if(ierro == 3 || erlarg <= ertest) goto l140;
//
//           the smallest interval has the largest error.
//           before bisecting decrease the sum of the errors over
//           the larger intervals (erlarg) and perform extrapolation.
//
        id = nrmax;
        jupbnd = last;
        if(last > (2+limit/2)) jupbnd = limit+3-last;
        for (k=id; k<=jupbnd; k++) { //do 130 k = id,jupbnd
          maxerr = iord[nrmax];
          errmax = elist[maxerr];
// ***jump out of do-loop
          if(level[maxerr]+1 <= levmax) goto l160;
          nrmax = nrmax+1;
l130:   ;}
//
//           perform extrapolation.
//
l140:   numrl2 = numrl2+1;
        rlist2[numrl2] = area;
        if(numrl2 <= 2) goto l155;
        qelg!Real(numrl2,rlist2.ptr,reseps,abseps,res3la.ptr,nres);
        ktmin = ktmin+1;
        if(ktmin > 5 && abserr < 0.1e-2*errsum) ier = 5;
        if(abseps >= abserr) goto l150;
        ktmin = 0;
        abserr = abseps;
        result = reseps;
        correc = erlarg;
        ertest = max(epsabs,epsrel*fabs(reseps));
// ***jump out of do-loop
        if(abserr < ertest) goto l170;
//
//           prepare bisection of the smallest interval.
//
l150:   if(numrl2 == 1) noext = true;
        if(ier >= 5) goto l170;
l155:   maxerr = iord[1];
        errmax = elist[maxerr];
        nrmax = 1;
        extrap = false;
        levmax = levmax+1;
        erlarg = errsum;
l160: ;}
//
//           set the final result.
//           ---------------------
//
//
l170: if(abserr == oflow) goto l190;
      if((ier+ierro) == 0) goto l180;
      if(ierro == 3) abserr = abserr+correc;
      if(ier == 0) ier = 3;
      if(result != 0.0 && area != 0.0)goto l175;
      if(abserr > errsum)goto l190;
      if(area == 0.0) goto l210;
      goto l180;
l175: if(abserr/fabs(result) > errsum/fabs(area))goto l190;
//
//           test on divergence.
//
l180: if(ksgn == (-1) && max(fabs(result),fabs(area)) <= 
        resabs*0.1e-1) goto l210;
      if(0.1e-1 > (result/area) || (result/area) > 0.1e3 || 
        errsum > fabs(area)) ier = 6;
      goto l210;
//
//           compute global integral sum.
//
l190: result = 0.0;
      for (k=1; k<=last; k++) { //do 200 k = 1,last
        result = result+rlist[k];
l200: ;}
      abserr = errsum;
l210: if(ier > 2) ier = ier-1;
      result = result*sign;
l999: return;
}

unittest
{
    alias qagpe!(float, float delegate(float)) fqagpe;
    alias qagpe!(double, double delegate(double)) dqagpe;
    alias qagpe!(double, double function(double)) dfqagpe;
    alias qagpe!(real, real delegate(real)) rqagpe;
}
