/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qawfe;


import std.algorithm: max;
import std.conv;
import std.math;

import scid.common.fortran;
import scid.ports.quadpack.qagie;
import scid.ports.quadpack.qawoe;
import scid.ports.quadpack.qelg;




///
void qawfe(Real, Func)(Func f, Real a, Real omega, int integr, Real epsabs,
    int limlst, int limit, int maxp1, out Real result, out Real abserr,
    out int neval, out int ier, Real* rslst_, Real* erlst_, int* ierlst_,
    out int lst, Real* alist_, Real* blist_, Real* rlist_, Real* elist_,
    int* iord_, int* nnlog_, Real* chebmo_)
{
//***begin prologue  dqawfe
//***date written   800101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a3a1
//***keywords  automatic integrator, special-purpose,
//             fourier integrals,
//             integration between zeros with dqawoe,
//             convergence acceleration with dqelg
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           dedoncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  the routine calculates an approximation result to a
//            given fourier integal
//            i = integral of f(x)*w(x) over (a,infinity)
//            where w(x)=cos(omega*x) or w(x)=sin(omega*x),
//            hopefully satisfying following claim for accuracy
//            abs(i-result).le.epsabs.
//***description
//
//        computation of fourier integrals
//        standard fortran subroutine
//        double precision version
//
//        parameters
//         on entry
//            f      - double precision
//                     function subprogram defining the integrand
//                     function f(x). the actual name for f needs to
//                     be declared e x t e r n a l in the driver program.
//
//            a      - double precision
//                     lower limit of integration
//
//            omega  - double precision
//                     parameter in the weight function
//
//            integr - integer
//                     indicates which weight function is used
//                     integr = 1      w(x) = cos(omega*x)
//                     integr = 2      w(x) = sin(omega*x)
//                     if integr.ne.1.and.integr.ne.2, the routine will
//                     end with ier = 6.
//
//            epsabs - double precision
//                     absolute accuracy requested, epsabs.gt.0
//                     if epsabs.le.0, the routine will end with ier = 6.
//
//            limlst - integer
//                     limlst gives an upper bound on the number of
//                     cycles, limlst.ge.1.
//                     if limlst.lt.3, the routine will end with ier = 6.
//
//            limit  - integer
//                     gives an upper bound on the number of subintervals
//                     allowed in the partition of each cycle, limit.ge.1
//                     each cycle, limit.ge.1.
//
//            maxp1  - integer
//                     gives an upper bound on the number of
//                     chebyshev moments which can be stored, i.e.
//                     for the intervals of lengths abs(b-a)*2**(-l),
//                     l=0,1, ..., maxp1-2, maxp1.ge.1
//
//         on return
//            result - double precision
//                     approximation to the integral x
//
//            abserr - double precision
//                     estimate of the modulus of the absolute error,
//                     which should equal or exceed abs(i-result)
//
//            neval  - integer
//                     number of integrand evaluations
//
//            ier    - ier = 0 normal and reliable termination of
//                             the routine. it is assumed that the
//                             requested accuracy has been achieved.
//                     ier.gt.0 abnormal termination of the routine. the
//                             estimates for integral and error are less
//                             reliable. it is assumed that the requested
//                             accuracy has not been achieved.
//            error messages
//                    if omega.ne.0
//                     ier = 1 maximum number of  cycles  allowed
//                             has been achieved., i.e. of subintervals
//                             (a+(k-1)c,a+kc) where
//                             c = (2*int(abs(omega))+1)*pi/abs(omega),
//                             for k = 1, 2, ..., lst.
//                             one can allow more cycles by increasing
//                             the value of limlst (and taking the
//                             according dimension adjustments into
//                             account).
//                             examine the array iwork which contains
//                             the error flags on the cycles, in order to
//                             look for eventual local integration
//                             difficulties. if the position of a local
//                             difficulty can be determined (e.g.
//                             singularity, discontinuity within the
//                             interval) one will probably gain from
//                             splitting up the interval at this point
//                             and calling appropriate integrators on
//                             the subranges.
//                         = 4 the extrapolation table constructed for
//                             convergence acceleration of the series
//                             formed by the integral contributions over
//                             the cycles, does not converge to within
//                             the requested accuracy. as in the case of
//                             ier = 1, it is advised to examine the
//                             array iwork which contains the error
//                             flags on the cycles.
//                         = 6 the input is invalid because
//                             (integr.ne.1 and integr.ne.2) or
//                              epsabs.le.0 or limlst.lt.3.
//                              result, abserr, neval, lst are set
//                              to zero.
//                         = 7 bad integrand behaviour occurs within one
//                             or more of the cycles. location and type
//                             of the difficulty involved can be
//                             determined from the vector ierlst. here
//                             lst is the number of cycles actually
//                             needed (see below).
//                             ierlst(k) = 1 the maximum number of
//                                           subdivisions (= limit) has
//                                           been achieved on the k th
//                                           cycle.
//                                       = 2 occurrence of roundoff error
//                                           is detected and prevents the
//                                           tolerance imposed on the
//                                           k th cycle, from being
//                                           achieved.
//                                       = 3 extremely bad integrand
//                                           behaviour occurs at some
//                                           points of the k th cycle.
//                                       = 4 the integration procedure
//                                           over the k th cycle does
//                                           not converge (to within the
//                                           required accuracy) due to
//                                           roundoff in the
//                                           extrapolation procedure
//                                           invoked on this cycle. it
//                                           is assumed that the result
//                                           on this interval is the
//                                           best which can be obtained.
//                                       = 5 the integral over the k th
//                                           cycle is probably divergent
//                                           or slowly convergent. it
//                                           must be noted that
//                                           divergence can occur with
//                                           any other value of
//                                           ierlst(k).
//                    if omega = 0 and integr = 1,
//                    the integral is calculated by means of dqagie
//                    and ier = ierlst(1) (with meaning as described
//                    for ierlst(k), k = 1).
//
//            rslst  - double precision
//                     vector of dimension at least limlst
//                     rslst(k) contains the integral contribution
//                     over the interval (a+(k-1)c,a+kc) where
//                     c = (2*int(abs(omega))+1)*pi/abs(omega),
//                     k = 1, 2, ..., lst.
//                     note that, if omega = 0, rslst(1) contains
//                     the value of the integral over (a,infinity).
//
//            erlst  - double precision
//                     vector of dimension at least limlst
//                     erlst(k) contains the error estimate corresponding
//                     with rslst(k).
//
//            ierlst - integer
//                     vector of dimension at least limlst
//                     ierlst(k) contains the error flag corresponding
//                     with rslst(k). for the meaning of the local error
//                     flags see description of output parameter ier.
//
//            lst    - integer
//                     number of subintervals needed for the integration
//                     if omega = 0 then lst is set to 1.
//
//            alist, blist, rlist, elist - double precision
//                     vector of dimension at least limit,
//
//            iord, nnlog - integer
//                     vector of dimension at least limit, providing
//                     space for the quantities needed in the subdivision
//                     process of each cycle
//
//            chebmo - double precision
//                     array of dimension at least (maxp1,25), providing
//                     space for the chebyshev moments needed within the
//                     cycles
//
//***references  (none)
//***routines called  d1mach,dqagie,dqawoe,dqelg
//***end prologue  dqawfe
//
      Real abseps,correc,cycle,
        c1,c2,dl,dla,drl,ep,eps,epsa,
        errsum,fact,p,p1,reseps,
        uflow;
      Real[52] psum_;
      Real[3] res3la_;
      int ktmin,l,last,ll,
          momcom,nev,nres,numrl2;
//
      auto alist = dimension(alist_, limit);
      auto blist = dimension(blist_, limit);
      auto chebmo = dimension(chebmo_, maxp1, 25);
      auto elist = dimension(elist_, limit);
      auto erlst = dimension(erlst_, limlst);
      auto ierlst = dimension(ierlst_, limlst);
      auto iord = dimension(iord_, limit);
      auto nnlog = dimension(nnlog_, limit);
      auto psum = dimension(psum_.ptr, 52);
      auto res3la = dimension(res3la_.ptr, 3);
      auto rlist = dimension(rlist_, limit);
      auto rslst = dimension(rslst_, limlst);
//
//
//            the dimension of  psum  is determined by the value of
//            limexp in subroutine dqelg (psum must be of dimension
//            (limexp+2) at least).
//
//           list of major variables
//           -----------------------
//
//           c1, c2    - end points of subinterval (of length cycle)
//           cycle     - (2*int(abs(omega))+1)*pi/abs(omega)
//           psum      - vector of dimension at least (limexp+2)
//                       (see routine dqelg)
//                       psum contains the part of the epsilon table
//                       which is still needed for further computations.
//                       each element of psum is a partial sum of the
//                       series which should sum to the value of the
//                       integral.
//           errsum    - sum of error estimates over the subintervals,
//                       calculated cumulatively
//           epsa      - absolute tolerance requested over current
//                       subinterval
//           chebmo    - array containing the modified chebyshev
//                       moments (see also routine dqc25f)
//
      p = 0.9;
//
//           test on validity of parameters
//           ------------------------------
//
//***first executable statement  dqawfe
      result = 0.0;
      abserr = 0.0;
      neval = 0;
      lst = 0;
      ier = 0;
      if((integr != 1 && integr != 2) || epsabs <= 0.0 || 
        limlst < 3) ier = 6;
      if(ier == 6) goto l999;
      if(omega != 0.0) goto l10;
//
//           integration by dqagie if omega is zero
//           --------------------------------------
//
      if(integr == 1) qagie!(Real,Func)(f,0.0,1,epsabs,0.0,limit,
        result,abserr,neval,ier,alist.ptr,blist.ptr,rlist.ptr,elist.ptr,iord.ptr,last);
      rslst[1] = result;
      erlst[1] = abserr;
      ierlst[1] = ier;
      lst = 1;
      goto l999;
//
//           initializations
//           ---------------
//
 l10: l = to!int(fabs(omega));
      dl = 2*l+1;
      cycle = dl*PI/fabs(omega);
      ier = 0;
      ktmin = 0;
      neval = 0;
      numrl2 = 0;
      nres = 0;
      c1 = a;
      c2 = cycle+a;
      p1 = 0.1e1-p;
      uflow = Real.min_normal;
      eps = epsabs;
      if(epsabs > uflow/p1) eps = epsabs*p1;
      ep = eps;
      fact = 0.1e1;
      correc = 0.0;
      abserr = 0.0;
      errsum = 0.0;
//
//           main do-loop
//           ------------
//
      for (lst=1; lst<=limlst; lst++) { //do 50 lst = 1,limlst
//
//           integrate over current subinterval.
//
        dla = lst;
        epsa = eps*fact;
        qawoe!(Real,Func)(f,c1,c2,omega,integr,epsa,0.0,limit,lst,maxp1,
          rslst[lst],erlst[lst],nev,ierlst[lst],last,alist.ptr,blist.ptr,rlist.ptr,
          elist.ptr,iord.ptr,nnlog.ptr,momcom,chebmo.ptr);
        neval = neval+nev;
        fact = fact*p;
        errsum = errsum+erlst[lst];
        drl = 0.5e2*fabs(rslst[lst]);
//
//           test on accuracy with partial sum
//
        if((errsum+drl) <= epsabs && lst >= 6) goto l80;
        correc = max(correc,erlst[lst]);
        if(ierlst[lst] != 0) eps = max(ep,correc*p1);
        if(ierlst[lst] != 0) ier = 7;
        if(ier == 7 && (errsum+drl) <= correc*0.1e2 && 
          lst > 5) goto l80;
        numrl2 = numrl2+1;
        if(lst > 1) goto l20;
        psum[1] = rslst[1];
        goto l40;
 l20:   psum[numrl2] = psum[ll]+rslst[lst];
        if(lst == 2) goto l40;
//
//           test on maximum number of subintervals
//
        if(lst == limlst) ier = 1;
//
//           perform new extrapolation
//
        qelg!Real(numrl2,psum.ptr,reseps,abseps,res3la.ptr,nres);
//
//           test whether extrapolated result is influenced by roundoff
//
        ktmin = ktmin+1;
        if(ktmin >= 15 && abserr <= 0.1e-2*(errsum+drl)) ier = 4;
        if(abseps > abserr && lst != 3) goto l30;
        abserr = abseps;
        result = reseps;
        ktmin = 0;
//
//           if ier is not 0, check whether direct result (partial sum)
//           or extrapolated result yields the best integral
//           approximation
//
        if((abserr+0.1e2*correc) <= epsabs || 
          (abserr <= epsabs && 0.1e2*correc >= epsabs)) goto l60;
 l30:   if(ier != 0 && ier != 7) goto l60;
 l40:   ll = numrl2;
        c1 = c2;
        c2 = c2+cycle;
 l50: ;}
//
//         set final result and error estimate
//         -----------------------------------
//
 l60: abserr = abserr+0.1e2*correc;
      if(ier == 0) goto l999;
      if(result != 0.0 && psum[numrl2] != 0.0) goto l70;
      if(abserr > errsum) goto l80;
      if(psum[numrl2] == 0.0) goto l999;
 l70: if(abserr/fabs(result) > (errsum+drl)/fabs(psum[numrl2]))
        goto l80;
      if(ier >= 1 && ier != 7) abserr = abserr+drl;
      goto l999;
 l80: result = psum[numrl2];
      abserr = errsum+drl;
l999: return;
}


unittest
{
    alias qawfe!(float, float delegate(float)) fqawfe;
    alias qawfe!(double, double delegate(double)) dqawfe;
    alias qawfe!(double, double function(double)) dfqawfe;
    alias qawfe!(real, real delegate(real)) rqawfe;
}
