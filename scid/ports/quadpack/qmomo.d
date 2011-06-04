/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qmomo;


import std.math;
import scid.common.fortran;



///
void qmomo(Real)(Real alfa, Real beta, Real* ri_, Real* rj_, Real* rg_,
    Real* rh_, int integr)
{
//***begin prologue  dqmomo
//***date written   820101   (yymmdd)
//***revision date  830518   (yymmdd)
//***category no.  h2a2a1,c3a2
//***keywords  modified chebyshev moments
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  this routine computes modified chebsyshev moments. the k-th
//            modified chebyshev moment is defined as the integral over
//            (-1,1) of w(x)*t(k,x), where t(k,x) is the chebyshev
//            polynomial of degree k.
//***description
//
//        modified chebyshev moments
//        standard fortran subroutine
//        double precision version
//
//        parameters
//           alfa   - double precision
//                    parameter in the weight function w(x), alfa.gt.(-1)
//
//           beta   - double precision
//                    parameter in the weight function w(x), beta.gt.(-1)
//
//           ri     - double precision
//                    vector of dimension 25
//                    ri(k) is the integral over (-1,1) of
//                    (1+x)**alfa*t(k-1,x), k = 1, ..., 25.
//
//           rj     - double precision
//                    vector of dimension 25
//                    rj(k) is the integral over (-1,1) of
//                    (1-x)**beta*t(k-1,x), k = 1, ..., 25.
//
//           rg     - double precision
//                    vector of dimension 25
//                    rg(k) is the integral over (-1,1) of
//                    (1+x)**alfa*log((1+x)/2)*t(k-1,x), k = 1, ..., 25.
//
//           rh     - double precision
//                    vector of dimension 25
//                    rh(k) is the integral over (-1,1) of
//                    (1-x)**beta*log((1-x)/2)*t(k-1,x), k = 1, ..., 25.
//
//           integr - integer
//                    input parameter indicating the modified
//                    moments to be computed
//                    integr = 1 compute ri, rj
//                           = 2 compute ri, rj, rg
//                           = 3 compute ri, rj, rh
//                           = 4 compute ri, rj, rg, rh
//
//***references  (none)
//***routines called  (none)
//***end prologue  dqmomo
//
      Real alfp1,alfp2,an,anm1,betp1,betp2,ralf,rbet;
      int i,im1;
//
      auto rg = dimension(rg_, 25);
      auto rh = dimension(rh_, 25);
      auto ri = dimension(ri_, 25);
      auto rj = dimension(rj_, 25);
//
//
//***first executable statement  dqmomo
      alfp1 = alfa+0.1e1;
      betp1 = beta+0.1e1;
      alfp2 = alfa+0.2e1;
      betp2 = beta+0.2e1;
      ralf = (cast(Real)0.2e1)^^alfp1;
      rbet = (cast(Real)0.2e1)^^betp1;
//
//           compute ri, rj using a forward recurrence relation.
//
      ri[1] = ralf/alfp1;
      rj[1] = rbet/betp1;
      ri[2] = ri[1]*alfa/alfp2;
      rj[2] = rj[1]*beta/betp2;
      an = 0.2e1;
      anm1 = 0.1e1;
      for (i=3; i<=25; i++) { //do 20 i=3,25
        ri[i] = -(ralf+an*(an-alfp2)*ri[i-1])/(anm1*(an+alfp1));
        rj[i] = -(rbet+an*(an-betp2)*rj[i-1])/(anm1*(an+betp1));
        anm1 = an;
        an = an+0.1e1;
 l20: ;}
      if(integr == 1) goto l70;
      if(integr == 3) goto l40;
//
//           compute rg using a forward recurrence relation.
//
      rg[1] = -ri[1]/alfp1;
      rg[2] = -(ralf+ralf)/(alfp2*alfp2)-rg[1];
      an = 0.2e1;
      anm1 = 0.1e1;
      im1 = 2;
      for (i=3; i<=25; i++) { //do 30 i=3,25
        rg[i] = -(an*(an-alfp2)*rg[im1]-an*ri[im1]+anm1*ri[i])/
          (anm1*(an+alfp1));
        anm1 = an;
        an = an+0.1e1;
        im1 = i;
 l30: ;}
      if(integr == 2) goto l70;
//
//           compute rh using a forward recurrence relation.
//
 l40: rh[1] = -rj[1]/betp1;
      rh[2] = -(rbet+rbet)/(betp2*betp2)-rh[1];
      an = 0.2e1;
      anm1 = 0.1e1;
      im1 = 2;
      for (i=3; i<=25; i++) { //do 50 i=3,25
        rh[i] = -(an*(an-betp2)*rh[im1]-an*rj[im1]+
        anm1*rj[i])/(anm1*(an+betp1));
        anm1 = an;
        an = an+0.1e1;
        im1 = i;
 l50: ;}
      for (i=2; i<=25; i+=2) { //do 60 i=2,25,2
        rh[i] = -rh[i];
 l60: ;}
 l70: for (i=2; i<=25; i+=2) { //do 80 i=2,25,2
        rj[i] = -rj[i];
 l80: ;}
 l90: return;
}

unittest
{
    alias qmomo!float fqmomo;
    alias qmomo!double dqmomo;
    alias qmomo!real rqmomo;
}
