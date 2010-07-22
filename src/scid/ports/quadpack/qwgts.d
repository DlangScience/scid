/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qwgts;


import std.math;


///
Real qwgts(Real)(Real x, Real a, Real b, Real alfa, Real beta, int integr)
{
//***begin prologue  dqwgts
//***refer to dqk15w
//***routines called  (none)
//***revision date  810101   (yymmdd)
//***keywords  weight function, algebraico-logarithmic
//             end-point singularities
//***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. & progr. div. - k.u.leuven
//***purpose  this function subprogram is used together with the
//            routine dqaws and defines the weight function.
//***end prologue  dqwgts
//
      Real bmx,xma;
//***first executable statement  dqwgts
      xma = x-a;
      bmx = b-x;
      Real temp = (xma^^alfa)*(bmx^^beta);
      switch (integr)
      {
        case 1: return temp;
        case 2: return temp * log(xma);
        case 3: return temp * log(bmx);
        case 4: return temp * log(xma)*log(bmx);
        default: assert(0);
      }
}


unittest
{
    alias qwgts!float fqwgts;
    alias qwgts!double dqwgts;
    alias qwgts!real rqwgts;
}
