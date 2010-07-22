/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.quadpack.qwgtf;


import std.math: cos, sin;




///
Real qwgtf(Real)(Real x, Real omega, Real p2, Real p3, Real p4, int integr)
{
//***begin prologue  dqwgtf
//***refer to   dqk15w
//***routines called  (none)
//***revision date 810101   (yymmdd)
//***keywords  cos or sin in weight function
//***author  piessens,robert, appl. math. & progr. div. - k.u.leuven
//           de doncker,elise,appl. math. * progr. div. - k.u.leuven
//***end prologue  dqwgtf
//
//***first executable statement  dqwgtf
      switch (integr)
      {
        case 1:  return cos(omega*x);
        case 2:  return sin(omega*x);
        default: return Real.nan;
      }
}


unittest
{
    alias qwgtf!float fqwgtf;
    alias qwgtf!double dqwgtf;
    alias qwgtf!real rqwgtf;
}
