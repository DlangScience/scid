/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.ports.napack.addchg;


import std.math;




//
//      ________________________________________________________
//     |                                                        |
//     |       ADD CHG TO NEW AND EVALUATE DIF AND SIZE         |
//     |________________________________________________________|
///
void addchg(Real)(ref Real dif, ref Real size, Real[] new_, Real[] chg, int n)
{
    dif = 0.0;
    size = 0.0;
    foreach (i; 0 .. n)
    {
        dif += fabs(chg[i]);
        new_[i] += chg[i];
        size += fabs(new_[i]);
    }
}


unittest
{
    alias addchg!float faddchg;
    alias addchg!double daddchg;
    alias addchg!real raddchg;
}
