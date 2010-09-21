/** Special functions.

    The functions in this module are ported (with permission) to D from
    $(LINK2 http://www.moshier.net/,Stephen L. Moshier)'s
    $(LINK2 http://www.netlib.org/cephes/,Cephes) library.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.functions;



import std.math;

import scid.core.testing;




/*  Evaluate Chebyshev series.

    Coefficients are stored in reverse order, i.e. the
    zeroth-order term is the last in the array.
*/
private real chebEval(real x, const real[] coefficients)
{
    real
        b0 = coefficients[0],
        b1 = 0,
        b2;

    foreach (c; coefficients[1 .. $])
    {
        b2 = b1;
        b1 = b0;
        b0 = x * b1 - b2 + c;
    }

    return 0.5 * (b0 - b2);
}




/** Modified Bessel function of the first kind, order 1. */
real besselI1(real x)
{
    // Chebyshev coefficients for exp(-x) I1(x) / x
    // in the interval [0,8].
    immutable real[29] chebCoeffsA = [
        2.77791411276104639959e-18L,
       -2.11142121435816608115e-17L,
        1.55363195773620046921e-16L,
       -1.10559694773538630805e-15L,
        7.60068429473540693410e-15L,
       -5.04218550472791168711e-14L,
        3.22379336594557470981e-13L,
       -1.98397439776494371520e-12L,
        1.17361862988909016308e-11L,
       -6.66348972350202774223e-11L,
        3.62559028155211703701e-10L,
       -1.88724975172282928790e-9L,
        9.38153738649577178388e-9L,
       -4.44505912879632808065e-8L,
        2.00329475355213526229e-7L,
       -8.56872026469545474066e-7L,
        3.47025130813767847674e-6L,
       -1.32731636560394358279e-5L,
        4.78156510755005422638e-5L,
       -1.61760815825896745588e-4L,
        5.12285956168575772895e-4L,
       -1.51357245063125314899e-3L,
        4.15642294431288815669e-3L,
       -1.05640848946261981558e-2L,
        2.47264490306265168283e-2L,
       -5.29459812080949914269e-2L,
        1.02643658689847095384e-1L,
       -1.76416518357834055153e-1L,
        2.52587186443633654823e-1L  ];

    // Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
    // in the inverted interval [8,infinity].
    immutable real[25] chebCoeffsB = [
        7.51729631084210481353e-18L,
        4.41434832307170791151e-18L,
       -4.65030536848935832153e-17L,
       -3.20952592199342395980e-17L,
        2.96262899764595013876e-16L,
        3.30820231092092828324e-16L,
       -1.88035477551078244854e-15L,
       -3.81440307243700780478e-15L,
        1.04202769841288027642e-14L,
        4.27244001671195135429e-14L,
       -2.10154184277266431302e-14L,
       -4.08355111109219731823e-13L,
       -7.19855177624590851209e-13L,
        2.03562854414708950722e-12L,
        1.41258074366137813316e-11L,
        3.25260358301548823856e-11L,
       -1.89749581235054123450e-11L,
       -5.58974346219658380687e-10L,
       -3.83538038596423702205e-9L,
       -2.63146884688951950684e-8L,
       -2.51223623787020892529e-7L,
       -3.88256480887769039346e-6L,
       -1.10588938762623716291e-4L,
       -9.76109749136146840777e-3L,
        7.78576235018280120474e-1L  ];


    auto z = fabs(x);
    if (z <= 8)
    {
        immutable y = z/2 - 2;
        z = chebEval(y, chebCoeffsA[]) * z * exp(z);
    }
    else
    {
        z = exp(z) * chebEval(32/z - 2, chebCoeffsB[]) / sqrt(z);
    }
    if (x < 0) z = -z;

    return z;
}


unittest
{
    check (approxEqual(besselI1( 1), 5.6515910399248502721e-1, 1.2e-16));
    check (approxEqual(besselI1(10), 2.6709883037012546543e3,  1.2e-16));
}




/** Modified Bessel function of the second kind, order 1. */
real besselK1(real x)
in { assert(x > 0, "Can only calculate K1 for positive argument"); }
body {
    // Chebyshev coefficients for x(K1(x) - log(x/2) I1(x))
    // in the interval [0,2].
    immutable real[11] chebCoeffsA = [
       -7.02386347938628759343e-18L,
       -2.42744985051936593393e-15L,
       -6.66690169419932900609e-13L,
       -1.41148839263352776110e-10L,
       -2.21338763073472585583e-8L,
       -2.43340614156596823496e-6L,
       -1.73028895751305206302e-4L,
       -6.97572385963986435018e-3L,
       -1.22611180822657148235e-1L,
       -3.53155960776544875667e-1L,
        1.52530022733894777053e0L  ];

    // Chebyshev coefficients for exp(x) sqrt(x) K1(x)
    // in the interval [2,infinity].
    immutable real[25] chebCoeffsB = [
       -5.75674448366501715755e-18L,
        1.79405087314755922667e-17L,
       -5.68946255844285935196e-17L,
        1.83809354436663880070e-16L,
       -6.05704724837331885336e-16L,
        2.03870316562433424052e-15L,
       -7.01983709041831346144e-15L,
        2.47715442448130437068e-14L,
       -8.97670518232499435011e-14L,
        3.34841966607842919884e-13L,
       -1.28917396095102890680e-12L,
        5.13963967348173025100e-12L,
       -2.12996783842756842877e-11L,
        9.21831518760500529508e-11L,
       -4.19035475934189648750e-10L,
        2.01504975519703286596e-9L,
       -1.03457624656780970260e-8L,
        5.74108412545004946722e-8L,
       -3.50196060308781257119e-7L,
        2.40648494783721712015e-6L,
       -1.93619797416608296024e-5L,
        1.95215518471351631108e-4L,
       -2.85781685962277938680e-3L,
        1.03923736576817238437e-1L,
        2.72062619048444266945e0L   ];


    real z = 0.5 * x;
    if (x <= 2.0)
    {
        real y = x^^2 - 2;
        y = log(z) * besselI1(x) + chebEval(y, chebCoeffsA[]) / x;
        return y;
    }
    
    return exp(-x) * chebEval(8.0/x - 2.0, chebCoeffsB[]) / sqrt(x);
}


unittest
{
    check (approxEqual(besselK1( 1), 6.0190723019723457474e-1, 8.9e-17));
    check (approxEqual(besselK1(10), 1.8648773453825584597e-5, 8.9e-17));
}

