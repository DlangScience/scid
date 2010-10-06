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



version (LittleEndian)
{
    static if (real.sizeof == 12)
    {
        version = LEReal80;
        private enum short PAD = 0;
    }
}




/*  Stirling's formula for the gamma function. */
private real stirlingGamma(real x)
{
    auto w = 1.0L/x;

    // For large z, use rational coefficients from the analytical expansion.
    if (x > 1024.0L)
    {
        w = (((((
            6.97281375836585777429e-5L * w
          + 7.84039221720066627474e-4L) * w
          - 2.29472093621399176955e-4L) * w
          - 2.68132716049382716049e-3L) * w
          + 3.47222222222222222222e-3L) * w
          + 8.33333333333333333333e-2L) * w
          + 1.0L;
    }
    else
    {
        version (LEReal80)
        {
            immutable short[9*6] coeffsBytes = cast(short[]) [
                0xa1d5,0xaaaa,0xaaaa,0xaaaa,0x3ffb, PAD,
                0xc3c9,0x906e,0x38e3,0xe38e,0x3ff6, PAD,
                0x3a1c,0x5ac8,0x3478,0xafb9,0xbff6, PAD,
                0xbef3,0x7023,0x6a08,0xf09e,0xbff2, PAD,
                0x30b7,0x1a21,0x98b2,0xcd87,0x3ff4, PAD,
                0x5704,0x1a39,0xb11d,0x9293,0x3ff1, PAD,
                0xba6f,0x7c59,0x5e47,0x9bfb,0xbff4, PAD,
                0xc395,0x0295,0x4443,0xc64b,0xbfef, PAD,
                0x6ede,0x69f7,0x54e3,0xbb5d,0x3ff4, PAD ];
            auto coeffs = cast(immutable(real[])) coeffsBytes[];
        }
        else
        {
            immutable real[9] coeffs = [
                8.333333333333331800504e-2L,
                3.472222222230075327854e-3L,
               -2.681327161876304418288e-3L,
               -2.294719747873185405699e-4L,
                7.840334842744753003862e-4L,
                6.989332260623193171870e-5L,
               -5.950237554056330156018e-4L,
               -2.363848809501759061727e-5L,
                7.147391378143610789273e-4L ];
        }

        w = 1.0L + w * poly(w, coeffs[]);
    }

    auto y = exp(x);

    enum maxStirling = 1024.0L;
    if (x > maxStirling)
    {
        // Avoid overflow in pow()
        auto v = x^^(0.5L*x - 0.25L);
        y = v * (v/y);
    }
    else
    {
        y = x^^(x - 0.5L) / y;
    }

    enum sqrt2Pi = sqrt(2*PI);
    return sqrt2Pi * y * w;
}


unittest
{
    // The Cephes docs claims this accuracy:
    enum eps = 9.44e-21L;
    check (approxEqual(stirlingGamma(  13.0L), 4.79001600000000000000e8L,    eps));

    // ...but we don't always reach it:
    check (approxEqual(stirlingGamma( 100.0L), 9.33262154439441526817e155L,  1e-16));
    check (approxEqual(stirlingGamma(1000.0L), 4.02387260077093773544e2564L, 1e-15L));
}




/** Gamma function. */
real gamma(real x)
{

    
    // Handle NaN and infinities
    if (isNaN(x)) return real.nan;
    if (isInfinity(x))
    {
        if (x > 0) return real.infinity;
        else return real.nan;
    }

    int sign = 1;
    auto xAbs = fabs(x);

    // For large numbers, use Stirling's formula
    if (xAbs > 13.0L)
    {
        enum maxGamma = 1755.455L;
        if (xAbs > maxGamma) return sign * real.infinity;

        real z;

        // Reflection formula for negative numbers
        if (x < 0)
        {
            auto xInt = floor(xAbs);

            // Check whether x is a negative integer
            if (xInt == xAbs) return real.nan;

            // If p is odd, the gamma function is negative
            int i = cast(int) xInt;
            if ((i & 1) == 0) sign = -1;

            z = xAbs - xInt;
            if (z > 0.5L)
            {
                xInt += 1;
                z = xAbs - xInt;
            }

            z = xAbs * sin(PI * z);
            z = fabs(z) * stirlingGamma(xAbs);
            if (z <= PI/real.max)  return sign * real.infinity;

            z = PI/z;
        }
        else
        {
            z = stirlingGamma(x);
        }

        return sign * z;
    }

    real z = 1;

    // For x in the ranges [-13,-0.03125] and [3,13], use recursion formula
    while (x >= 3)
    {
        x -= 1;
        z *= x;
    }

    while (x < -0.03125L)
    {
        z /= x;
        x += 1;
    }

    // Small arguments
    if (x <= 0.03125L)
    {
        if (x == 0)  return real.nan;
        else if (x < 0)
        {
            version (LEReal80)
            {
                immutable short[9*6] coeffsSmallNegShort = cast(short[]) [
                    0x0000,0x0000,0x0000,0x8000,0xbfff, PAD,
                    0xc7aa,0x7db0,0x67e3,0x93c4,0x3ffe, PAD,
                    0x5e26,0x57d1,0xa013,0xa7e7,0x3ffe, PAD,
                    0x7f64,0x1234,0xf47d,0xac0a,0xbffa, PAD,
                    0x7a5b,0xd76d,0x1905,0xaa89,0xbffc, PAD,
                    0x783f,0x41dd,0x87d1,0xacd7,0xbffa, PAD,
                    0x2ca1,0x18f0,0x386f,0x9da5,0x3ff8, PAD,
                    0x989b,0xdd68,0xc5f1,0xec9c,0x3ff7, PAD,
                    0x5dd1,0x02de,0xb9f7,0x948d,0x3ff5, PAD,
                    ];
                auto coeffsSmallNeg =
                    cast(immutable(real[])) coeffsSmallNegShort[];
            }
            else
            {
                immutable real[9] coeffsSmallNeg = [
                   -1.000000000000000000000e0L,
                    5.772156649015328608727e-1L,
                    6.558780715202536547116e-1L,
                   -4.200263503402112910504e-2L,
                   -1.665386113944413519335e-1L,
                   -4.219773343731191721664e-2L,
                    9.621911155035976733706e-3L,
                    7.220837261893170325704e-3L,
                    1.133374167243894382010e-3L,
                    ];
            }
            x = -x;
            return z / (x * poly(x, coeffsSmallNeg));
        }
        else
        {
            version (LEReal80)
            {
                immutable short[6*9] coeffsSmallPosShort = cast(short[]) [
                    0x0000,0x0000,0x0000,0x8000,0x3fff, PAD,
                    0xc7a9,0x7db0,0x67e3,0x93c4,0x3ffe, PAD,
                    0x7bf6,0x57d1,0xa013,0xa7e7,0xbffe, PAD,
                    0xf183,0x126b,0xf47d,0xac0a,0xbffa, PAD,
                    0x6b8d,0x7515,0x1905,0xaa89,0x3ffc, PAD,
                    0x10b0,0xec17,0x87dc,0xacd7,0xbffa, PAD,
                    0x9225,0xdfef,0xb0e9,0x9da5,0xbff8, PAD,
                    0xfe9a,0xceb4,0xc74e,0xec9a,0x3ff7, PAD,
                    0xbaeb,0xd6d3,0x25e5,0x9c7e,0xbff5, PAD,
                    ];
                auto coeffsSmallPos =
                    cast(immutable(real[])) coeffsSmallPosShort[];
        }
            else
            {
                immutable real[9] coeffsSmallPos = [
                    1.000000000000000000000E0L,
                    5.772156649015328608253E-1L,
                   -6.558780715202540684668E-1L,
                   -4.200263503403344054473E-2L,
                    1.665386113720805206758E-1L,
                   -4.219773360705915470089E-2L,
                   -9.622023360406271645744E-3L,
                    7.220599478036909672331E-3L,
                   -1.193945051381510095614E-3L,
                    ];
            }
            return z / (x * poly(x, coeffsSmallPos));
        }
        assert(0);
    }

    while (x < 2)
    {
        z /= x;
        x += 1;
    }

    if (x == 2)  return z;

    x -= 2;

    version (LEReal80)
    {
        static immutable short[8*6] coeffsNumerShort = cast(short[]) [
            0x0000,0x0000,0x0000,0x8000,0x3fff, PAD,
            0x29cf,0x19b3,0x16c8,0xd67a,0x3ffe, PAD,
            0x8d75,0x23af,0xc8e4,0xb9d4,0x3ffd, PAD,
            0x9549,0x8eb5,0x8c3a,0xe3f4,0x3ffb, PAD,
            0x7f43,0x5196,0xb166,0xc368,0x3ff9, PAD,
            0xbe6c,0x3757,0xc717,0x861b,0x3ff7, PAD,
            0xf5aa,0xe82f,0x335b,0xee2e,0x3ff3, PAD,
            0x434a,0x3f22,0x2bda,0xb0b2,0x3ff0, PAD,
            ];
        auto coeffsNumer = cast(immutable(real[])) coeffsNumerShort[];
        immutable short[9*6] coeffsDenomShort = cast(short[]) [
            0x0000,0x0000,0x0000,0x8000,0x3fff, PAD,
            0xe458,0x2ec7,0xfd57,0xd47c,0x3ffd, PAD,
            0x75ef,0x3ab7,0x4ad3,0xe5bc,0xbffc, PAD,
            0x3295,0x3698,0xd580,0xbdcd,0xbffa, PAD,
            0x0417,0x7989,0xd7bc,0xe338,0x3ff9, PAD,
            0x296e,0x7cb1,0x5dfd,0xd08f,0xbff4, PAD,
            0xbeed,0x1853,0xa691,0xa23d,0xbff5, PAD,
            0x334b,0xc2f0,0xa2dd,0xf60e,0x3ff2, PAD,
            0x5473,0x2de8,0x1268,0xea67,0xbfee, PAD,
            ];
        auto coeffsDenom = cast(immutable(real[])) coeffsDenomShort[];
    }
    else
    {
        immutable real[8] coeffsNumer = [
            1.000000000000000000009E0L,
            8.378004301573126728826E-1L,
            3.629515436640239168939E-1L,
            1.113062816019361559013E-1L,
            2.385363243461108252554E-2L,
            4.092666828394035500949E-3L,
            4.542931960608009155600E-4L,
            4.212760487471622013093E-5L,
            ];
        immutable real[9] coeffsDenom = [
            9.999999999999999999908E-1L,
            4.150160950588455434583E-1L,
           -2.243510905670329164562E-1L,
           -4.633887671244534213831E-2L,
            2.773706565840072979165E-2L,
           -7.955933682494738320586E-4L,
           -1.237799246653152231188E-3L,
            2.346584059160635244282E-4L,
           -1.397148517476170440917E-5L,
            ];
    }
    auto p = poly(x, coeffsNumer);
    auto q = poly(x, coeffsDenom);
    return z * p / q;
}


unittest
{
    // Tests for 13 < |x| < 40
    enum epsSmall = 3.6e-19L;
    check (approxEqual(gamma(-39.9L),  1.80391216710092201696e-47L, epsSmall));
    check (approxEqual(gamma(-24.6L), -7.75453165607905780208e-25L, epsSmall));
    check (approxEqual(gamma(-13.1L),  1.25800832362047759474e-09L, epsSmall));
    check (approxEqual(gamma(- 6.2L), -5.08874226047190522657e-03L, epsSmall));
    check (approxEqual(gamma(- 1e-3L),-1.00057820562935864799e003L, epsSmall));
    check (approxEqual(gamma(  1e-3L), 9.99423772484595466115e002L, epsSmall));
    check (gamma(1) == 1);
    check (gamma(2) == 1);
    check (gamma(3) == 2);
    check (gamma(4) == 6);
    check (approxEqual(gamma(  6.7L),  4.13407516765270695563e002L, epsSmall));
    check (approxEqual(gamma( 14.0L),  6.22702080000000000000e009L, epsSmall));
    check (approxEqual(gamma( 20.0L),  1.21645100408832000000e017L, epsSmall));
    check (approxEqual(gamma( 35.3L),  8.55223625497227026058e038L, 1e-17));

    // Tests for 40 < |x| < 1755
    enum epsLarge = 4.8e-18L;
    check (approxEqual(gamma(-1234.5L), -1.75011710987893708501e-3282L, epsLarge));
    check (approxEqual(gamma(- 123.4L),  3.95868650220624836363e-0206L, epsLarge));
    check (approxEqual(gamma(-  40.1L), -8.60443315835293399214e-0048L, epsLarge));
    check (approxEqual(gamma(   40.9L),  5.63575520473871710921e00047L, 1e-17));
    check (approxEqual(gamma(  246.8L),  2.81818959867730772329e00482L, 1e-16));
    check (approxEqual(gamma( 1754.9L),  9.37754511960745549608e04929L, 1e-15));
    

}




/*  Evaluate Chebyshev series.

    Coefficients are stored in reverse order, i.e. the
    zeroth-order term is the last in the array.
*/
private real chebEval(real x, const real[] coefficients) @safe pure nothrow
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




/** Modified Bessel function of the first kind, order 0. */
real besselI0(real x) pure
{
    // Chebyshev coefficients for exp(-x) I0(x)
    // in the interval [0,8].
    //
    // lim(x->0){ exp(-x) I0(x) } = 1.
    immutable real[30] chebCoeffsA = [
       -4.41534164647933937950e-18L,
        3.33079451882223809783e-17L,
       -2.43127984654795469359e-16L,
        1.71539128555513303061e-15L,
       -1.16853328779934516808e-14L,
        7.67618549860493561688e-14L,
       -4.85644678311192946090e-13L,
        2.95505266312963983461e-12L,
       -1.72682629144155570723e-11L,
        9.67580903537323691224e-11L,
       -5.18979560163526290666e-10L,
        2.65982372468238665035e-9L,
       -1.30002500998624804212e-8L,
        6.04699502254191894932e-8L,
       -2.67079385394061173391e-7L,
        1.11738753912010371815e-6L,
       -4.41673835845875056359e-6L,
        1.64484480707288970893e-5L,
       -5.75419501008210370398e-5L,
        1.88502885095841655729e-4L,
       -5.76375574538582365885e-4L,
        1.63947561694133579842e-3L,
       -4.32430999505057594430e-3L,
        1.05464603945949983183e-2L,
       -2.37374148058994688156e-2L,
        4.93052842396707084878e-2L,
       -9.49010970480476444210e-2L,
        1.71620901522208775349e-1L,
       -3.04682672343198398683e-1L,
        6.76795274409476084995e-1L  ];

    // Chebyshev coefficients for exp(-x) sqrt(x) I0(x)
    // in the inverted interval [8,infinity].
    //
    // lim(x->inf){ exp(-x) sqrt(x) I0(x) } = 1/sqrt(2pi).
    immutable real[25] chebCoeffsB = [
       -7.23318048787475395456e-18L,
       -4.83050448594418207126e-18L,
        4.46562142029675999901e-17L,
        3.46122286769746109310e-17L,
       -2.82762398051658348494e-16L,
       -3.42548561967721913462e-16L,
        1.77256013305652638360e-15L,
        3.81168066935262242075e-15L,
       -9.55484669882830764870e-15L,
       -4.15056934728722208663e-14L,
        1.54008621752140982691e-14L,
        3.85277838274214270114e-13L,
        7.18012445138366623367e-13L,
       -1.79417853150680611778e-12L,
       -1.32158118404477131188e-11L,
       -3.14991652796324136454e-11L,
        1.18891471078464383424e-11L,
        4.94060238822496958910e-10L,
        3.39623202570838634515e-9L,
        2.26666899049817806459e-8L,
        2.04891858946906374183e-7L,
        2.89137052083475648297e-6L,
        6.88975834691682398426e-5L,
        3.36911647825569408990e-3L,
        8.04490411014108831608e-1L  ];

    auto y = fabs(x);

    if (y < 8)
        return exp(y) * chebEval(y/2 - 2, chebCoeffsA);
    else
        return exp(y) * chebEval(32/y - 2, chebCoeffsB) / sqrt(y);
}


unittest
{
    check (approxEqual(besselI0( 1), 1.2660658777520083356e0, 8.2e-17));
    check (approxEqual(besselI0(10), 2.8157166284662544715e3, 8.2e-17));
    check (besselI0(- 1) == besselI0( 1));
    check (besselI0(-10) == besselI0(10));
}




/** Modified Bessel function of the first kind, order 1. */
real besselI1(real x) pure
{
    /* Chebyshev coefficients for exp(-x) I1(x) / x
     * in the interval [0,8].
     *
     * lim(x->0){ exp(-x) I1(x) / x } = 1/2.
     */
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

    /* Chebyshev coefficients for exp(-x) sqrt(x) I1(x)
     * in the inverted interval [8,infinity].
     *
     * lim(x->inf){ exp(-x) sqrt(x) I1(x) } = 1/sqrt(2pi).
     */
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
    check (besselI1(- 1) == -besselI1( 1));
    check (besselI1(-10) == -besselI1(10));
}




/** Modified Bessel function of the second kind, order 0. */
real besselK0(real x) pure
in { assert(x > 0, "Can only calculate K0 for positive argument"); }
body {
    /* Chebyshev coefficients for K0(x) + log(x/2) I0(x)
     * in the interval [0,2].  The odd order coefficients are all
     * zero; only the even order coefficients are listed.
     * 
     * lim(x->0){ K0(x) + log(x/2) I0(x) } = -EUL.
     */
    immutable real[10] chebCoeffsA = [
        1.37446543561352307156e-16L,
        4.25981614279661018399e-14L,
        1.03496952576338420167e-11L,
        1.90451637722020886025e-9L,
        2.53479107902614945675e-7L,
        2.28621210311945178607e-5L,
        1.26461541144692592338e-3L,
        3.59799365153615016266e-2L,
        3.44289899924628486886e-1L,
       -5.35327393233902768720e-1L  ];

    /* Chebyshev coefficients for exp(x) sqrt(x) K0(x)
     * in the inverted interval [2,infinity].
     * 
     * lim(x->inf){ exp(x) sqrt(x) K0(x) } = sqrt(pi/2).
     */
    immutable real[25] chebCoeffsB = [
        5.30043377268626276149e-18L,
       -1.64758043015242134646e-17L,
        5.21039150503902756861e-17L,
       -1.67823109680541210385e-16L,
        5.51205597852431940784e-16L,
       -1.84859337734377901440e-15L,
        6.34007647740507060557e-15L,
       -2.22751332699166985548e-14L,
        8.03289077536357521100e-14L,
       -2.98009692317273043925e-13L,
        1.14034058820847496303e-12L,
       -4.51459788337394416547e-12L,
        1.85594911495471785253e-11L,
       -7.95748924447710747776e-11L,
        3.57739728140030116597e-10L,
       -1.69753450938905987466e-9L,
        8.57403401741422608519e-9L,
       -4.66048989768794782956e-8L,
        2.76681363944501510342e-7L,
       -1.83175552271911948767e-6L,
        1.39498137188764993662e-5L,
       -1.28495495816278026384e-4L,
        1.56988388573005337491e-3L,
       -3.14481013119645005427e-2L,
        2.44030308206595545468e0L   ];

    if (x <= 2)
        return chebEval(x^^2 - 2, chebCoeffsA) - log(x/2) * besselI0(x);
    else
        return exp(-x) * chebEval(8/x - 2, chebCoeffsB) / sqrt(x);
}


unittest
{
    check (approxEqual(besselK0( 1), 4.2102443824070833334e-1, 1.3e-16));
    check (approxEqual(besselK0(10), 1.7780062316167651811e-5, 1.3e-16));
}




/** Modified Bessel function of the second kind, order 1. */
real besselK1(real x) pure
in { assert(x > 0, "Can only calculate K1 for positive argument"); }
body {
    /* Chebyshev coefficients for x(K1(x) - log(x/2) I1(x))
     * in the interval [0,2].
     * 
     * lim(x->0){ x(K1(x) - log(x/2) I1(x)) } = 1.
     */
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

    /* Chebyshev coefficients for exp(x) sqrt(x) K1(x)
     * in the interval [2,infinity].
     *
     * lim(x->inf){ exp(x) sqrt(x) K1(x) } = sqrt(pi/2).
     */
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
        return log(z) * besselI1(x) + chebEval(x^^2 - 2, chebCoeffsA) / x;
    else
        return exp(-x) * chebEval(8.0/x - 2.0, chebCoeffsB) / sqrt(x);
}


unittest
{
    check (approxEqual(besselK1( 1), 6.0190723019723457474e-1, 8.9e-17));
    check (approxEqual(besselK1(10), 1.8648773453825584597e-5, 8.9e-17));
}

