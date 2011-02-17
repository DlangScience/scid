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

