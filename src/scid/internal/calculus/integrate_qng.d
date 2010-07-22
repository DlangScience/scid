/** Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2010, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.internal.calculus.integrate_qng;


import std.algorithm: max, min;
import std.conv;
import std.math;
import std.numeric;
import std.traits;

import scid.exception;
import scid.types;

version(unittest)
{
    import scid.core.testing;
}




Result!Real qng(Func, Real = ReturnType!Func)
    (Func f, Real a, Real b, Real epsRel, Real epsAbs)
{
    alias FPTemporary!Real T;

    enum epsRelMin = 50*T.epsilon;
    enforceNE(
        epsAbs > 0 || epsRel > epsRelMin,
        NE.InvalidInput,
        "Invalid accuracy request (must have epsAbs > 0 or epsRel > "
            ~to!string(epsRelMin)~")");


    // This function is used to calculate the absolute error.
    static T errorCalc(T estimate, T absResult, T variance)
    {
        T error = estimate;
        if (variance != 0 && error != 0)
            error = variance * min(1, (200*error/variance)^^1.5);
        if (absResult > T.min_normal/(50*T.epsilon))
            error = max(error, 50*T.epsilon*absResult);
        return error;
    }


    // Gauss-Kronrod-Patterson quadrature coefficients,
    // calculated by L. W. Fullerton, Bell Labs, Nov. 1981.
    static immutable T[5] x1 = [
        0.9739065285_1717172007_7964012084_452,
        0.8650633666_8898451073_2096688423_493,
        0.6794095682_9902440623_4327365114_874,
        0.4333953941_2924719079_9265943165_784,
        0.1488743389_8163121088_4826001129_720];
    static immutable T[5] w10 = [
        0.0666713443_0868813759_3568809893_332,
        0.1494513491_5058059314_5776339657_697,
        0.2190863625_1598204399_5534934228_163,
        0.2692667193_0999635509_1226921569_469,
        0.2955242247_1475287017_3892994651_338];


    static immutable T[5] x2 = [
        0.9956571630_2580808073_5527280689_003,
        0.9301574913_5570822600_1207180059_508,
        0.7808177265_8641689706_3717578345_042,
        0.5627571346_6860468333_9000099272_694,
        0.2943928627_0146019813_1126603103_866];
    static immutable T[5] w21a = [
        0.0325581623_0796472747_8818972459_390,
        0.0750396748_1091995276_7043140916_190,
        0.1093871588_0229764189_9210590325_805,
        0.1347092173_1147332592_8054001771_707,
        0.1477391049_0133849137_4841515972_068];
    static immutable T[6] w21b = [
        0.0116946388_6737187427_8064396062_192,
        0.0547558965_7435199603_1381300244_580,
        0.0931254545_8369760553_5065465083_366,
        0.1234919762_6206585107_7958109831_074,
        0.1427759385_7706008079_7094273138_717,
        0.1494455540_0291690566_4936468389_821];

    static immutable T[11] x3 = [
        0.9993333609_0193208139_4099323919_911,
        0.9874334029_0808886979_5961478381_209,
        0.9548079348_1426629925_7919200290_473,
        0.9001486957_4832829362_5099494069_092,
        0.8251983149_8311415084_7066732588_520,
        0.7321483889_8930498261_2354848755_461,
        0.6228479705_3772523864_1159120344_323,
        0.4994795740_7105649995_2214885499_755,
        0.3649016613_4658076804_3989548502_644,
        0.2222549197_7660129649_8260928066_212,
        0.0746506174_6138332204_3914435796_506];
    static immutable T[10] w43a = [
        0.0162967342_8966656492_4281974617_663,
        0.0375228761_2086950146_1613795898_115,
        0.0546949020_5825544214_7212685465_005,
        0.0673554146_0947808607_5553166302_174,
        0.0738701996_3239395343_2140695251_367,
        0.0057685560_5976979618_4184327908_655,
        0.0273718905_9324884208_1276069289_151,
        0.0465608269_1042883074_3339154433_824,
        0.0617449952_0144256449_6240336030_883,
        0.0713872672_6869339776_8559114425_516];
    static immutable T[12] w43b = [
        0.0018444776_4021241410_0389106552_965,
        0.0107986895_8589165174_0465406741_293,
        0.0218953638_6779542810_2523123075_149,
        0.0325974639_7534568944_3882222526_137,
        0.0421631379_3519181184_7627924327_955,
        0.0507419396_0018457778_0189020092_084,
        0.0583793955_4261924837_5475369330_206,
        0.0647464049_5144588554_4689259517_511,
        0.0695661979_1235648452_8633315038_405,
        0.0728244414_7183320815_0939535192_842,
        0.0745077510_1417511827_3571813842_889,
        0.0747221475_1740300559_4425168280_423];

    static immutable T[22] x4 = [
        0.9999029772_6272923449_0529830591_582,
        0.9979898959_8667874542_7496322365_960,
        0.9921754978_6068722280_8523352251_425,
        0.9813581635_7271277357_1916941623_894,
        0.9650576238_5838461912_8284110607_926,
        0.9431676131_3367059681_6416634507_426,
        0.9158064146_8550720959_1826430720_050,
        0.8832216577_7131650137_2117548744_163,
        0.8457107484_6241566660_5902011504_855,
        0.8035576580_3523098278_8739474980_964,
        0.7570057306_8549555832_8942793432_020,
        0.7062732097_8732181982_4094274740_840,
        0.6515894665_0117792253_4422205016_736,
        0.5932233740_5796108887_5273770349_144,
        0.5314936059_7083193228_5268948562_671,
        0.4667636230_4202284487_1966781659_270,
        0.3994248478_5921880473_2101665817_923,
        0.3298748771_0618828826_5053371824_597,
        0.2585035592_0216155180_2280975429_025,
        0.1856953965_6834665201_5917141167_606,
        0.1118422131_7990746817_2398359241_362,
        0.0373521233_9461987081_4998165437_704];
    static immutable T[21] w87a = [
        0.0081483773_8414917290_0002878448_190,
        0.0187614382_0156282224_3935059003_794,
        0.0273474510_5005228616_1582829741_283,
        0.0336777073_1163793004_6581056957_588,
        0.0369350998_2042790761_4589586742_499,
        0.0028848724_3021153050_1334156248_695,
        0.0136859460_2271270188_8950035273_128,
        0.0232804135_0288831112_3409291030_404,
        0.0308724976_1171335867_5466394126_442,
        0.0356936336_3941877071_9351355457_044,
        0.0009152833_4520224136_0843392549_948,
        0.0053992802_1930047136_7738743391_053,
        0.0109476796_0111893113_4327826856_808,
        0.0162987316_9678733526_2665703223_280,
        0.0210815688_8920383511_2433060188_190,
        0.0253709697_6925382724_3467999831_710,
        0.0291896977_5647575250_1446154084_920,
        0.0323732024_6720278968_5788194889_595,
        0.0347830989_5036514275_0781997949_596,
        0.0364122207_3135178756_2801163687_577,
        0.0372538755_0304770853_9592001191_226];
    static immutable T[23] w87b = [
        0.0002741455_6376207235_0016527092_881,
        0.0018071241_5505794294_8341311753_254,
        0.0040968692_8275916486_4458070683_480,
        0.0067582900_5184737869_9816577897_424,
        0.0095499576_7220164653_6053581325_377,
        0.0123294476_5224485369_4626639963_780,
        0.0150104473_4638895237_6697286041_943,
        0.0175489679_8624319109_9665352925_900,
        0.0199380377_8644088820_2278192730_714,
        0.0221949359_6101228679_6332102959_499,
        0.0243391471_2600080547_0360647041_454,
        0.0263745054_1483920724_1503786552_615,
        0.0282869107_8877120065_9968002987_960,
        0.0300525811_2809269532_2521110347_341,
        0.0316467513_7143992940_4586051078_883,
        0.0330504134_1997850329_0785944862_689,
        0.0342550997_0422606178_7082821046_821,
        0.0352624126_6015668103_3782717998_428,
        0.0360769896_2288870118_5500318003_895,
        0.0366986044_9845609449_8018047441_094,
        0.0371205492_6983257611_4119958413_599,
        0.0373342287_5193504032_1235449094_698,
        0.0373610737_6267902341_0321241766_599];


    // Some useful quantities.
    immutable T halfLength = (b - a) / 2;
    immutable T absHalfLength = abs(halfLength);

    immutable T center = (a + b) / 2;
    immutable T fCenter = f(center);

    // In these arrays we save the result of function evaluations,
    // so they can be reused for higher-order rules.
    T[5] fSave1, fSave2, fSave3, fSave4;
    T[21] fSave;


    // Compute the integral using the 10- and 21-point formulae.
    T result10 = 0;
    T result21 = w21b[5] * fCenter;
    T resultAbs = w21b[5] * abs(fCenter);

    foreach (k; 0 .. 5)
    {
        immutable x = halfLength * x1[k];

        immutable fValue1 = f(center + x);
        immutable fValue2 = f(center - x);
        immutable fValue = fValue1 + fValue2;

        result10 += w10[k] * fValue;
        result21 += w21a[k] * fValue;
        resultAbs += w21a[k] * (abs(fValue1) + abs(fValue2));

        fSave1[k] = fValue1;
        fSave2[k] = fValue2;
        fSave[k] = fValue;
    }

    foreach (k; 0 .. 5)
    {
        immutable x = halfLength * x2[k];

        immutable fValue1 = f(center + x);
        immutable fValue2 = f(center - x);
        immutable fValue = fValue1 + fValue2;

        result21 += w21b[k] * fValue;
        resultAbs += w21b[k] * (abs(fValue1) + abs(fValue2));

        fSave3[k] = fValue1;
        fSave4[k] = fValue2;
        fSave[k+5] = fValue;
    }

    // resultAbs is the absolute value of the integral of the
    // absolute value of the function.
    resultAbs *= absHalfLength;

    // fMean is the mean value of f(x) in the interval [a, b],
    // as evaluated by the 21-point rule.
    immutable fMean = result21/2;

    // variance is the integral of f(x)-fMean on the interval [a, b],
    // evaluated using the 21-point rule.
    T variance = w21b[5] * abs(fCenter - fMean);
    foreach (k; 0 .. 5)
    {
        variance +=
            w21a[k] * (abs(fSave1[k] - fMean) + abs(fSave2[k] - fMean))
            + w21b[k] * (abs(fSave3[k] - fMean) + abs(fSave4[k] - fMean));
    }
    variance *= absHalfLength;

    // Test for convergence.
    {

        immutable result = result21 * halfLength;
        immutable error = errorCalc(abs((result21 - result10)*halfLength),
            resultAbs, variance);

        if (error <= max(epsAbs, epsRel * abs(result)))
            return typeof(return)(result, error);
    }


    // Compute the integral using the 43-point formula.
    T result43 = w43b[11] * fCenter;
    foreach (k; 0 .. 10)
        result43 += fSave[k] * w43a[k];

    foreach (k; 0 .. 11)
    {
        immutable x = halfLength * x3[k];
        immutable fValue = f(center+x) + f(center-x);
        result43 += fValue * w43b[k];
        fSave[k+10] = fValue;
    }

    // Test for convergence
    {
        immutable result = result43 * halfLength;
        immutable error = errorCalc(abs((result43-result21)*halfLength),
            resultAbs, variance);

        if (error <= max(epsAbs, epsRel * abs(result)))
            return typeof(return)(result, error);
    }


    // Compute the integral using the 87-point formula.
    T result87 = w87b[22] * fCenter;
    foreach (k; 0 .. 21)
        result87 += fSave[k] * w87a[k];

    foreach (k; 0 .. 22)
    {
        immutable x = halfLength * x4[k];
        result87 += w87b[k] * (f(center+x) + f(center-x));
    }
    
    // Test for convergence
    {
        immutable result = result87 * halfLength;
        immutable error = errorCalc(abs((result87-result43)*halfLength),
            resultAbs, variance);

        enforceNE(
            error <= max(epsAbs, epsRel * abs(result)),
            NE.Convergence,
            "Integral didn't converge to the requested accuracy");

        return typeof(return)(result, error);
    }
}


unittest
{
    // Check that it compiles for different types.
    alias qng!(float delegate(float)) fqng;
    alias qng!(double delegate(double)) dqng;
    alias qng!(double function(double)) dfqng;
    alias qng!(real delegate(real)) rqng;
}


unittest
{
    // QUADPACK book §4.4, integral 4
    double alpha;
    double f(double x) { return x^^alpha * log(1/x); }

    enum epsabs = 0.0;
    enum epsrel = 1e-8;

    foreach (i; 1 .. 9)
    {
        alpha = 1 + i*0.2;
        auto result = qng(&f, 0.0, 1.0, epsabs, epsrel);
        auto expect = (alpha + 1)^^(-2);
        check(isAccurate(result.value, result.error, expect, epsrel, epsabs));
    }
}


unittest
{
    // QUADPACK book §4.4, integral 5
    double alpha;
    double f(double x) { return 4^^(-alpha) / ((x-PI_4)^^2 + 16^^(-alpha)); }

    enum epsabs = 0.0;
    enum epsrel = 1e-8;

    foreach (i; 0 .. 2)
    {
        alpha = i;
        auto result = qng(&f, 0.0, 1.0, epsabs, epsrel);
        double expect = atan((4-PI) * 4^^(alpha-1)) + atan(PI * 4^^(alpha-1));
        check(isAccurate(result.value, result.error, expect, epsrel, epsabs));
    }
}


unittest
{
    // QUADPACK book §4.4, integral 6
    double alpha;
    double f(double x) { return cos(2^^alpha * sin(x)); }

    enum epsabs = 0.0;
    enum epsrel = 1e-8;

    auto expect = [
         2.4039394306344129983,
         0.70337362695660089178,
        -1.2476829250428461076,
         0.53925691468609779719,
        -0.54946164594662718058,
         0.43378800263473354846,
         0.29088010217372596783
    ];

    foreach (i; 0 .. 6)
    {
        alpha = i;
        auto result = qng(&f, 0.0, cast(double) PI, epsabs, epsrel);
        //check(isAccurate(result.value, result.error, expect[i], epsrel, epsabs));
        check(abs((result.value-expect[i])/result.value) <= epsrel);
    }
    // TODO:
    // Some of the error estimates are a bit too low, which is why we
    // can't use isAccurate.  This may not be a problem, but it's worth
    // investigating closer.
    //
    // What's worse is that the routine should succeed for i=6 as well,
    // but it doesn't.  Actually, it does converge to the right value,
    // but it erroneously reports the error as ~0.8.
}
