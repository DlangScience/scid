/** Exceptions and error codes for errors encountered in numerical
    calculations.

    Authors:    Lars Tandle Kyllingstad
    Copyright:  Copyright (c) 2009, Lars T. Kyllingstad. All rights reserved.
    License:    Boost License 1.0
*/
module scid.exception;


import scid.core.testing;




/** Error codes for various numerics errors. */
enum NE
{
    /** An unknown error has occurred. */
    Unknown,

    /** The maximum number of iterations or function evaluations
        was reached.
    */
    Limit,

    /** Failed to reach the requested tolerance/accuracy.
        If possible, try a lower accuracy requirement.
    */
    Accuracy,

    /** The algorithm either converges too slowly, or it diverges. */
    Convergence,

    /** The supplied input to an algorithm was invalid. */
    InvalidInput,

    /** The algorithm fails because of roundoff error. */
    Roundoff,

    /** A function passed to the algorithm behaves badly
        at one or more points.
    */
    Behaviour,

    /** The problem has no solution. */
    NoSolution,
}




/** If the test evaluates to false, throw a NumericsException with
    the specified code and/or message.
*/
void enforceNE(string file = __FILE__, int line = __LINE__)
    (bool test, NE code, lazy string msg)
{
    if (test) return;
    throw new NumericsException(code, msg, file, line);
}


/// ditto
void enforceNE(string file = __FILE__, int line = __LINE__)
    (bool test, NE code = NE.Unknown)
{
    if (test) return;
    throw new NumericsException(code, file, line);
}


/// ditto
void enforceNE(string file = __FILE__, int line = __LINE__)
    (bool test, lazy string msg)
{
    if (test) return;
    throw new NumericsException(msg, file, line);
}


unittest
{
    enforceNE(true);

    try
    {
        enforceNE(false, NE.Limit, "fail");
        check(false);
    }
    catch (NumericsException e)
    {
        check(e.code == NE.Limit);
        check(e.msg == "fail");
    }

    try
    {
        enforceNE(false, NE.Limit);
        check(false);
    }
    catch (NumericsException e)
    {
        check(e.code == NE.Limit);
    }

    try
    {
        enforceNE(false, "fail");
        check(false);
    }
    catch (NumericsException e)
    {
        check(e.code == NE.Unknown);
        check(e.msg == "fail");
    }
}




/** General exception class for numerics errors. */
class NumericsException : Exception
{
    /** A code signifying what kind of error has occurred. */
    immutable NE code;


    /** Constructors. */
    this (string file = null, int line = 0)
    {
        this(NE.Unknown, file, line);
    }

    /// ditto
    this (NE code, string file = null, int line = 0)
    {
        string msg;
        final switch (code)
        {
            case NE.Limit:
                msg = "The maximum number of iterations or function "
                    ~ "evaluations was reached.";
                break;

            case NE.Accuracy:
                msg = "Failed to reach the requested tolerance/accuracy.";
                break;

            case NE.Convergence:
                msg = "Algorithm converges too slowly, or it diverges.";
                break;

            case NE.InvalidInput:
                msg = "Invalid input parameters.";
                break;

            case NE.Roundoff:
                msg = "Roundoff error detected.";
                break;

            case NE.Behaviour:
                msg = "Bad function behaviour detected.";
                break;

            case NE.NoSolution:
                msg = "The problem has no solution.";
                break;

            case NE.Unknown:
                msg = "An unknown error has occurred.";
                break;
        }

        super (msg, file, line);
        this.code = code;
    }

    /// ditto
    this (string msg, string file = null, int line = 0)
    {
        this (NE.Unknown, msg, file, line);
    }

    /// ditto
    this (NE code, string msg, string file = null, int line = 0)
    {
        super (msg, file, line);
        this.code = code;
    }
}
