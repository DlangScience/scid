/** Common error messages to use with assertions.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.internal.assertmessages;

import std.string : format;

/** Error messages for array-like structs. */
mixin template ArrayErrorMessages() {
	/** Prefix for all error messages. */
	private enum msgPrefix_ = typeof(this).stringof ~ ": ";
	
	/** Out of bounds error message. */
	string boundsMsg_( size_t i ) const {
		return format( msgPrefix_ ~ "Out of bounds %d vs. %d.", i, length );
	}

	/** Invalid slice indices error message. */
	string sliceMsg_( size_t s, size_t e ) const {
		return format( msgPrefix_ ~ "Invalid slice indices [%d .. %d] vs. %d.", s, e, length );
	}
	
	/** Invalid slice indices in slice assignment error message. */
	string sliceAssignMsg_( size_t s, size_t e, size_t rl ) const {
		return format( msgPrefix_ ~ "Slice assignment length mismatch [%d .. %d] with length %d vs %d.", s, e, e-s, rl );
	}
	
	/** (pop)front/back called on empty array. */
	string emptyMsg_( string func ) const pure {
		return msgPrefix_ ~ func ~ " called on empty.";
	}
	
	/** Length mismatch in vector operation. */
	string lengthMismatch_( size_t len, string op = "" ) const {
		return format("%sLength mismatch in vector operation '%s': %d vs. %d.", msgPrefix_, op, length, len );
	}
	
	/** Constructed with length = 0. */
	enum zeroDimMsg_ = msgPrefix_ ~ "Zero length in constructor.";
}

/** Error messages for array-like structs (matrix literal, storages and containers). */
mixin template MatrixErrorMessages() {
	/** Prefix for all error messages. */
	enum msgPrefix_ = typeof(this).stringof ~ ": ";
	
	/** Out of bounds error message. */
	string boundsMsg_( size_t i, size_t j ) const {
		return format( msgPrefix_ ~ "Out of bounds [%d, %d] vs. [%d, %d].", i, j, rows, columns );
	}
	
	/** Invalid slice indices error message. */
	string sliceMsg_( size_t a, size_t b, size_t c, size_t d ) const {
		return format( msgPrefix_ ~ "Invalid slice indices [%d .. %d][%d .. %d] vs. [%d, %d] .", a,c,b,d, rows, columns );
	}
	
	/** Invalid initializer error message. */
	string initMsg_( size_t r, ElementType[] i ) const {
		return format( msgPrefix_ ~ "Invalid initializer (%d, %s).", r, i );
	}
	
	/** Dimension mismatch in matrix operation. */
	string dimMismatch_( size_t r, size_t c, string op="" ) const {
		return format("%sDimension mismatch in matrix operation '%s': (%d,%d) vs. (%d,%d).", msgPrefix_, op, rows, columns, r, c );
	}
	
	/** Zero dimension eror message */
	string zeroDimMsg_( size_t r, size_t c ) const {
		return format( msgPrefix_ ~ "Zero dimension in [%d, %d].", r, c );
	}
}
