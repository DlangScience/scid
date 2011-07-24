/** Common error messages to use with assertions.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.internal.assertmessages;

import std.format : formattedWrite;
import std.range  : appender;

/** Error messages for array-like structs. */
mixin template ArrayErrorMessages() {
	/** Prefix for all error messages. */
	private enum msgPrefix_ = typeof(this).stringof ~ ": ";
	
	/** Out of bounds error message. */
	string boundsMsg_( size_t i ) const {
		return formatString( msgPrefix_ ~ "Out of bounds %d vs. %d.", i, length );
	}

	/** Invalid slice indices error message. */
	string sliceMsg_( size_t s, size_t e ) const {
		return formatString( msgPrefix_ ~ "Invalid slice indices [%d .. %d] vs. %d.", s, e, length );
	}
	
	/** Invalid slice indices in slice assignment error message. */
	string sliceAssignMsg_( size_t s, size_t e, size_t rl ) const {
		return formatString( msgPrefix_ ~ "Slice assignment length mismatch [%d .. %d] with length %d vs %d.", s, e, e-s, rl );
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
		return formatString( msgPrefix_ ~ "Out of bounds [%d, %d] vs. [%d, %d].", i, j, rows, columns );
	}
	
	/** Invalid slice indices error message. */
	string sliceMsg_( size_t a, size_t b, size_t c, size_t d ) const {
		return formatString( msgPrefix_ ~ "Invalid slice indices [%d .. %d][%d .. %d] vs. [%d, %d] .", a,b,c,d, rows, columns );
	}
	
	/** Invalid initializer error message. */
	string initMsg_( size_t r, ElementType[] i ) const {
		return formatString( msgPrefix_ ~ "Invalid initializer (%d, %s).", r, i );
	}
	
	/** Zero dimension eror message */
	string zeroDimMsg_( size_t r, size_t c ) const {
		return formatString( msgPrefix_ ~ "Zero dimension in [%d, %d].", r, c );
	}
}

// short hand for formattedWrite with a new appender.
string formatString( Char, A... )( in Char[] f, in A args ) {
	auto app = appender!string();
	formattedWrite( app, f, args );
	return app.data;
}
