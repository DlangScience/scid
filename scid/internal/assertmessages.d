module scid.internal.assertmessages;

import std.format, std.range;

string fmt( Char, A... )( in Char[] f, in A args ) {
	auto app = appender!string();
	formattedWrite( app, f, args );
	return app.data;
}

mixin template SliceIndexMessages() {
	enum msgPrefix_ = typeof(this).stringof ~ ": ";
	
	string boundsMsg_( size_t i ) const {
		return fmt( msgPrefix_ ~ "Out of bounds %d vs. %d.", i, length );
	}

	string sliceMsg_( size_t s, size_t e ) const {
		return fmt( msgPrefix_ ~ "Invalid slice indices [%d .. %d] vs. %d.", s, e, length );
	}
	
	string sliceAssignMsg_( size_t s, size_t e, size_t rl ) const {
		return fmt( msgPrefix_ ~ "Slice assignment length mismatch [%d .. %d] with length %d vs %d.", s, e, e-s, rl );
	}
	
	string zeroDimMsg_() const {
		return fmt( msgPrefix_ ~ "Zero dimension." );
	}
}

mixin template SliceIndex2dMessages() {
	enum msgPrefix_ = typeof(this).stringof ~ ": ";
	
	string boundsMsg_( size_t i, size_t j ) const {
		return fmt( msgPrefix_ ~ "Out of bounds [%d, %d] vs. [%d, %d].", i, j, rows, columns );
	}

	string sliceMsg_( size_t a, size_t b, size_t c, size_t d ) const {
		return fmt( msgPrefix_ ~ "Invalid slice indices [%d .. %d][%d .. %d] vs. [%d, %d] .", a,b,c,d, rows, columns );
	}
	
	string initMsg_( size_t r, ElementType[] i ) const {
		return fmt( msgPrefix_ ~ "Invalid initializer (%d, %s).", r, i );
	}
	
	string zeroDimMsg_( size_t r, size_t c ) const {
		return fmt( msgPrefix_ ~ "Zero dimension in [%d, %d].", r, c );
	}
}