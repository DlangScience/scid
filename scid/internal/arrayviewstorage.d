module scid.internal.arrayviewstorage;

import scid.internal.arraystorage;
import std.conv;


struct ArrayViewStorage( T ) {
	alias typeof( this ) View;
	alias T              ElementType;
	alias ArrayStorage!T Array;
	
	/** Create a view of a given ArrayStorage. */
	this( Array array, size_t a, size_t b, size_t stride = 1 ) {
		assign( array, a, b, stride );
	}
	
	/** Create a view of another view. */
	this( ref View view, size_t a, size_t b, size_t stride = 1 ) {
		assign( view, a, b, stride );
	}
	
	this( View view ) {
		assign( view );
	}
	
	/** Returns true if the view isn't pointing at anything. */
	bool isNull() const { return stride_ == 0; }

	/** Assign to a view of a given ArrayStorage. */
	void assign( Array array, size_t a, size_t b, size_t stride = 1 ) {
		array_      = array;
		firstIndex_ = a;
		length_     = (b - a) / stride;
		stride_     = stride;
	}
	
	/** Assign to a view of another view. */
	void assign( ref View view, size_t a, size_t b, size_t stride = 1 ) {
		array_      = view.array_;
		firstIndex_ = a * view.stride_ + view.firstIndex_;
		stride_     = stride * view.stride_;
		length_     = (b - a) / stride;
	}
	
	/** Storage concept. Assign to another view. */
	void assign( ref View view ) {
		array_      = view.array_;
		firstIndex_ = view.firstIndex_;
		length_     = view.length_;
		stride_     = view.stride_;   
	}
	
	/** Storage concept. Assign to another array. */
	void assign( Array array ) {
		array_      = array;
		firstIndex_ = 0;
		length_     = array.length;
		stride_     = 1;
	}
	
	/** Storage concept. Returns the ith element. */
	T at( size_t i ) const
	in {
		assert( i < length_ );
		assert( !isNull() );
	} body {
		return array_.at( i * stride_ + firstIndex_ );
	}
	
	/** Storage concept. Sets the ith element. */
	void set( T rhs, size_t i )
	in {
		assert( i < length_ );
		assert( !isNull() );
	} body {
		array_.set( rhs, i * stride_ + firstIndex_ );
	}
	
	/** Storage concept. Returns a slice of the view. */
	View slice( size_t a, size_t b, size_t stride = 1 ) {
		return View( this, a, b, stride );
	}
	
	/** Storage concept. Copies a range to a  slice of the view. */
	void sliceAssign( S )( S rhs, size_t a, size_t b ) if( isInputRange!S ) {
		static if( __traits( compiles, rhs.length ) )
			assert( rhs.length == b - a );
		
		int end = firstIndex_ + b * stride_;
		for( int i = firstIndex_ + a * stride_; i < end; i += stride_ ) {
			array_[ i ] = rhs.front;
			rhs.popFront();
		}
	}
	
	/** Storage concept. Fills a slice of the view with a value. */
	void sliceAssign( S )( S rhs, size_t a, size_t b ) if( is( T == S ) ) {
		int end = firstIndex_ + b * stride_;
		for( int i = firstIndex_ + a * stride_; i < end; i += stride_ )
			array_[ i ] = rhs;
	}
	
	/** Return the underlying ArrayStorage. */
	@property Array   array()          { return array_; }
	
	/** Return the stride of the view. */
	@property size_t  stride()         { return stride_; }
	
	/** Return the first index of the view. */
	@property size_t  firstIndex()     { return firstIndex_; }
	
	/** Return the length of the view. */
	@property size_t  length()   const { return length_; }
	
	/** InputRange method. Returns the first element in the array. */
	@property T       front()    const { return at( 0 ); }
	
	/** InputRange method. Returns whether the array is empty. */
	@property bool    empty()    const { return length_ == 0; }
	
	/** InputRange method. Returns the last element in the array. */
	@property T       back()     const { return at( length_ - 1 ); }
	
	/** InputRange method. Removes the first element. */
	void popFront()
	in {
		assert( !empty );
	} body {
		firstIndex_ += stride_;
		length_--;
	}
	
	/** InputRange method. Removes the last element. */
	void popBack()
	in {
		assert( !empty );
	} body {
		length_--;
	}
	
	string toString() {
		auto v = this;
		if( v.empty )
			return "[]";
		
		string r = '[' ~ to!string( v.front );
		v.popFront();
		foreach( x; v )
			r ~= ", " ~ to!string( x );
		r ~= ']';
		return r;
	}
	
	private size_t         firstIndex_, length_, stride_;
	private ArrayStorage!T array_;
}