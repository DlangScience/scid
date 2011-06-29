module scid.vectorview;

import scid.internal.cowarray;
import std.stdio;
import std.conv;
import std.range;

import std.typecons : RefCounted;

/** Implements a strided view of a Vector. Can be used directly with BLAS, with increment == stride_. 
  * Has reference semantics.
  */
struct VectorView( Wrapped ) {
	alias RefCounted!Wrapped  Array;
	alias Wrapped.ElementType ElementType;
	
	/** Create a view of a given ArrayStorage. */
	this( Array array, size_t a, size_t b, size_t stride = 1 ) {
		reset( array, a, b, stride );
	}
	
	/** Create a view of another view. */
	this( ref typeof(this) view, size_t a, size_t b, size_t stride = 1 ) {
		reset( view, a, b, stride );
	}
	
	/** Explictly construct a view from an existing one. */
	this( typeof(this) view ) {
		this = view;
	}
	
	/** Reset to a view of a given Array. */
	void reset( Array array, size_t a, size_t b, size_t stride = 1 ) {
		array_      = array;
		firstIndex_ = a;
		length_     = (b - a) / stride;
		if( ( b - a ) % stride )
			length_ ++;
		stride_     = stride;
	}
	
	/** Reset to a view of another view. */
	void reset( ref typeof(this) view, size_t a, size_t b, size_t stride = 1 ) {
		array_      = view.array_;
		firstIndex_ = a * view.stride_ + view.firstIndex_;
		stride_     = stride * view.stride_;
		length_     = (b - a) / stride;
		if( ( b - a ) % stride )
			length_ ++;
	}
	
	/** Assign to another view. */
	ref typeof(this) opAssign( typeof(this) view ) {
		array_      = view.array_;
		firstIndex_ = view.firstIndex_;
		length_     = view.length_;
		stride_     = view.stride_;   
		return this;
	}
	
	/** Element access. */
	ElementType opIndex( size_t i ) const
	in {
		assert( i < length_, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string(length) );
	} body {
		return array_.refCountedPayload()[ map_( i ) ];
	}
	
	/// ditto
	void opIndexAssign( ElementType rhs, size_t i )
	in {
		assert( i < length_, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string(length) );
	} body {
		array_.refCountedPayload()[ map_( i ) ] = rhs;
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t i )
	in {
		assert( i < length_, "Out of bounds " ~ to!string(i) ~ " >= " ~ to!string(length) );
	} body {
		mixin( "array[ map_( i ) ] " ~ op ~"= rhs;" );
	}
	
	/** Returns a slice of the view. */
	typeof(this) view( size_t a, size_t b, size_t stride = 1 ) {
		return typeof(this)( this, a, b, stride );
	}
	
	/// ditto
	typeof(this) opSlice() {
		return view( 0, length, 1 );
	}

	/// ditto
	typeof(this) opSlice( size_t a, size_t b ) {
		return view( a, b, 1 );
	}
	
	/** Copies a range to a slice of the view. */
	void opSliceAssign( S )( S rhs, size_t a, size_t b ) if( isInputRange!S ) {
		
		static if( __traits( compiles, rhs.length ) )
			assert( rhs.length == b - a );
		
		size_t end = map_( b );
		for( size_t i = map_( a ); i < end; i += stride_ ) {
			array_.refCountedPayload()[ i ] = rhs.front();
			rhs.popFront();
		}
	}
	
	/** Copies a range to a slice of the view. */
	void opSliceAssign( S )( S rhs ) if( isInputRange!S ) {
		static if( __traits( compiles, rhs.length ) )
			assert( rhs.length == length_ );
		
		int end = map_(length_);
		for( int i = firstIndex_; i < end; i += stride_ ) {
			array_.refCountedPayload()[ i ] = rhs.front();
			rhs.popFront();
		}
		
		
	}
	
	/** Fills a slice of the view with a value. */
	void opSliceAssign( S )( S rhs ) if( is( ElementType == S ) ) {
		int end = map_(length_);
		for( int i = firstIndex_; i < end; i += stride_ )
			array_.refCountedPayload()[ i ] = rhs;
	}
	
	/** Fills a slice of the view with a value. */
	void opSliceAssign( S )( S rhs, size_t a, size_t b ) if( is( ElementType == S ) ) {
		int end = map_( b );
		for( int i = firstIndex_; i < end; i += stride_ )
			array_.refCountedPayload()[ i ] = rhs;
	}
	
	/** Return the wrapped memory block. Writable. */
	@property ElementType[] data() {
		assert(array_.RefCounted.isInitialized());
		return array_.data[ map_(0) .. map_(length_) ];
	}
	
	/** Return the wrapped memory block. Read-only. */
	@property const(ElementType) [] cdata() const {
		assert(array_.RefCounted.isInitialized());
		return array_.cdata[ map_(0) .. map_(length_) ];
	}
	
	/** Return the underlying Array. */
	@property Array       array()          { return array_; }
	
	/** Return the stride of the view. */
	@property size_t      stride()         { return stride_; }
	
	/** Return the first index of the view. */
	@property size_t      firstIndex()     { return firstIndex_; }
	
	/** Return the length of the view. */
	@property size_t      length()   const { return length_; }
	
	/** BidirectionalRange method. Returns the first element in the array. */
	@property ElementType front()    const { return this[ 0 ]; }
	
	/** BidirectionalRange method. Returns whether the array is empty. */
	@property bool        empty()    const { return length_ == 0; }
	
	/** BidirectionalRange method. Returns the last element in the array. */
	@property ElementType back()     const { return this[ length_ - 1 ]; }
	
	/** BidirectionalRange method. Removes the first element. */
	void popFront()
	in {
		assert( !empty );
	} body {
		firstIndex_ += stride_;
		length_--;
	}
	
	/** BidirectionalRange method. Removes the last element. */
	void popBack()
	in {
		assert( !empty );
	} body {
		length_--;
	}
	
	// Map index to the wrapped array.
	private size_t map_( size_t i ) const {
		return i * stride_ + firstIndex_;
	}
	
	private size_t firstIndex_, length_, stride_;
	private Array  array_;
}