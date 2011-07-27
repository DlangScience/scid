module scid.storage.extarrayview;

import scid.vector;
import scid.internal.assertmessages;
import scid.ops.common;
import scid.ops.expression;
import scid.common.meta;
import scid.common.traits;

// import blas functions with the blas. prefix to avoid clashes with
// similarly named methods
static import scid.bindings.blas.dblas;
alias scid.bindings.blas.dblas blas;

import std.string, std.array;

struct ExternalArrayView( T, VectorType vectorType_ = VectorType.Column ) {
	alias typeof( this ) ContainerRef;
	alias T ElementType;
	alias vectorType_ vectorType;
	alias ExternalArrayView!( T, transposeVectorType!vectorType ) Transposed;
	alias ExternalArrayView!( T, vectorType  ) StridedView;
	alias typeof( this ) View;
	
	this( Allocator )( size_t len, Allocator *allocator )
			if( isAllocator!Allocator ) {
		array_ = ( cast(T*) allocator.allocate( len * T.sizeof ) )[ 0 .. len ];
		stride_ = 1;
	}
	
	this()( T[] data, size_t stride )
	in {
		assert( !data.empty, format("%sNull/empty data.", msgPrefix_) );
		assert( data.length % stride == 0,
			format("%sLength of data (%d) is not divisible by the stride (%d).", msgPrefix_, data.length, stride ) );
	} body {
		array_  = data;
		stride_ = stride;
	}
	
	this()( T* data, size_t firstIndex, size_t length, size_t stride = 1 )
	in {
		assert( stride != 0, msgPrefix_~"Zero stride." );
		assert( data != null, msgPrefix_ ~ "Null data pointer." );
		assert( length != 0, msgPrefix_ ~ "Zero length." );
	} body {
		array_  = data[ firstIndex .. length * stride ];
		stride_ = stride;
	}
	
	void forceRefSharing( ref typeof(this) rhs ) {
		array_  = rhs.array_;
		stride_ = rhs.stride_;
	}
	
	ref typeof( this ) opAssign( typeof(this) rhs ) {
		array_  = rhs.array_;
		stride_ = rhs.stride_;
		return this;
	}
	
	/** Element access forwarded to the container. Part of the VectorStorage concept. */
	ElementType index( size_t i ) const
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		return array_[ i ];
	}
	
	/** This method provides both simple opIndexAssign and opIndexOpAssign-like functionality.
	    If the operator template parameter is left empty then it performs a simple assignment.
	    Forwarded to the container.
	*/
	void indexAssign( string op="" )( ElementType rhs, size_t i )
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		mixin("array_[ i * stride_  ] " ~ op ~ "= rhs;");
	}
	
	typeof(this) slice( size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
	} body {
		assert(false);
	}
	
	/** Returns a another contiguous view of the array. Part of the VectorStorage concept. */
	View view( size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_(start, end) );
	} body {
		assert(false);
	}
	
	/** Returns another strided view of the array. Part of the VectorStorage concept. */
	StridedView view( size_t start, size_t end, size_t newStride )
	in {
		assert( newStride != 0, msgPrefix_ ~ "Zero stride." );
		assert( start < end && end <= length, sliceMsg_(start, end) );
	} body {
		assert(false);
	}
	
	/** Check if the length is the same as the one given and zero out all the elements. Part of the VectorStorage concept. */
	void resize( size_t rlength ) {
		resize( rlength, null );
		blas.scal( length, Zero!ElementType, data, stride );
	}
	
	/** Check if the length is the same as the one given. Part of the VectorStorage concept. */
	void resize( size_t rlength, void* ) {
		assert( length == rlength, "Length mismatch in vector operation." );
	}
	
	/** Copy specialization. */
	void copy( Transpose tr, S )( auto ref S rhs ) if( isStridedVectorStorage!(S, ElementType) ) {
		stridedCopy( rhs, this );
	}
	
	/** Use the common scale(), scaledAddition() and dot() methods for strided storages. */
	mixin StridedScalingAdditionDot;
	
	@property {
		size_t        stride() const { return stride_; }
		size_t        length() const { return array_.length; }
		bool          empty()  const { return array_.empty; }
		T             front()  const { return array_.front; }
		T             back()   const { return array_.back; }
		const(T)*     cdata()  const { return array_.ptr; }
		T*            data()         { return array_.ptr; }
		typeof(this)* ptr()          { return &this; }
		enum          firstIndex     = 0;
	}
	
	void popFront() {
		array_.popFront();
	}
	
	void popBack() {
		array_.popBack();
	}
	
private:
	mixin ArrayErrorMessages;

	T[]    array_;
	size_t stride_;
}