/** Implementation of the BasicArrayViewStorage vector storage, which is the view type associated with
    BasicArrayStorage, the default storage type for vectors.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.storage.arrayview;

import scid.vector;

import scid.ops.eval, scid.ops.common;
import scid.storage.array, scid.storage.cowarray;
import scid.common.meta;
import scid.common.storagetraits;
import std.traits, std.range, std.algorithm;
import scid.internal.assertmessages;

// import blas functions with the blas. prefix to avoid clashes with
// similarly named methods
static import scid.bindings.blas.dblas;
alias scid.bindings.blas.dblas blas;


/** Enumeration that specifies if the elements in an array view ar contiguous or not. */
enum ArrayViewType {
	Interval,
	Strided
}

/** Template that aliases to a contiguous (ArrayViewType.Interval) BasicArrayViewStorage. If it's given a floating
    point type it defaults to using CowArrayRef as container. If passed a container type is uses that.
*/
template ArrayViewStorage( ElementOrArray, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray ) {
		// if the given type is a floating point no. use CowArrayRef as container
		alias BasicArrayViewStorage!(
			CowArrayRef!ElementOrArray,
			ArrayViewType.Interval,
			vectorType
		) ArrayViewStorage;
	} else {
		// if the given type is a container use that
		alias BasicArrayViewStorage!(
			ElementOrArray,
			ArrayViewType.Interval,
			vectorType
		) ArrayViewStorage;
	}
}

/** Template that aliases to a strided (ArrayViewType.Strided) BasicArrayViewStorage. If it's given a floating
    point type it defaults to using CowArrayRef as container. If passed a container type is uses that.
*/
template StridedArrayViewStorage( ElementOrArray, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray ) {
		alias BasicArrayViewStorage!(
			CowArrayRef!ElementOrArray,
			ArrayViewType.Strided,
			vectorType
		) StridedArrayViewStorage;
	} else {
		alias BasicArrayViewStorage!(
			ElementOrArray,
			ArrayViewType.Strided,
			vectorType
		) StridedArrayViewStorage;
	}
}

/** This is the storage type used by ArrayStorage as view. It is the storage type for VectorView and StridedVectorView.
    Like all storage types it wraps a custom container reference type which is by default a ref-counted CowArray. This
    means that by default views are not invalidated if the original ArrayStorage goes out of scope, but the behaviour
    is of course governed by the container type.
*/
struct BasicArrayViewStorage( ContainerRef_, ArrayViewType strided_, VectorType vectorType_ ) {
	/** The wrapped container reference type. */
	alias ContainerRef_ ContainerRef;
	
	/** The type of the vector elements. */
	alias BaseElementType!ContainerRef ElementType;
	
	/** The orientation of the vector (row/column). */
	alias vectorType_ vectorType;
	
	/** Conceptually an ArrayViewStorage is a reference to a ArrayStorage. This type ensures that when evaluating an
	    expression involving views, the resulting vector uses ArrayStorage as storage, rather than a view. The reason
	    for using ArrayTypeOf is explained in the comment for that template, defined after this struct.
	*/
	alias BasicArrayStorage!( ArrayTypeOf!ContainerRef, vectorType ) Referenced;
	
	/** The result of evaluating a transposition operation on this. */
	alias BasicArrayStorage!( ArrayTypeOf!ContainerRef, transposeVectorType!vectorType ) Transposed;
	
	/** The type returned by view(). */
	alias typeof( this ) View;
	
	/** The type returned by view() when passed a stride. */
	alias BasicArrayViewStorage!( ContainerRef, ArrayViewType.Strided, vectorType ) StridedView;
	
	/** Whether this view is strided (as opposed to contiguous). */
	enum isStrided = (strided_ == ArrayViewType.Strided );
	
	static if( isStrided ) {
		/** Construct a view from a container reference, a start index, length and a stride. */
		this( ref ContainerRef containerRef, size_t iFirstIndex, size_t iLength, size_t iStride )
		in {
			assert( iStride != 0, strideMsg_ );
		} body {
			containerRef_ = containerRef;
			setParams_( iFirstIndex, iLength, iStride );
		}
	} else {
		/** Construct a view from a container reference, a start index and length. */
		this( ref ContainerRef containerRef, size_t iFirstIndex, size_t iLength ) {
			containerRef_ = containerRef;
			setParams_( iFirstIndex, iLength );
		}
	}
	
	/** Forces reference sharing with another Array. This will cause two storages to refer to the same array.
	    This is usually a bad idea - it is used internally for proxy objects. */
	void forceRefSharing( ref typeof(this) rhs ){
		this = rhs;
	}
	
	/** Assignment has reference semantics. */
	ref typeof( this ) opAssign( typeof( this ) rhs ) {
		firstIndex_ = rhs.firstIndex;
		length_     = rhs.length_;
		move( rhs.containerRef_, containerRef_ );
		static if( isStrided )
			stride_ = rhs.stride_;
		return this;
	}
	
	/** Element access forwarded to the container. Part of the VectorStorage concept. */
	ElementType index( size_t i ) const
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		return containerRef_.cdata[ map_( i ) ];
	}
	
	/** This method provides both simple opIndexAssign and opIndexOpAssign-like functionality.
	    If the operator template parameter is left empty then it performs a simple assignment.
	    Forwarded to the container.
	*/
	void indexAssign( string op="" )( ElementType rhs, size_t i )
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		mixin("containerRef_.data[ map_( i ) ] " ~ op ~ "= rhs;");
	}
	
	/** Returns a slice of the array. Part of the VectorStorage concept. */
	typeof(this) slice( size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
	} body {
		static if( isStrided )
			return typeof( this )( containerRef_, map_(start), end - start, stride_ );
		else
			return typeof( this )( containerRef_, map_(start), end - start );
	}
	
	/** Returns a another contiguous view of the array. Part of the VectorStorage concept. */
	View view( size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_(start, end) );
	} body {
		static if( isStrided )
			return typeof( return )( containerRef_, start, end-start, stride_ );
		else
			return typeof( return )( containerRef_, start, end-start );
	}
	
	/** Returns another strided view of the array. Part of the VectorStorage concept. */
	StridedView view( size_t start, size_t end, size_t newStride )
	in {
		assert( newStride != 0, strideMsg_ );
		assert( start < end && end <= length, sliceMsg_(start, end) );
	} body {
		size_t len = (end - start);
		len = len / newStride + ( (len % newStride) != 0 );
		static if( isStrided )
			newStride *= stride_;
		return typeof( return )( containerRef_, start, len, newStride );
	}
	
	/** Check if the length is the same as the one given and zero out all the elements. Part of the VectorStorage concept. */
	void resize( size_t rlength ) {
		resize( rlength, null );
		blas.scal( length_, Zero!ElementType, data, stride );
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
	
	/** Forward range methods to the wrapped container. */
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront() on empty." );
	} body {
		
		static if( isStrided ) {
			firstIndex_ += stride_;
		} else
			++ firstIndex_;
		-- length_;
	}
	
	/// ditto
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ ": popBack() on empty." );
	} body {
		-- length_;
	}
	
	@property {
		/** Get a const pointer to the memory used by this storage. */
		ElementType* data() {
			return containerRef_.data + firstIndex_;
		}
		
		/** Get a mutable pointer to the memory used by this storage. */
		const(ElementType)* cdata() const {
			return containerRef_.cdata + firstIndex_;
		}
		
		/** The index in the array at which this view starts. */
		size_t firstIndex() const {
			return firstIndex_;
		}
		
		/** Forward range methods to the wrapped container, checking that the reference is initialized. */
		size_t length() const {
			return length_;
		}
		
		/// ditto
		bool empty() const {
			return length_ == 0;
		}
		
		/// ditto
		void front( ElementType newValue )
		in {
			assert( !empty, msgPrefix_ ~ "front assign on empty." );
		} body {
			indexAssign( newValue, 0 );
		}
		
		/// ditto
		void back( ElementType newValue  )
		in {
			assert( !empty, msgPrefix_ ~ "back assign on empty." );
		} body {
			indexAssign( newValue, length - 1 );
		}
		
		/// ditto
		ElementType front() const
		in {
			assert( !empty, msgPrefix_ ~ "front get on empty." );
		} body {
			return this.index( 0 );
		}
		
		/// ditto
		ElementType back() const
		in {
			assert( !empty, msgPrefix_ ~ "back get on empty." );
		} body {
			return index( length_ - 1 );
		}
		
		/** The stride of the view i.e. the index difference between two consecutive elements. */
		static if( isStrided ) {
			size_t stride() const { return stride_; }
		} else {
			enum stride = 1;
		}
	}
	
	/** Promotions for this type are inherited from ArrayStorage */
	template Promote( Other ) {
		alias Promotion!( BasicArrayStorage!(ArrayTypeOf!ContainerRef, vectorType), Other ) Promote;
	}
	
private:
	mixin ArrayErrorMessages;

	enum strideMsg_ = msgPrefix_ ~ "Zero stride is impossible.";

	static if( isStrided ) {		
		void setParams_( size_t f, size_t l, size_t s ) {
			firstIndex_ = f;
			length_     = l;
			stride_ = s;
		}
		
		size_t map_( size_t i ) const {
			return i * stride_ + firstIndex_;	
		}
		
		size_t stride_;
	} else {
		void setParams_( size_t f, size_t l ) {
			firstIndex_ = f;
			length_     = l;
		}
		
		size_t map_( size_t i ) const {
			return i + firstIndex_;
		}
	}
	
	size_t firstIndex_, length_;
	ContainerRef containerRef_;
}

unittest {
	// TODO: Write tests for ArrayViewStorage.
}