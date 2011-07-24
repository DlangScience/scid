/** Implementation of the BasicArrayStorage vector storage, which is the default storage type for vectors.

    Authors:    Cristian Cobzarenco
    Copyright:  Copyright (c) 2011, Cristian Cobzarenco. All rights reserved.
    License:    Boost License 1.0
*/
module scid.storage.array;

import scid.vector;

import scid.storage.cowarray;
import scid.storage.arrayview;

import scid.common.traits : isFortranType, BaseElementType;
import std.algorithm      : swap;

import scid.ops.eval;
import scid.ops.common;
import scid.internal.assertmessages;


/** Template that aliases to BasicArrayStorage. If it's given a floating point type it defaults to
    using CowArrayRef as container. Otherwise it uses the container type provided.
*/
template ArrayStorage( ElementOrArray, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias BasicArrayStorage!( CowArrayRef!ElementOrArray, vectorType ) ArrayStorage;
	else
		alias BasicArrayStorage!( ElementOrArray, vectorType )             ArrayStorage;
}

/** This is the storage type used by vectors by default. It is a thin wrapper over a reference to its container
    (the first template parameter, CowArrayRef by default). This provides an extra level of indirection that allows
    ArrayViewStorages to share the reference. This ensures that ArrayViewStorages remain valid even after the
    original ArrayStorage went out of scope.
    The second parameter specifies whether the storage represents a row vector or a column vector.
*/
struct BasicArrayStorage( ContainerRef_, VectorType vectorType_ ) {
	/** The wrapped container reference. */
	alias ContainerRef_ ContainerRef;
	
	/** The type of the elements of the array. */
	alias BaseElementType!ContainerRef ElementType;
	
	/** The vector type i.e. Row or Column. */
	alias vectorType_ vectorType;
	
	/** The type returned by view(). */
	alias BasicArrayViewStorage!( ContainerRef, ArrayViewType.Interval, vectorType ) View;
	
	/** The type returned by view() when passed a stride. */
	alias BasicArrayViewStorage!( ContainerRef, ArrayViewType.Strided, vectorType ) StridedView;
	
	/** The type that represents the transposition of this ArrayStorage - simply negating the vector type. */
	alias BasicArrayStorage!( ContainerRef, transposeVectorType!vectorType ) Transposed;
	
	/** Define a stride and a first index to implement the StridedVectorStorage concept. */
	enum  stride = 1;
	
	/// ditto
	enum  firstIndex = 0;
	
	/** Create a new array of given length. Initialize with zero. */
	this( size_t newLength ) {
		containerRef_ = ContainerRef( newLength );
	}
	
	/** Create a new array of given length. Uninitialized. */
	this( size_t newLength, void* ) {
		containerRef_ = ContainerRef( newLength, null );
	}
	
	/** Create a new array and initialize with a given built-in array. */
	this( ElementType[] initializer ) {
		containerRef_ = ContainerRef( initializer );
	}
	
	/** Postblit c-tor creates a new container initialized with a copy of the source.
	    This ensures copy semantics on assignment.
	*/
	this( this ) {
		if( containerRef_.RefCounted.isInitialized() )
			containerRef_ = ContainerRef( containerRef_.ptr );
	}
	
	/** Forces reference sharing with another Array. This will cause two storages to refer to the same array.
	    This is usually a bad idea - it is used internally for proxy objects. */
	void forceRefAssign( ref typeof(this) rhs ) {
		containerRef_ = rhs.containerRef_;	
	}
		
	/** Implements copy on assignment semnatics. */
	ref typeof(this) opAssign( typeof(this) rhs ) {
		swap( rhs.containerRef_, containerRef_ );
		return this;
	}
	
	/** Element access forwarded to the container. Part of the VectorStorage concept. */
	ElementType index( size_t i ) const
	in {
		// this also silently checks for reference initialization through the call to length().
		assert(	i < length, boundsMsg_( i ) );
	} body {
		return containerRef_.index( i );
	}
	
	/** This method provides both simple opIndexAssign and opIndexOpAssign-like functionality.
	    If the operator template parameter is left empty then it performs a simple assignment.
	    Forwarded to the container.
	*/
	void indexAssign( string op = "" )( ElementType rhs, size_t i )
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		containerRef_.indexAssign!op( rhs, i ) ;
	}
	
	/** Returns a slice of the array. Part of the VectorStorage concept. */
	typeof( this ) slice( size_t start, size_t end ) {
		return typeof(this)( ContainerRef( containerRef_.slice( start, end ).ptr )  );
	}
	
	/** Returns a contiguous view of the array. Part of the VectorStorage concept. */
	View view( size_t start, size_t end ) {
		return typeof( return )( containerRef_, start, end - start );
	}
	
	/** Returns a strided view of the array. Part of the VectorStorage concept. */
	StridedView view( size_t start, size_t end, size_t stride ) {
		size_t len = (end - start);
		len = len / stride + (len % stride != 0);
		return typeof( return )( containerRef_, start, len, stride );
	}
	
	/** Resizes the array to a given size and sets all the elements to zero.
	    If the length is unchanged then no new memory is allocated.
	    Note: The call actually gets forwarded to the container.
	*/
	void resize( size_t rlength ) {
		if( containerRef_.RefCounted.isInitialized() ) 
			containerRef_.resize( rlength );
		else
			containerRef_ = ContainerRef( rlength );
	}
	
	/** Resize the array to a given size and leaves the memory uninitialized.
	    If the length is unchanged then the function does nothing.
		Note: The call actually gets forwarded to the container.
	*/
	void resize( size_t rlength, void* ) {
		if( containerRef_.RefCounted.isInitialized() ) 
			containerRef_.resize( rlength, null );
		else
			containerRef_ = ContainerRef( rlength, null );
	}
	
	/** This storage type uses the common scale(), scaledAddition() and dot() methods for all strided storages.
	    copy() is specialized to forward the call to the container to allow for copy-on-write implementations
	    of the container.
	*/
	void copy( Transpose tr, Source )( auto ref Source rhs ) if( isStridedVectorStorage!(Source,ElementType) ) {
		static if( is( Source : ArrayStorage!( ContainerRef, transposeVectorType!(vectorType,tr) ) ) ) {
			containerRef_.RefCounted.ensureInitialized();
			this.containerRef_ = ContainerRef( rhs.containerRef_.ptr );	
		} else static if( is( Source : ArrayViewStorage!( ContainerRef, transposeVectorType!(vectorType,tr) ) ) ) {
			containerRef_.RefCounted.ensureInitialized();
			this.containerRef_ = rhs.array.slice( rhs.firstIndex, rhs.firstIndex + rhs.length );	
		} else {
			stridedCopy( rhs, this );
		}
	}
	
	@property {
		/** Get a const pointer to the memory used by this storage. */
		const(ElementType)* cdata() const {
			return isInitd_() ? containerRef_.cdata  : null;
		}
		
		/** Get a mutable pointer to the memory used by this storage. */
		ElementType* data() {
			return isInitd_() ? containerRef_.data : null;
		}
		
		/** Forward range methods to the wrapped container, checking that the reference is initialized. */
		bool empty() const {
			return isInitd_() ? containerRef_.empty  : true;
		}
		
		/// ditto
		size_t length() const {
			return isInitd_() ? containerRef_.length : 0;
		}
		
		/// ditto
		ElementType front() const
		in {
			assert( !empty, msgPrefix_ ~ "front() on empty." );
		} body {
			return containerRef_.front;
		}
		
		/// ditto
		ElementType back() const
		in {
			assert( !empty, msgPrefix_ ~ "back() on empty." );
		} body {
			return containerRef_.back;
		}
	}
	
	/// ditto
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront() on empty." );
	} body {
		containerRef_.popFront();
	}
	
	/// ditto
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popBack() on empty." );
	} body {
		containerRef_.popBack();
	}
	
	/** Use the common axpy(), scal() and dot() for all strided storages. */
	mixin StridedScalingAdditionDot;
	
private:
	// The mixin includes some common error messages to use in assertions.
	mixin ArrayErrorMessages;

	// Checks if the wrapped reference is initialized.	
	bool isInitd_() const {
		// TODO: This assumes the reference type is RefCounted. Provide a more general impl.
		return containerRef_.RefCounted.isInitialized();
	}
	
	// Constructor for slice()
	this( ContainerRef containerRef ) {
		swap( containerRef, containerRef_ );	
	}
	
	// The wrapped reference.
	ContainerRef containerRef_;
}

unittest {
	// TODO: Write tests for ArrayStorage.
}

