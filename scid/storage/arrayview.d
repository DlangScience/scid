module scid.storage.arrayview;

static import scid.bindings.blas.dblas;
alias scid.bindings.blas.dblas blas;

import scid.internal.assertmessages;
import scid.internal.hlblas;
import scid.storage.array;
import scid.storage.cowarray;
import scid.common.traits, scid.common.meta;
import scid.vector;
import std.traits, std.range, std.algorithm;

template ArrayViewStorage( ElementOrArray, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias BasicArrayViewStorage!( CowArrayRef!ElementOrArray, ArrayViewType.Interval, vectorType ) ArrayViewStorage;
	else
		alias BasicArrayViewStorage!( ElementOrArray, ArrayViewType.Interval, vectorType )             ArrayViewStorage;
}

template StridedArrayViewStorage( ElementOrArray, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias BasicArrayViewStorage!( CowArrayRef!ElementOrArray, ArrayViewType.Strided, vectorType ) StridedArrayViewStorage;
	else
		alias BasicArrayViewStorage!( ElementOrArray, ArrayViewType.Strided, vectorType )             StridedArrayViewStorage;
}

enum ArrayViewType {
	Interval,
	Strided
}

struct BasicArrayViewStorage( ArrayRef_, ArrayViewType strided_, VectorType vectorType_ ) {
	alias ArrayRef_                                                                                 ArrayRef;
	alias BaseElementType!ArrayRef                                                                  ElementType;	
	alias typeof( this )                                                                            View;
	alias vectorType_                                                                               vectorType;
	// alias BasicArrayStorage!( ArrayRef, vectorType )                                                Referenced;
	alias BasicArrayViewStorage!( ArrayRef, ArrayViewType.Strided, vectorType )                     StridedView;
	alias BasicArrayViewStorage!( ArrayRef, ArrayViewType.Strided, transposeVectorType!vectorType ) Transposed;
	
	enum isStrided = (strided_ == ArrayViewType.Strided );
	
	static if( isStrided ) {	
		this( ref ArrayRef arrayRef, size_t iFirstIndex, size_t iLength, size_t iStride )
		in {
			assert( iStride != 0, strideMsg_ );
		} body {
			array_ = arrayRef;
			setParams_( iFirstIndex, iLength, iStride );
		}
	} else {
		this( ref ArrayRef arrayRef, size_t iFirstIndex, size_t iLength ) {
			array_ = arrayRef;
			setParams_( iFirstIndex, iLength );
		}
	}
	
	void forceRefSharing( ref typeof(this) rhs ){
		this = rhs;
	}
	
	ref typeof( this ) opAssign( typeof( this ) rhs ) {
		firstIndex_ = rhs.firstIndex;
		length_     = rhs.length_;
		move( rhs.array_, array_ );
		static if( isStrided )
			stride_ = rhs.stride_;
		return this;
	}
	
	ElementType index( size_t i ) const
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		return array_.cdata[ map_( i ) ];
	}
	
	void indexAssign( string op="" )( ElementType rhs, size_t i )
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		mixin("array_.data[ map_( i ) ] " ~ op ~ "= rhs;");
	}
	
	typeof(this) slice( size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
	} body {
		static if( isStrided )
			return typeof( this )( array_, map_(start), end - start, stride_ );
		else
			return typeof( this )( array_, map_(start), end - start );
	}
	
	View view( size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_(start, end) );
	} body {
		static if( isStrided )
			return typeof( return )( array_, start, end-start, stride_ );
		else
			return typeof( return )( array_, start, end-start );
	}
	
	StridedView view( size_t start, size_t end, size_t newStride )
	in {
		assert( newStride != 0, strideMsg_ );
		assert( start < end && end <= length, sliceMsg_(start, end) );
	} body {
		size_t len = (end - start);
		len = len / newStride + ( (len % newStride) != 0 );
		static if( isStrided )
			newStride *= stride_;
		return typeof( return )( array_, start, len, newStride );
	}
	
	void resizeOrClear( size_t rlength ) {
		resizeOrClear( rlength, null );
		blas.scal( length_, Zero!ElementType, data, stride );
	}
	
	void resizeOrClear( size_t rlength, void* ) {
		assert( length == rlength, "Length mismatch in vector operation." );
	}
	
	void copy( Transpose tr, S )( ref S rhs ) if( isStridedStorage!(S, ElementType) ) {
		hlStridedCopy( rhs, this );
	}
	
	mixin StridedAxpyScalDot;
	
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
	
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ ": popBack() on empty." );
	} body {
		-- length_;
	}
	
	@property {
		ref ArrayRef         array()            { return array_; }
		ElementType*         data()             { return array_.data + firstIndex_; }
		const(ElementType)*  cdata()      const { return array_.cdata + firstIndex_; }
		size_t               length()     const { return length_; }
		bool                 empty()      const { return length_ == 0; }
		size_t               firstIndex() const { return firstIndex_; }
		
		void front( ElementType newValue )
		in {
			assert( !empty, msgPrefix_ ~ "front assign on empty." );
		} body {
			indexAssign( newValue, 0 );
		}
		
		void back( ElementType newValue  )
		in {
			assert( !empty, msgPrefix_ ~ "back assign on empty." );
		} body {
			indexAssign( newValue, length - 1 );
		}
			
		ElementType front() const
		in {
			assert( !empty, msgPrefix_ ~ "front get on empty." );
		} body {
			return this.index( 0 );
		}
		
		ElementType back() const
		in {
			assert( !empty, msgPrefix_ ~ "back get on empty." );
		} body {
			return index( length_ - 1 );
		}
		
		static if( isStrided ) {
			size_t stride() const { return stride_; }
		} else {
			enum stride = 1;
		}
	}
	
private:
	mixin SliceIndexMessages;

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
	
	size_t   firstIndex_, length_;
	ArrayRef array_;
}



unittest {
	// TODO: Write tests for ArrayViewStorage.
}