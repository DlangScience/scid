module scid.storage.arrayview;

import scid.bindings.blas.dblas;

import scid.internal.assertmessages;
import scid.storage.cowarray;
import scid.common.traits;
import std.traits, std.range, std.algorithm;



template ArrayViewStorage( ElementOrArray )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias BasicArrayViewStorage!( CowArrayRef!ElementOrArray, ArrayViewType.Interval ) ArrayViewStorage;
	else
		alias BasicArrayViewStorage!( ElementOrArray, ArrayViewType.Interval )             ArrayViewStorage;
}

template StridedArrayViewStorage( ElementOrArray )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias BasicArrayViewStorage!( CowArrayRef!ElementOrArray, ArrayViewType.Strided ) StridedArrayViewStorage;
	else
		alias BasicArrayViewStorage!( ElementOrArray, ArrayViewType.Strided )             StridedArrayViewStorage;
}

enum ArrayViewType {
	Interval,
	Strided
}

struct BasicArrayViewStorage( ArrayRef_, ArrayViewType strided_ ) {
	alias ArrayRef_                                                 ArrayRef;
	alias BaseElementType!ArrayRef                                  ElementType;	
	alias typeof( this )                                            View;
	alias BasicArrayViewStorage!( ArrayRef, ArrayViewType.Strided ) StridedView;
	
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
	
	void copyLeft( S )( ref S rhs ) 
	in {
		static assert( hasLength!S );
		assert( rhs.length == length, sliceAssignMsg_( 0, length, rhs.length ) );
	} body {
		static if( is( S == ElementType[] ) ) {
			// assignment to built-in array
			copy( length_, rhs.ptr, 1, array_.data.ptr + firstIndex_, stride );
		} else static if( isConcreteStorageOf!(S, ElementType) ) {
			// assignment to another view or ArrayStorage
			static if( isStrided || (!is( S == ArrayStorage!ArrayRef ) && !is( S == View )) ) {
				// view = view / view = storage
				array_ = rhs.array_;
				setParams(rhs.firstIndex, rhs.length);
			} else {
				copy( length_, rhs.cdata.ptr, rhs.stride, array_.data.ptr + firstIndex_, stride );
			}
		} else static if( isInputRange!S && hasLength!S ) {
			// assignment to a range
			T[] d = data;
			foreach( x ; rhs ) {
				d[ 0 ] = x;
				d = d[ stride .. $ ];
			}
		} else
			static assert( false, "Invalid type '" ~ typeof(rhs).stringof ~ "' for assignment to '" ~ typeof(this).stringof ~ "'." );
	}
	
	void slicedCopyLeft( S )( ref S rhs, size_t start, size_t end )
	in {
		static assert( hasLength!S );
		assert( start < end && end <= length, sliceMsg_( start, end ) );
		assert( end-start == rhs.length, sliceAssignMsg_( start, end, rhs.length )  );
	} body {
		size_t n = end - start;
		static if( is( S == ElementType[] ) ) {
			// assignment to built-in array
			copy( n, rhs.ptr, 1, array_.data.ptr + map_(start), stride );
		} else static if( isConcreteStorageOf!(S, ElementType) ) {
			// assignment to any concrete storage
			copy( n, rhs.cdata.ptr, rhs.stride, array_.data.ptr + map_(start), stride );
		} else static if( isInputRange!S && hasLength!S ) {
			// assignment to a range
			T[] d = array_.data[ map_(start) .. map_(end) ];
			foreach( x ; rhs ) {
				d[ 0 ] = x;
				d = d[ stride .. $ ];
			}
		} else {
			// invalid assignment
			static assert( false, "Invalid type '" ~ typeof(rhs).stringof ~ "' for assignment to '" ~ typeof(this).stringof ~ "'." );
		}
	}
	
	void axpyLeft( S )( ref S rhs, ElementType alpha )
	in {
		static assert( hasLength!S );
		assert( rhs.length == length, fmt( msgPrefix_ ~ "axpy length mismatch %d vs %d", length, rhs.length ) ); 
	} body {
		static if( is( S == ElementType[] ) ) {
			// built-in array
			axpy( length, alpha, rhs.ptr, 1, array_.data.ptr + firstIndex_, stride );
		} else static if( isConcreteStorageOf!(S, ElementType) ) {
			// any concrete storage
			axpy( length, alpha, rhs.cdata.ptr, rhs.stride, array_.data.ptr + firstIndex_, stride );
		} else static if( isInputRange!S && hasLength!S ) {
			// input range
			T[] d = data;
			foreach( x ; rhs ) {
				d[ 0 ] += x * alpha;
				d = d[ stride .. $ ];
			}
		} else
			static assert( false, "Invalid type '" ~ typeof(rhs).stringof ~ "' for axpy to '" ~ typeof(this).stringof ~ "'." );
	}
	
	void slicedAxpyLeft( S )( ref S rhs, ElementType alpha, size_t start, size_t end )
	in {
		static assert( hasLength!S );
		assert( start < end && end <= length, sliceMsg_( start, end ) );
		assert( end-start == rhs.length, sliceAssignMsg_( start, end, rhs.length )  );
	} body {
		size_t n = end - start;
		static if( is( S == ElementType[] ) ) {
			// built-in array
			axpy( n, alpha, rhs.ptr, 1, array_.data.ptr + map_(start), stride );
		} else static if( isConcreteStorageOf!(S, ElementType) ) {
			// any concrete storage
			axpy( n, alpha, rhs.cdata.ptr, rhs.stride, array_.data.ptr + map_(start), stride );
		} else static if( isInputRange!S && hasLength!S ) {
			// input range
			T[] d = array_.data[ map_(start) .. map_(end) ];
			foreach( x ; rhs ) {
				d[ 0 ] += x * alpha;
				d = d[ stride .. $ ];
			}
		} else
			static assert( false, "Invalid type '" ~ typeof(rhs).stringof ~ "' for axpy to '" ~ typeof(this).stringof ~ "'." );
	}
	
	void scalLeft( ElementType rhs ) {
		scal( length, rhs, array.data.ptr + firstIndex_, stride );	
	}
	
	void slicedScalLeft( ElementType rhs, size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
	} body {
		scal( end-start, rhs, array.data.ptr + map_(start), stride );
	}
	
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
		ArrayRef              array()            { return array_; }
		ElementType[]         data()             { return array_.data[ firstIndex_ .. map_(length_) ]; }
		const(ElementType[])  cdata()      const { return array_.cdata[ firstIndex_ .. map_(length_) ]; }
		size_t                length()     const { return length_; }
		bool                  empty()      const { return length_ == 0; }
		size_t                firstIndex() const { return firstIndex_; }
		
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