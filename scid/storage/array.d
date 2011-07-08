module scid.storage.array;

import scid.bindings.blas.dblas;

import scid.internal.assertmessages;
import scid.storage.cowarray;
import scid.storage.arrayview;
import scid.common.traits;
import std.traits, std.range, std.algorithm;

template isConcreteStorageOf( S, E ) {
	static if( is( S == struct ) &&                          // it's a struct
			   is( typeof(S.init.stride) : size_t ) &&       // it has a stride
			   is( typeof(S.init.cdata) : const(E[]) ) &&    // it has const access to the data
			   is( typeof(S.init.length) : size_t ) )        // it has a length
		enum isConcreteStorageOf = true;
	else
		enum isConcreteStorageOf = false;
}

template ArrayStorage( ElementOrArray )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias BasicArrayStorage!( CowArrayRef!ElementOrArray ) ArrayStorage;
	else
		alias BasicArrayStorage!( ElementOrArray )             ArrayStorage;
}

struct BasicArrayStorage( ArrayRef_ ) {
	alias ArrayRef_                ArrayRef;
	alias BaseElementType!ArrayRef ElementType;	
	
	alias ArrayView!( ArrayRef )        View;
	alias StridedArrayView!( ArrayRef ) StridedView;
	alias array_                        this;
	
	enum  stride     = 1;
	enum  firstIndex = 0;
	
	this( size_t newLength, ElementType initWith = ElementType.init ) {
		array_ = ArrayRef( newLength, initWith );	
	}
	
	this( size_t newLength, void* ) {
		array_ = ArrayRef( newLength, null );
	}
	
	this( ElementType[] initializer ) {
		array_ = ArrayRef( initializer );
	}
	
	this( this ) {
		array_ = ArrayRef( array_.ptr );
	}
	
	void forceRefAssign( ref typeof(this) rhs ) {
		array_ = rhs.array;	
	}
		
	ref typeof(this) opAssign( typeof(this) rhs ) {
		array_ = ArrayRef( rhs.array_.ptr );
		return this;
	}
	
	typeof( this ) slice( size_t start, size_t end ) {
		return typeof(this)( ArrayRef(array_.slice( start, end )) );
	}
	
	View view( size_t start, size_t end ) {
		return typeof( return )( array_, start, end - start );
	}
	
	StridedView view( size_t start, size_t end, size_t stride ) {
		size_t len = (end - start);
		len = len / stride + (len % stride != 0);
		return typeof( return )( array_, start, len, stride );
	}
	
	@property {
		ArrayRef     array()       { return array_; }
	}
	
	// Workaround: alias this doesn't make this a range //
	@property {
		bool         empty()  const { return array_.empty; }
		ElementType  front()  const { return array_.front; }
		ElementType  back()   const { return array_.back; }
		size_t       length() const { return array_.length; }
	}
	
	void popFront() { array_.popFront(); }
	void popBack()  { array_.popBack(); }
	// Workaround ends. //
	
	void copyLeft( S )( ref S rhs ) {
		static if( is( S == ElementType[] ) ) {
			// assignment to built-in array
			auto n = rhs.length;
			if( n != length ) {
				array_ = ArrayRef( rhs );
			} else {
				copy( n, rhs.ptr, 1, array_.data.ptr, 1 );
			}
		} else static if( is( S == typeof( this ) ) ) {
			// assignment to another ArrayStorage
			this = rhs;
		} else static if( is( S == View ) )  {
			// assignment to a contiguous view of an ArrayStorage
			this = rhs.array.slice( rhs.firstIndex, rhs.firstIndex + rhs.length );
		} else static if( isConcreteStorageOf!(S, ElementType) ) {
			// assignment to any concrete storage
			auto n = rhs.length;
			if( n != length ) 
				array_ = ArrayRef( n, null );
			copy( n, rhs.cdata.ptr, rhs.stride, data.ptr, 1 );
		} else static if( isInputRange!S && hasLength!S ) {
			// assignment to a range
			auto n = rhs.length;
			
			if( n != length ) 
				array_ = ArrayRef( n, null );
			
			S rangeCopy = rhs;
			foreach( ref x; array_.data ) {
				x = rangeCopy.front;
				rangeCopy.popFront();
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
			copy( n, rhs.ptr, 1, array_.data.ptr, 1 );
		} else static if( isConcreteStorageOf!(S, ElementType) ) {
			copy( n, rhs.cdata.ptr, rhs.stride, array_.data.ptr + start, 1 );
		} else static if( isInputRange!S && hasLength!S ) {
			// assignment to a range
			auto n = rhs.length;
			
			if( n != length ) 
				array_ = ArrayRef( n, null );
			
			S rangeCopy = rhs;
			foreach( ref x; array_.data[ start .. end ] ) {
				x = rangeCopy.front;
				rangeCopy.popFront();
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
			axpy( length, alpha, rhs.ptr, 1, array_.data.ptr, 1 );
		} else static if( isConcreteStorageOf!(S, ElementType) ) {
			// any concrete storage
			axpy( length, alpha, rhs.cdata.ptr, rhs.stride, array_.data.ptr, 1 );
		} else static if( isInputRange!S && hasLength!S ) {
			// input range
			S rangeCopy = rhs;
			foreach( ref x; array_.data ) {
				x += rangeCopy.front * alpha;
				rangeCopy.popFront();
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
			axpy( n, alpha, rhs.ptr, 1, array_.data.ptr + start, 1 );
		} else static if( isConcreteStorageOf!(S, ElementType) ) {
			// any concrete storage
			axpy( n, alpha, rhs.cdata.ptr, rhs.stride, array_.data.ptr + start, 1 );
		} else static if( isInputRange!S && hasLength!S ) {
			// input range
			S rangeCopy = rhs;
			foreach( ref x; array_.data[ start .. end ] ) {
				x += rangeCopy.front * alpha;
				rangeCopy.popFront();
			}
		} else
			static assert( false, "Invalid type '" ~ typeof(rhs).stringof ~ "' for axpy to '" ~ typeof(this).stringof ~ "'." );
	}
	
	void scalLeft( ElementType rhs ) {
		scal( length, rhs, array.data.ptr, 1 );	
	}
	
	void slicedScalLeft( ElementType rhs, size_t start, size_t end ) {
		scal( end-start, rhs, array.data.ptr + start, 1 );
	}
	
	ArrayRef array_;
	
private:
	mixin SliceIndexMessages;

	this( ArrayRef arrayRef ) {
		array_ = arrayRef;
	}
}

unittest {
	// TODO: Write tests for ArrayStorage.
}