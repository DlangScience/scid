module scid.storage.cowarray;

import scid.internal.assertmessages;
import scid.storage.arraydata;
import scid.bindings.blas.dblas;
import scid.common.meta;

import std.typecons, std.algorithm, std.array;
import std.conv;

version(unittest) {
    //import scid.common.testing;
	import std.stdio;
}

/** Copy-on-write array implementation. Wrapped by Vectors & Matrices for block storage. */
struct CowArray( ElementType_ ) {
	alias ElementType_          ElementType;
	alias ArrayData!ElementType Data;
	
	this( size_t newLength ) {
		this( newLength, null );
		scal( newLength, Zero!ElementType, ptr_, 1 );
	}
	
	this( size_t newLength, void* ) {
		data_.reset( newLength );
		ptr_    = data_.ptr;
		length_ = newLength;
	}
	
	this( ref Data data, ElementType* ptr, size_t length ) {
		data_   = data;
		ptr_    = ptr;
		length_ = length;
	}
	
	this( ElementType[] initializer ) {
		data_.reset( initializer );
		ptr_    = data_.ptr;
		length_ = initializer.length;
	}
	
	this( typeof(this)* other ) {
		data_   = other.data_;
		ptr_    = other.ptr_;
		length_ = other.length_;
	}
	
	void resizeOrClear( size_t len ) {
		resizeOrClear( len, null );
		scal( len, Zero!ElementType, ptr_, 1 );
	}
	
	void resizeOrClear( size_t len, void* ) {
		if( len != length ) {
			data_.reset( len );
			ptr_ = data_.ptr;
			length_ = len;
		}
	}
	
	ref typeof(this) opAssign( CowArray rhs ) {
		move( rhs.data_, data_ );
		ptr_    = rhs.ptr_;
		length_ = rhs.length_;
		return this;
	}
	
	ElementType index( size_t i ) const
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		return ptr_[ i ];
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i )
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		unshareData_();
		mixin( "ptr_[ i ]" ~ op ~ "= rhs;" );
	}
	
	typeof( this ) slice( size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
	} body {
		CowArray r;
		r.data_ = data_;
		r.length_ = end - start;
		r.ptr_    = ptr_ + start;
		return r;
	}
	
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront() on empty." );
	} body {
		++ ptr_;
		-- length_;
	}
	
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popBack() on empty." );
	} body {
		-- length_;
	}
	
	@property {
		ElementType*         data()         { unshareData_(); return ptr_; }
		const(ElementType)*  cdata()  const { return ptr_; }
		size_t               length() const { return length_; }
		bool                 empty()  const { return length_ == 0; }
		typeof(this)*        ptr()          { return &this; }
		
		void front( ElementType newValue )
		in {
			assert( !empty, msgPrefix_ ~ "front assign on empty." );
		} body {
			unshareData_();
			*ptr_ = newValue;
		}
		
		void back( ElementType newValue  )
		in {
			assert( !empty, msgPrefix_ ~ "back assign on empty." );
		} body {
			unshareData_();
			*(ptr_ + length_ - 1) = newValue;
		}
			
		ElementType front() const
		in {
			assert( !empty, msgPrefix_ ~ "front get on empty." );
		} body {
			return *ptr_;
		}
		
		ElementType back() const
		in {
			assert( !empty, msgPrefix_ ~ "back get on empty." );
		} body {
			return *(ptr_ + length_ - 1);
		}
	}
	
private:
	mixin SliceIndexMessages;
	
	void unshareData_() {
		if( data_.refCount() == 1 )
			return;
		
		if( ptr_ == data_.ptr && length_ == data_.length )
			data_.unshare();
		else
			data_.reset( ptr_[ 0 .. length_ ] );
		
		ptr_ = data_.ptr;
	}
	
	Data         data_;
	ElementType* ptr_;
	size_t       length_;
}

template CowArrayRef( T ) {
	alias RefCounted!(CowArray!T, RefCountedAutoInitialize.no )
			CowArrayRef;
}

unittest {
	// TODO: CowArray tests.
}