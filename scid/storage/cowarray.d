module scid.storage.cowarray;

import scid.internal.assertmessages;
import scid.storage.arraydata;

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
	
	this( size_t newLength, ElementType initWith = ElementType.init ) {
		this( newLength, null );
		slice_[] = initWith;
	}
	
	this( size_t newLength, void* ) {
		data_.reset( newLength );
		slice_ = data_.array;
	}
	
	this( ref Data data, ElementType[] slice ) {
		data_  = data;
		slice_ = slice;
	}
	
	this( ElementType[] initializer ) {
		data_.reset( initializer );
		slice_ = data_.array;
	}
	
	this( typeof(this)* other ) {
		data_  = other.data_;
		slice_ = data_.array;
	}
	
	ref typeof(this) opAssign( CowArray rhs ) {
		move( rhs.data_, data_ );
		return this;
	}
	
	ElementType index( size_t i ) const
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		return slice_[ i ];
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i )
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		unshareData_();
		mixin( "slice_[ i ]" ~ op ~ "= rhs;" );
	}
	
	typeof( this ) slice( size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
	} body {
		CowArray r;
		r.data_  = data_;
		r.slice_ = slice_[ start .. end ];
		
		return r;
	}
	
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront() on empty." );
	} body {
		slice_.popFront();	
	}
	
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popBack() on empty." );
	} body {
		slice_.popFront();	
	}
	
	string toString() const {
		return to!string( cdata );
	}
	
	@property {
		ElementType[]         data()         { unshareData_(); return slice_; }
		const(ElementType[])  cdata()  const { return slice_; }
		size_t                length() const { return slice_.length; }
		bool                  empty()  const { return slice_.empty; }
		typeof(this)*         ptr()          { return &this; }
		
		void front( ElementType newValue )
		in {
			assert( !empty, msgPrefix_ ~ "front assign on empty." );
		} body {
			unshareData_();
			slice_.front = newValue;
		}
		
		void back( ElementType newValue  )
		in {
			assert( !empty, msgPrefix_ ~ "back assign on empty." );
		} body {
			unshareData_();
			slice_.back = newValue;
		}
			
		ElementType front() const
		in {
			assert( !empty, msgPrefix_ ~ "front get on empty." );
		} body {
			return slice_.front;
		}
		
		ElementType back() const
		in {
			assert( !empty, msgPrefix_ ~ "back get on empty." );
		} body {
			return slice_.back;
		}
	}
	
private:
	mixin SliceIndexMessages;
	
	void unshareData_() {
		if( data_.refCount() == 1 )
			return;
		
		if( slice_ == data_.array )
			data_.unshare();
		else
			data_.reset( slice_ );
		
		slice_ = data_.array;
	}
	
	Data          data_;
	ElementType[] slice_;
}

template CowArrayRef( T ) {
	alias RefCounted!(CowArray!T, RefCountedAutoInitialize.yes )
			CowArrayRef;
}

unittest {
	// TODO: CowArray tests.
}