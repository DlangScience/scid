module scid.storage.packedsubvec;

import scid.internal.assertmessages;
import scid.ops.common;
import scid.vector, scid.matrix;
import scid.common.storagetraits;
import scid.common.meta;
import scid.ops.eval;
import std.algorithm, std.traits, std.range, std.conv;
import scid.storage.array;


struct PackedSubVectorStorage( ContainerRef_, VectorType vtype_ ) {
	alias vtype_                                                           vectorType;
	alias ContainerRef_                                                       ContainerRef;
	alias BaseElementType!ContainerRef                                        ElementType;
	alias PackedSubVectorStorage!( ContainerRef, transposeVectorType!vtype_ ) Transposed;
	// alias ArrayStorage!(ArrayTypeOf!ContainerRef)                             Referenced
	alias typeof(this)                                                     Slice;
	alias typeof(this)                                                     View;
	alias typeof(this)                                                     StridedView;
	
	enum bool isRow   = ( vectorType       == VectorType.Row );
	enum bool isUpper = ( containerRef_.triangle == MatrixTriangle.Upper );
	
	this( ref ContainerRef containerRef, size_t fixed, size_t start, size_t len ) {
		containerRef_ = containerRef;
		fixed_  = fixed;
		start_  = start;
		length_ = len;
	}
	
	void forceRefAssign( ref typeof(this) rhs ) {
		this = rhs;	
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		move( rhs.containerRef_, containerRef_ );	
		fixed_  = rhs.fixed_;
		start_  = rhs.start_;
		length_ = rhs.length_;
		
		return this;
	}
	
	void resize( size_t newLength ) {
		resize( newLength, null );
		scale( Zero!ElementType );
	}
	
	void resize( size_t newLength, void* ) {
		assert( newLength == length_, msgPrefix_ ~ "Cannot change length of packed subvector." );
	}
	
	Slice slice( size_t start, size_t end )
	in {
		assert( start < end && end <= length_, sliceMsg_( start, end ) );
	} body {
		return typeof( return )( containerRef_, fixed_, start + start_, end - start );	
	}
	
	View view( size_t start, size_t end ) {
		return slice( start, end );
	}
	
	View view( size_t start, size_t end, size_t stride ) {
		assert( false, "Strided views of packed sub vectors is not implemented." );
	}
	
	ElementType index( size_t i ) const
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		static if( isRow ) return containerRef_.index( fixed_, i + start_ );
		else               return containerRef_.index( i + start_, fixed_ );
	}
	
	void indexAssign( string op="" )( ElementType rhs, size_t i )
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		static if( isRow ) containerRef_.indexAssign!op( rhs, fixed_, i + start_ );
		else               containerRef_.indexAssign!op( rhs, i + start_, fixed_ );
	}
	
	void copy( Transpose tr = Transpose.yes, S )( auto ref S rhs ) if( isInputRange!(Unqual!S) && hasLength!(Unqual!S) )
	in {
		assert( rhs.length == length, sliceAssignMsg_( 0, length, rhs.length ) );
	} body {
		slicedCopy_!tr( rhs, 0, length_ );
	}
	
	void scaledAddition( Transpose tr = Transpose.yes, S )( ElementType alpha, auto ref S rhs ) if( isInputRange!(Unqual!S) && hasLength!(Unqual!S) )
	in {
		assert( rhs.length <= length, msgPrefix_ ~ "axpy length mismatch " ~ to!string(length) ~ " vs. " ~ to!string(rhs.length) ); 
	} body {
		slicedScaledAddition_!tr( alpha, rhs, 0, length_ );
	}
	
	void scale( ElementType rhs ) {
		slicedScale_( rhs, 0, length_ );
	}
	
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront on empty." );
	} body {
		++ start_; -- length_;
	}
	
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popBack on empty." );
	} body {
		-- length_;
	}
	
	@property {
		ref ContainerRef        matrix()       { return containerRef_; }
		size_t               start()  const { return start_; }
		size_t               fixed()  const { return fixed_; }
		size_t               length() const { return length_; }
		bool                 empty()  const { return length_ == 0; }
		
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
			indexAssign( newValue, length_ - 1 );
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
	}
	
	/** Promotions for this type are inherited from ArrayStorage */
	template Promote( Other ) {
		alias Promotion!( BasicArrayStorage!(ArrayTypeOf!ContainerRef, vectorType), Other ) Promote;
	}
	
private:
	mixin ArrayErrorMessages;

	void slicedCopy_( Transpose tr = Transpose.yes, S )( auto ref S rhs, size_t start, size_t end ) if( isInputRange!(Unqual!S) && hasLength!(Unqual!S) )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
		assert( end-start == rhs.length, sliceAssignMsg_( start, end, rhs.length )  );
	} body {
		size_t i = realStart_( start );
		size_t e = realEnd_( end );
		popFrontN(rhs,  i );
		for( ; i < e ; ++i ) {
			static if( tr && isComplex!ElementType )
				indexAssign( gconj(rhs.front), i );
			else
				indexAssign( rhs.front, i );
			rhs.popFront();
		}
	}

	void slicedScale_( ElementType rhs, size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
	} body {
		size_t i = realStart_( start );
		size_t e = realEnd_( end );
		for( ; i < e ; ++i )
			indexAssign!"*"( rhs, i );
	}
	
	void slicedScaledAddition_( Transpose tr = Transpose.yes, S )( ElementType alpha, auto ref S rhs, size_t start, size_t end ) if( isInputRange!(Unqual!S) && hasLength!(Unqual!S) )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
		assert( end-start >= rhs.length, sliceAssignMsg_( start, end, rhs.length )  );
	} body {
		size_t i = realStart_( start );
		size_t e = realEnd_( end );
		popFrontN( rhs, i );
		for( ; i < e ; ++i ) {
			static if( tr && isComplex!ElementType )
				indexAssign!"+"( gconj(rhs.front) * alpha, i );
			else
				indexAssign!"+"( rhs.front * alpha, i );
			rhs.popFront();
		}
	}

	size_t realStart_( size_t fakeStart ) {
		static if( !isUpper ^ !isRow ) return fakeStart;
		else {
			if( fixed_ < (start_ + fakeStart) )
				return fakeStart;
			else
				return fixed_ - start_;
		}
	}
	
	size_t realEnd_( size_t fakeEnd ) {
		static if( isUpper ^ !isRow ) return fakeEnd;
		else {
			if( start_ <= fixed_ )
				return min( fakeEnd, fixed_ - start_ + 1 );
			else
				return 0;
		}
	}

	ContainerRef containerRef_;
	size_t    fixed_, start_, length_;
}
