module scid.storage.packedsubvec;

import scid.internal.assertmessages;
import scid.vector, scid.matrix;
import scid.common.traits;
import std.algorithm, std.traits, std.range;

struct PackedSubVectorStorage( MatrixRef_, VectorType vtype_ ) {
	alias vtype_                    vectorType;
	alias MatrixRef_                MatrixRef;
	alias BaseElementType!MatrixRef ElementType;
	alias typeof(this)              Slice;
	alias typeof(this)              View;
	alias typeof(this)              StridedView;
	
	enum bool isRow   = ( vectorType       == VectorType.Row );
	enum bool isUpper = ( matrix_.triangle == MatrixTriangle.Upper );
	
	this( ref MatrixRef matrixRef, size_t fixed, size_t start, size_t len ) {
		matrix_ = matrixRef;
		fixed_  = fixed;
		start_  = start;
		length_ = len;
	}
	
	void forceRefAssign( ref typeof(this) rhs ) {
		this = rhs;	
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		move( rhs.matrix_, matrix_ );	
		fixed_  = rhs.fixed_;
		start_  = rhs.start_;
		length_ = rhs.length_;
		
		return this;
	}
	
	Slice slice( size_t start, size_t end )
	in {
		assert( start < end && end <= length_, sliceMsg_( start, end ) );
	} body {
		return typeof( return )( matrix_, fixed_, start + start_, end - start );	
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
		static if( isRow ) return matrix_.index( fixed_, i + start_ );
		else               return matrix_.index( i + start_, fixed_ );
	}
	
	void indexAssign( string op="" )( ElementType rhs, size_t i )
	in {
		assert( i < length, boundsMsg_( i ) );
	} body {
		static if( isRow ) matrix_.indexAssign!op( rhs, fixed_, i + start_ );
		else               matrix_.indexAssign!op( rhs, i + start_, fixed_ );
	}
	
	void copyLeft( S )( S rhs ) if( isInputRange!(Unqual!S) && hasLength!(Unqual!S) )
	in {
		assert( rhs.length == length, sliceAssignMsg_( 0, length, rhs.length ) );
	} body {
		slicedCopyLeft( rhs, 0, length_ );
	}
	
	void slicedCopyLeft( S )( ref S rhs, size_t start, size_t end ) if( isInputRange!(Unqual!S) && hasLength!(Unqual!S) )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
		assert( end-start == rhs.length, sliceAssignMsg_( start, end, rhs.length )  );
	} body {
		size_t i = realStart_( start );
		size_t e = realEnd_( end );
		popFrontN(rhs,  i );
		for( ; i < e ; ++i ) {
			indexAssign( rhs.front, i );
			rhs.popFront();
		}
	}
	
	void axpyLeft( S )( ref S rhs, ElementType alpha ) if( isInputRange!(Unqual!S) && hasLength!(Unqual!S) )
	in {
		assert( rhs.length <= length, fmt( msgPrefix_ ~ "axpy length mismatch %d vs %d", length, rhs.length ) ); 
	} body {
		slicedAxpyLeft( rhs, alpha, 0, length_ );
	}
	
	void slicedAxpyLeft( S )( ref S rhs, ElementType alpha, size_t start, size_t end ) if( isInputRange!(Unqual!S) && hasLength!(Unqual!S) )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
		assert( end-start >= rhs.length, sliceAssignMsg_( start, end, rhs.length )  );
	} body {
		size_t i = realStart_( start );
		size_t e = realEnd_( end );
		popFrontN( rhs, i );
		for( ; i < e ; ++i ) {
			indexAssign!"+"( rhs.front * x, i );
			rhs.popFront();
		}
	}
	
	void scalLeft( ElementType rhs ) {
		slicedScalLeft( rhs, 0, length_ );
	}
	
	void slicedScalLeft( ElementType rhs, size_t start, size_t end )
	in {
		assert( start < end && end <= length, sliceMsg_( start, end ) );
	} body {
		size_t i = realStart_( start );
		size_t e = realEnd_( end );
		for( ; i < e ; ++i )
			indexAssign!"*"( rhs, i );
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
		ElementType[]        data()         { return matrix_.data; }
		const(ElementType[]) cdata()  const { return matrix_.cdata; }
		MatrixRef_           matrix()       { return matrix_; }
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
	
private:
	mixin SliceIndexMessages;

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

	MatrixRef matrix_;
	size_t    fixed_, start_, length_;
}
