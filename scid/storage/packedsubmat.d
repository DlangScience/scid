module scid.storage.packedsubmat;

import scid.common.traits;
import scid.matrix, scid.vector;
import scid.internal.assertmessages;
import scid.storage.packedsubvec;
import std.algorithm;

enum SubMatrixType {
	View,
	Slice
}

struct PackedSubMatrixStorage( MatrixRef_, SubMatrixType type_ ) {
	alias MatrixRef_                                                         MatrixRef;
	alias BaseElementType!MatrixRef                                          ElementType;
	alias storageOrderOf!MatrixRef                                           storageOrder;
	alias typeof(this)                                                       Slice;
	alias typeof(this)                                                       View;
	alias Vector!( PackedSubVectorStorage!( MatrixRef, VectorType.Row ) )    RowView;
	alias Vector!( PackedSubVectorStorage!( MatrixRef, VectorType.Column ) ) ColumnView;
	alias ColumnView                                                         DiagonalView;
	
	enum isView     = type_ == SubMatrixType.View;
	enum isRowMajor = (storageOrder == StorageOrder.RowMajor );
	
	this( ref MatrixRef matrixRef, size_t rowStart, size_t numRows, size_t colStart, size_t numCols ) {
		assignMatrix_( matrixRef );
		rowStart_ = rowStart; rows_ = numRows;
		colStart_ = colStart; cols_ = numCols;
	}
	
	static if( !isView ) {
		this( this ) {
			
			assignMatrix_( matrix_ );
		}
	}
		
	void forceRefAssign( ref typeof(this) rhs ) {
		matrix_ = rhs.matrix_;
		rowStart_ = rhs.rowStart_; rows_ = rhs.rows_;
		colStart_ = rhs.colStart_; cols_ = rhs.cols_;
	}
	
	ref typeof( this ) opAssign( typeof(this) rhs ) {
		move( rhs.matrix_, matrix_ );
		rowStart_ = rhs.rowStart_; rows_ = rhs.rows_;
		colStart_ = rhs.colStart_; cols_ = rhs.cols_;
		return this;
	}
	
	ElementType index( size_t i, size_t j ) const 
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		return matrix_.index( i + rowStart_, j + colStart_ );
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		matrix_.indexAssign!op( rhs, i + rowStart_, j + colStart_ );
	}
	
	RowView row( size_t i )
	in {
		assert( i < rows_, sliceMsg_(i,0,i,cols_) );
	} body {
		return typeof( return )( matrix_, i + rowStart_, colStart_, cols_ );
	}
	
	ColumnView column( size_t j )
	in {
		assert( j < cols_, sliceMsg_(0,j,rows_,j) );
	} body {
		return typeof( return )( matrix_, j + colStart_, rowStart_, rows_ );
	}
	
	RowView rowSlice( size_t i, size_t start, size_t end )
	in {
		assert( i < rows_ && start < end && end <= cols_, sliceMsg_(i,start,i,end) );
	} body {
		return typeof( return )( matrix_, i + rowStart_, start + colStart_, end - start );
	}
	
	ColumnView columnSlice( size_t j, size_t start, size_t end )
	in {
		assert( j < cols_ && start < end && end <= rows_, sliceMsg_(start,j,end,j) );
	} body {
		return typeof( return )( matrix_, j + colStart_, start + rowStart_, end - start );
	}
	
	View view( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( matrix_, rowStart + rowStart_, rowEnd - rowStart, colStart + colStart_, colEnd - colStart );
	}
	
	Slice slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( matrix_, rowStart + rowStart_, rowEnd - rowStart, colStart + colStart_, colEnd - colStart );
	}
	
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront on empty." );
	} body {
		static if( isRowMajor ) { ++ rowStart_; -- rows_; }
		else                    { ++ colStart_; -- cols_; }
	}
	
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popBack on empty." );
	} body {
		-- major_;	
	}
	
	static if( isRowMajor ) {
		alias rows    major;
		alias columns minor;
		private {
			alias rows_ major_;
			alias cols_ minor_;
		}
	} else {
		alias rows    minor;
		alias columns major;
		private {
			alias rows_ minor_;
			alias cols_ major_;
		}
	}
	
	@property {
		MatrixRef            matrix()        { return matrix_; }
		ElementType[]        data()          { return matrix_.data; }
		const(ElementType[]) cdata()   const { return matrix_.cdata; }
		size_t               rows()    const { return rows_; }
		size_t               columns() const { return cols_; }
		bool                 empty()   const { return major_ != 0; }
		
		auto front() {
			static if( isRowMajor ) return row(0);
			else                    return column(0);
		}
		
		auto back() {
			static if( isRowMajor ) return row( rows_ - 1);
			else                    return column( cols_ - 1);
		}
	}
	
private:
	mixin SliceIndex2dMessages;

	void assignMatrix_( ref MatrixRef rhs ) {
		static if( isView )
			matrix_ = rhs;
		else
			matrix_ = MatrixRef( rhs.ptr );
	}

	MatrixRef matrix_;
	size_t    rowStart_, colStart_;
	size_t    rows_, cols_;
}
