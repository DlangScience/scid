module scid.storage.packedsubmat;

import scid.common.storagetraits;
import scid.matrix, scid.vector;
import scid.internal.assertmessages;
import scid.storage.packedsubvec;
import scid.ops.common;
import std.algorithm;

enum SubMatrixType {
	View,
	Slice
}

struct PackedSubMatrixStorage( ContainerRef_, SubMatrixType type_ ) {
	alias ContainerRef_                                                         ContainerRef;
	alias BaseElementType!ContainerRef                                          ElementType;
	alias storageOrderOf!ContainerRef                                           storageOrder;
	alias typeof(this)                                                          Slice;
	alias typeof(this)                                                          View;
	alias Vector!( PackedSubVectorStorage!( ContainerRef, VectorType.Row ) )    RowView;
	alias Vector!( PackedSubVectorStorage!( ContainerRef, VectorType.Column ) ) ColumnView;
	alias ColumnView                                                            DiagonalView;
	
	alias PackedSubMatrixStorage!( TransposedOf!ContainerRef, type_ ) Transposed;
	
	enum isView     = type_ == SubMatrixType.View;
	enum isRowMajor = (storageOrder == StorageOrder.RowMajor );
	enum storageType = MatrixStorageType.Virtual;
	
	this( ref ContainerRef containerRef, size_t rowStart, size_t numRows, size_t colStart, size_t numCols ) {
		assignMatrix_( containerRef );
		rowStart_ = rowStart; rows_ = numRows;
		colStart_ = colStart; cols_ = numCols;
	}
	
	static if( !isView ) {
		this( this ) {
			assignMatrix_( containerRef_ );
		}
	}
		
	void forceRefAssign( ref typeof(this) rhs ) {
		containerRef_ = rhs.containerRef_;
		rowStart_ = rhs.rowStart_; rows_ = rhs.rows_;
		colStart_ = rhs.colStart_; cols_ = rhs.cols_;
	}
	
	ref typeof( this ) opAssign( typeof(this) rhs ) {
		swap( rhs.containerRef_, containerRef_ );
		rowStart_ = rhs.rowStart_; rows_ = rhs.rows_;
		colStart_ = rhs.colStart_; cols_ = rhs.cols_;
		return this;
	}
	
	ElementType index( size_t i, size_t j ) const 
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		return containerRef_.index( i + rowStart_, j + colStart_ );
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		containerRef_.indexAssign!op( rhs, i + rowStart_, j + colStart_ );
	}
	
	RowView row( size_t i )
	in {
		assert( i < rows_, sliceMsg_(i,0,i,cols_) );
	} body {
		return typeof( return )( containerRef_, i + rowStart_, colStart_, cols_ );
	}
	
	ColumnView column( size_t j )
	in {
		assert( j < cols_, sliceMsg_(0,j,rows_,j) );
	} body {
		return typeof( return )( containerRef_, j + colStart_, rowStart_, rows_ );
	}
	
	RowView rowSlice( size_t i, size_t start, size_t end )
	in {
		assert( i < rows_ && start < end && end <= cols_, sliceMsg_(i,start,i,end) );
	} body {
		return typeof( return )( containerRef_, i + rowStart_, start + colStart_, end - start );
	}
	
	ColumnView columnSlice( size_t j, size_t start, size_t end )
	in {
		assert( j < cols_ && start < end && end <= rows_, sliceMsg_(start,j,end,j) );
	} body {
		return typeof( return )( containerRef_, j + colStart_, start + rowStart_, end - start );
	}
	
	View view( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( containerRef_, rowStart + rowStart_, rowEnd - rowStart, colStart + colStart_, colEnd - colStart );
	}
	
	Slice slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( containerRef_, rowStart + rowStart_, rowEnd - rowStart, colStart + colStart_, colEnd - colStart );
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
		ContainerRef           matrix()        { return containerRef_; }
		size_t              rows()    const { return rows_; }
		size_t              columns() const { return cols_; }
		bool                empty()   const { return major_ != 0; }
		
		auto front() {
			static if( isRowMajor ) return row(0);
			else                    return column(0);
		}
		
		auto back() {
			static if( isRowMajor ) return row( rows_ - 1);
			else                    return column( cols_ - 1);
		}
	}
	
	/** Promotions for this type are inherited either from its container or from general matrix. */
	template Promote( Other ) {
		private import scid.storage.generalmat;
		alias Promotion!( GeneralMatrixStorage!ElementType, Other ) Promote;
	}
	
private:
	mixin MatrixErrorMessages;

	void assignMatrix_( ref ContainerRef rhs ) {
		static if( isView )
			containerRef_ = rhs;
		else
			containerRef_ = ContainerRef( rhs.ptr );
	}

	ContainerRef containerRef_;
	size_t    rowStart_, colStart_;
	size_t    rows_, cols_;
}
