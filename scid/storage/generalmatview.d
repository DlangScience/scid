module scid.storage.generalmatview;

import scid.internal.assertmessages;
import scid.storage.cowmatrix;
import scid.storage.generalmat;
import scid.common.storagetraits;
import scid.common.meta;
import scid.matrix, scid.vector;
import std.algorithm;
import scid.ops.eval, scid.ops.common;
import scid.blas;

template GeneralMatrixViewStorage( ElementOrMatrix, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isScalar!(BaseElementType!ElementOrMatrix) ) {
	
	static if( isScalar!ElementOrMatrix )
		alias BasicGeneralMatrixViewStorage!( CowMatrixRef!(ElementOrMatrix,order_) ) GeneralMatrixViewStorage;
	else
		alias BasicGeneralMatrixViewStorage!( ElementOrMatrix ) GeneralMatrixViewStorage;
}

struct BasicGeneralMatrixViewStorage( ContainerRef_ ) {
	alias ContainerRef_                  ContainerRef;
	alias BaseElementType!ContainerRef   ElementType;
	alias typeof(this)                   Slice;
	alias storageOrderOf!ContainerRef    storageOrder;
	
	enum bool isRowMajor = ( storageOrder == StorageOrder.RowMajor );
	
	static if( isRowMajor ) {
		alias VectorView!(ContainerRef, VectorType.Row)           RowView;
		alias StridedVectorView!(ContainerRef, VectorType.Column) ColumnView;
		alias StridedVectorView!(ContainerRef, VectorType.Row)    DiagonalView;
		alias RowView                                             MajorView;
		alias ColumnView                                          MinorView;
	} else {
		alias VectorView!(ContainerRef, VectorType.Column)        ColumnView;
		alias StridedVectorView!(ContainerRef, VectorType.Row)    RowView;
		alias StridedVectorView!(ContainerRef, VectorType.Column) DiagonalView;
		alias ColumnView                                          MajorView;
		alias RowView                                             MinorView;
	}
	
	alias typeof( this )                                            View;
	
	alias BasicGeneralMatrixStorage!( TransposedOf!(MatrixTypeOf!ContainerRef) )
		Transposed;
	
	this()( ref ContainerRef containerRef, size_t rowStart, size_t numRows, size_t colStart, size_t numCols, size_t offset=0 ) {
		containerRef_ = containerRef;
		firstIndex_   = containerRef_.mapIndex( rowStart, colStart ) + offset;
		rows_         = numRows;
		cols_         = numCols;
		leading_      = containerRef_.leading;
	}
	
	this( A ... )( A args ) if( A.length > 0 && !is( A[ 0 ] : ContainerRef ) ) {
		containerRef_ = ContainerRef( args );
		rows_         = containerRef_.rows;
		cols_         = containerRef_.columns;
		leading_      = minor_;
		firstIndex_   = 0;
	}
	
	void forceRefAssign( ref typeof(this) rhs ) {
		this = rhs;
	}
	
	ref typeof( this ) opAssign( typeof( this ) rhs ) {
		move( rhs.containerRef_, containerRef_ );
		firstIndex_ = rhs.firstIndex_;
		rows_       = rhs.rows_;
		cols_       = rhs.cols_;
		leading_    = rhs.leading_;
		return this;
	}
	
	typeof( this ) slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		auto numRows = rowEnd - rowStart;
		auto numCols = colEnd - colStart;	
		return typeof( return )( containerRef_, rowStart, numRows, colStart, numCols, firstIndex_);
	}
	
	ElementType index( size_t i, size_t j ) const
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		return containerRef_.cdata[ map_( i, j ) ];
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < rows_ && j < cols_, boundsMsg_(i, j) );
	} body {
		mixin( "containerRef_.data[ map_( i, j ) ]" ~ op ~ "= rhs;" );
	}
	
	RowView row( size_t i )
	in {
		assert( i < rows, sliceMsg_(i,0,i,columns) );
	} body {
		static if( isRowMajor )
			return typeof( return )( containerRef_, firstIndex_ + i * leading, columns );
		else
			return typeof( return )( containerRef_, firstIndex_ + i, columns, leading );
	}
	
	ColumnView column( size_t j )
	in {
		assert( j < columns, sliceMsg_(0,j,rows,j) );
	} body {
		static if( isRowMajor )
			return typeof( return )( containerRef_, firstIndex_ + j, rows, leading );
		else
			return typeof( return )( containerRef_, firstIndex_ + j * leading, rows );
	}
	
	RowView rowSlice( size_t i, size_t start, size_t end )
	in {
		assert( i < rows_ && start < end && end <= cols_, sliceMsg_(i,start,i,end) );
	} body {
		static if( isRowMajor )
			return typeof( return )( containerRef_, i * leading + start + firstIndex_, end-start );
		else
			return typeof( return )( containerRef_, i + start * leading + firstIndex_, end-start, leading );
	}
	
	ColumnView columnSlice( size_t j, size_t start, size_t end )
	in {
		assert( j < cols_ && start < end && end <= rows_, sliceMsg_(start,j,end,j) );
	} body {
		static if( isRowMajor )
			return typeof( return )( containerRef_, j + start * leading + firstIndex_, end-start, leading );
		else
			return typeof( return )( containerRef_, j * leading + start + firstIndex_, end-start );
	}
	
	alias slice view;
	
	void resize( size_t rows, size_t columns, void* ) {
		assert( rows == rows_ && cols_ == columns, dimMismatch_( rows, columns ) );
	}
	
	void resize( size_t rows, size_t columns ) {
		resize( rows, columns, null );
		generalMatrixScaling!storageOrder( rows, columns, Zero!ElementType, this.data, leading );
	}
	
	void copy( Transpose tr = Transpose.no, Source )( auto ref Source source )
		if( is( Source M : BasicGeneralMatrixStorage!M ) || is( Source M : BasicGeneralMatrixViewStorage!M ) ) {
		resize( source.rows, source.columns, null );
		generalMatrixCopy!tr( source, this );
	}
	
	void popFront()
	in {
		assert( !empty, msgPrefix_ ~ "popFront on empty." );
	} body {
		-- major_;
		firstIndex_ += leading_;
	}
	
	void popBack()
	in {
		assert( !empty, msgPrefix_ ~ "popFront on empty." );
	} body {
		-- major_;
	}
	
	@property {
		ref ContainerRef     container()        { return containerRef_; }
		ElementType*         data()             { return containerRef_.data + firstIndex_ ; }
		const(ElementType)*  cdata()      const { return containerRef_.cdata + firstIndex_; }
		bool                 empty()      const { return major_ == 0; }
		size_t               length()     const { return major_; }
		size_t               rows()       const { return rows_; }
		size_t               columns()    const { return cols_; }
		size_t               major()      const { return major_; }
		size_t               minor()      const { return minor_; }
		size_t               leading()    const { return leading_; }
		size_t               firstIndex() const { return firstIndex_; }
		
		MajorView front() {
			
			return typeof( return )( containerRef_, firstIndex_, minor_ );
		}
		
		MajorView back() {
			return typeof( return )( containerRef_, firstIndex_ + (major_ - 1) * leading_, minor_ );
		}
	}
	
	/** Promotions for this type are inherited from GeneralMatrix */
	private import scid.storage.generalmat;
	template Promote( Other ) {
		alias Promotion!( BasicGeneralMatrixStorage!(MatrixTypeOf!ContainerRef), Other ) Promote;
	}
	
	mixin GeneralMatrixScalingAndAddition;
	
private:
	mixin MatrixErrorMessages;

	size_t map_( size_t i, size_t j ) const {
		static if( isRowMajor )
			return firstIndex_ + i * leading_ + j;
		else
			return firstIndex_ + j * leading_ + i;
	}
	
	static if( isRowMajor ) {
		alias rows_ major_;
		alias cols_ minor_;
	} else {
		alias rows_ minor_;
		alias cols_ major_;
	}
	
	size_t    rows_, cols_, firstIndex_, leading_;
	ContainerRef containerRef_;
}
