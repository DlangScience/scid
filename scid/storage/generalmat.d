module scid.storage.generalmat;

import scid.internal.assertmessages;
import scid.storage.cowmatrix;
import scid.storage.generalmatview;
import scid.common.storagetraits;
import scid.matrix;
import scid.vector;
import scid.ops.eval, scid.ops.common;
import std.algorithm;
import scid.storage.external;


template GeneralMatrixStorage( ElementOrMatrix, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isScalar!(BaseElementType!ElementOrMatrix) ) {
	
	static if( isScalar!ElementOrMatrix )
		alias BasicGeneralMatrixStorage!( CowMatrixRef!(ElementOrMatrix,order_) ) GeneralMatrixStorage;
	else
		alias BasicGeneralMatrixStorage!( ElementOrMatrix ) GeneralMatrixStorage;
}

struct BasicGeneralMatrixStorage( ContainerRef_ ) {
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
	
	alias GeneralMatrixViewStorage!ContainerRef_                View;
	alias containerRef_                                         this;
	alias BasicGeneralMatrixStorage!(TransposedOf!ContainerRef) Transposed;
	alias BasicGeneralMatrixViewStorage!( ExternalMatrix!(ElementType, storageOrder, ContainerRef) )
		Temporary;
	
	enum storageType = MatrixStorageType.General;
	
	this( A ... )( auto ref A containerArgs ) if( A.length > 0 && !is( A[0] : ContainerRef ) ) {
		containerRef_ = ContainerRef( containerArgs );	
	}
	
	this( this ) {
		containerRef_ = ContainerRef( containerRef_.ptr );
	}
	
	void resize( size_t rows, size_t columns, void* ) {
		if( containerRef_.RefCounted.isInitialized() ) 
			containerRef_.resize( rows, columns, null );
		else
			containerRef_ = ContainerRef( rows, columns, null );
	}
	
	void resize( size_t rows, size_t columns ) {
		if( containerRef_.RefCounted.isInitialized() ) 
			containerRef_.resize( rows, columns );
		else
			containerRef_ = ContainerRef( rows, columns );
	}
	
	void copy( Transpose tr = Transpose.no, Source )( auto ref Source source ) if( isGeneralMatrixStorage!Source ) {
		enum srcOrder = transposeStorageOrder!( Source.storageOrder, tr ) ;
		static if( (!tr || !isComplex!ElementType) && srcOrder == storageOrder && is( Source : BasicGeneralMatrixViewStorage!ContainerRef ) ) {
			containerRef_ = ContainerRef( source.matrix.ptr, source.firstIndex, source.rows, source.columns );
		} else static if( (!tr || !isComplex!ElementType) && srcOrder == storageOrder && is( Source : typeof( this ) ) ) {
			containerRef_ = ContainerRef( source.containerRef_.ptr );
		} else {
			resize( source.rows, source.columns, null );
			generalMatrixCopy!tr( source, this );
		}
	}
	
	void forceRefAssign( ref typeof( this ) rhs ) {
		containerRef_ = rhs.matrix;
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		move( rhs.containerRef_, containerRef_ );
		return this;
	}
	
	typeof( this ) slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( ContainerRef( containerRef_.ptr, rowStart, rowEnd - rowStart, colStart, colEnd - colStart ) );
	}
	
	RowView row( size_t i )
	in {
		assert( i < rows, sliceMsg_(i,0,i,columns) );
	} body {
		static if( isRowMajor )
			return typeof( return )( containerRef_, i * leading, columns );
		else
			return typeof( return )( containerRef_, i, columns, leading );
	}
	
	ColumnView column( size_t j )
	in {
		assert( j < columns, sliceMsg_(0,j,rows,j) );
	} body {
		static if( isRowMajor )
			return typeof( return )( containerRef_, j, rows, leading );
		else
			return typeof( return )( containerRef_, j * leading, rows );
	}
	
	RowView rowSlice( size_t i, size_t start, size_t end )
	in {
		assert( i < rows && start < end && end <= columns, sliceMsg_(i,start,i,end) );
	} body {
		static if( isRowMajor )
			return typeof( return )( containerRef_, i * leading + start, end-start );
		else
			return typeof( return )( containerRef_, i + start * leading, end-start, leading );
	}
	
	ColumnView columnSlice( size_t j, size_t start, size_t end )
	in {
		assert( j < columns && start < end && end <= rows, sliceMsg_(start,j,end,j) );
	} body {
		static if( isRowMajor )
			return typeof( return )( containerRef_, j + start * leading, end-start, leading );
		else
			return typeof( return )( containerRef_, j * leading + start, end-start );
	}
	
	View view( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( containerRef_, rowStart, rowEnd - rowStart, colStart, colEnd - colStart );	
	}
	
	void popFront() { containerRef_.popFront(); }
	void popBack()  { containerRef_.popBack(); }
	
	@property {
		ContainerRef matrix()     { return containerRef_; }
		bool      empty()   const { return isInitd_() ? containerRef_.empty  : true; }
		size_t    length()  const { return isInitd_() ? containerRef_.major  : 0;    }
		size_t    rows()    const { return isInitd_() ? containerRef_.rows   : 0;    }
		size_t    columns() const { return isInitd_() ? containerRef_.columns: 0;    }
		size_t    major()   const { return isInitd_() ? containerRef_.major  : 0;    }
		size_t    minor()   const { return isInitd_() ? containerRef_.minor  : 0;    }
		size_t    leading() const { return isInitd_() ? containerRef_.leading: 0;    }
		
		MajorView front() {
			return typeof(return)( containerRef_, 0, minor );
		}
		
		MajorView back() {
			return typeof(return)( containerRef_, (containerRef_.major - 1) * containerRef_.minor, minor );
		}
	}
	
	ContainerRef containerRef_;
	
	/** Promotions for this type */
	private import scid.storage.array;
	template Promote( Other ) {
		static if( isMatrixStorage!Other ) {
			alias BasicGeneralMatrixStorage!( Promotion!((typeof(this)).ContainerRef,Other.ContainerRef) ) Promote;
		} else static if( isScalar!Other ) {
			alias BasicGeneralMatrixStorage!( Promotion!((typeof(this)).ContainerRef, Other ) ) Promote;
		}
	}
	
	mixin GeneralMatrixScalingAndAddition;
private:
	mixin MatrixErrorMessages;

	this()( ContainerRef containerRef ) {
		containerRef_ = containerRef;
	}
	
	bool isInitd_() const {
		return containerRef_.RefCounted.isInitialized();
	}
}
