module scid.storage.generalmat;

import scid.internal.assertmessages;
import scid.storage.cowmatrix;
import scid.storage.generalmatview;
import scid.common.storagetraits;
import scid.matrix;
import scid.vector;
import scid.ops.eval, scid.ops.common;
import std.algorithm;


template GeneralMatrixStorage( ElementOrMatrix, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isFortranType!(BaseElementType!ElementOrMatrix) ) {
	
	static if( isFortranType!ElementOrMatrix )
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
	
	enum storageType = MatrixStorageType.General;
	
	this( size_t newRows, size_t newCols ) {
		containerRef_ = ContainerRef( newRows, newCols );
	}	
	
	this( size_t newRows, size_t newCols, void* ) {
		containerRef_ = ContainerRef( newRows, newCols, null );
	}	
	
	this( size_t newRows, ElementType[] initializer ) {
		containerRef_ = ContainerRef( newRows, initializer );
	}
	
	this( ElementType[][] initializer ) {
		if( initializer.length == 0 || initializer[0].length == 0)
			return;
		
		containerRef_ = ContainerRef( initializer.length, initializer[0].length, null );
		foreach( i ; 0 .. this.rows ) {
			assert( initializer[i].length == this.columns );
			this.row( i )[] = initializer[ i ].t;
		}
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
	
	void copy( Transpose tr, Source )( auto ref Source source ) {
		enum srcOrder = transposeStorageOrder!( Source.storageOrder, tr ) ;
		static if( srcOrder == storageOrder && is( Source : BasicGeneralMatrixViewStorage!ContainerRef ) ) {
			containerRef_ = ContainerRef( source.matrix.ptr, source.firstIndex, source.rows, source.columns );
		} else static if( srcOrder == storageOrder && is( Source : typeof( this ) ) ) {
			containerRef_ = ContainerRef( source.containerRef_.ptr );
		} else {
			resize( source.rows, source.columns, null );
			generalMatrixCopy!( srcOrder, storageOrder )(
				rows, columns,
				source.cdata, source.leading,
				this.data,
				leading
			);
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
		ContainerRef matrix()        { return containerRef_; }
		bool      empty()   const { return containerRef_.empty; }
		size_t    length()  const { return containerRef_.major; }
		
		size_t    rows()    const { return containerRef_.rows; }
		size_t    columns() const { return containerRef_.columns; }
		size_t    major()   const { return containerRef_.major; }
		size_t    minor()   const { return containerRef_.minor; }
		size_t    leading() const { return containerRef_.leading; }
		
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
		} else static if( isFortranType!Other ) {
			alias BasicGeneralMatrixStorage!( Promotion!((typeof(this)).ContainerRef, Other ) ) Promote;
		}
	}
	
	mixin GeneralMatrixScalingAndAddition;
private:
	mixin MatrixErrorMessages;

	this( ContainerRef containerRef ) {
		containerRef_ = containerRef;
	}
}
