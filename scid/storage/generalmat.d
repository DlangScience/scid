module scid.storage.generalmat;

import scid.internal.assertmessages;
import scid.storage.cowmatrix;
import scid.storage.generalmatview;
import scid.common.traits;
import scid.matrix;
import scid.vector;
import scid.internal.hlblas;
import std.algorithm;


template GeneralMatrixStorage( ElementOrMatrix, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isFortranType!(BaseElementType!ElementOrMatrix) ) {
	
	static if( isFortranType!ElementOrMatrix )
		alias BasicGeneralMatrixStorage!( CowMatrixRef!(ElementOrMatrix,order_) ) GeneralMatrixStorage;
	else
		alias BasicGeneralMatrixStorage!( ElementOrMatrix ) GeneralMatrixStorage;
}

struct BasicGeneralMatrixStorage( MatrixRef_ ) {
	mixin GeneralMatrixTypes!MatrixRef_;
	
	alias GeneralMatrixViewStorage!MatrixRef_ View;
	alias matrix_                             this;
	
	enum storageType = MatrixStorageType.General;
	
	this( size_t newRows, size_t newCols ) {
		matrix_ = MatrixRef( newRows, newCols );
	}	
	
	this( size_t newRows, size_t newCols, void* ) {
		matrix_ = MatrixRef( newRows, newCols, null );
	}	
	
	this( size_t newRows, ElementType[] initializer ) {
		matrix_ = MatrixRef( newRows, initializer );
	}
	
	this( this ) {
		matrix_ = MatrixRef( matrix_.ptr );
	}
	
	void resizeOrClear( size_t rows, size_t columns, void* ) {
		if( matrix_.RefCounted.isInitialized() ) 
			matrix_.resizeOrClear( rows, columns, null );
		else
			matrix_ = MatrixRef( rows, columns, null );
	}
	
	void resizeOrClear( size_t rows, size_t columns ) {
		if( matrix_.RefCounted.isInitialized() ) 
			matrix_.resizeOrClear( rows, columns );
		else
			matrix_ = MatrixRef( rows, columns );
	}
	
	void copy( Transpose tr, Source )( auto ref Source source ) {
		enum srcOrder = transposeStorageOrder!( Source.storageOrder, tr ) ;
		static if( srcOrder == storageOrder && is( Source : BasicGeneralMatrixViewStorage!MatrixRef ) ) {
			matrix_ = MatrixRef( source.matrix.ptr, source.firstIndex, source.rows, source.columns );
		} else static if( srcOrder == storageOrder && is( Source : typeof( this ) ) ) {
			matrix_ = MatrixRef( source.matrix_.ptr );
		} else {
			resizeOrClear( source.rows, source.columns, null );
			hlGeCopy!( srcOrder, storageOrder )(
				rows, columns,
				source.cdata, source.leading,
				this.data,
				leading
			);
		}
	}
	
	void forceRefAssign( ref typeof( this ) rhs ) {
		matrix_ = rhs.matrix;
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		move( rhs.matrix_, matrix_ );
		return this;
	}
	
	typeof( this ) slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( MatrixRef( matrix_.ptr, rowStart, rowEnd - rowStart, colStart, colEnd - colStart ) );
	}
	
	RowView row( size_t i )
	in {
		assert( i < rows, sliceMsg_(i,0,i,columns) );
	} body {
		static if( isRowMajor )
			return typeof( return )( matrix_, i * leading, columns );
		else
			return typeof( return )( matrix_, i, columns, leading );
	}
	
	ColumnView column( size_t j )
	in {
		assert( j < columns, sliceMsg_(0,j,rows,j) );
	} body {
		static if( isRowMajor )
			return typeof( return )( matrix_, j, rows, leading );
		else
			return typeof( return )( matrix_, j * leading, rows );
	}
	
	RowView rowSlice( size_t i, size_t start, size_t end )
	in {
		assert( i < rows && start < end && end <= columns, sliceMsg_(i,start,i,end) );
	} body {
		static if( isRowMajor )
			return typeof( return )( matrix_, i * leading + start, end-start );
		else
			return typeof( return )( matrix_, i + start * leading, end-start, leading );
	}
	
	ColumnView columnSlice( size_t j, size_t start, size_t end )
	in {
		assert( j < columns && start < end && end <= rows, sliceMsg_(start,j,end,j) );
	} body {
		static if( isRowMajor )
			return typeof( return )( matrix_, j + start * leading, end-start, leading );
		else
			return typeof( return )( matrix_, j * leading + start, end-start );
	}
	
	View view( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( matrix_, rowStart, rowEnd - rowStart, colStart, colEnd - colStart );	
	}
	
	void popFront() { matrix_.popFront(); }
	void popBack()  { matrix_.popBack(); }
	
	@property {
		MatrixRef matrix()        { return matrix_; }
		bool      empty()   const { return matrix_.empty; }
		size_t    length()  const { return matrix_.major; }
		
		size_t    rows()    const { return matrix_.rows; }
		size_t    columns() const { return matrix_.columns; }
		size_t    major()   const { return matrix_.major; }
		size_t    minor()   const { return matrix_.minor; }
		size_t    leading() const { return matrix_.leading; }
		
		MajorView front() {
			return typeof(return)( matrix_, 0, minor );
		}
		
		MajorView back() {
			return typeof(return)( matrix_, (matrix_.major - 1) * matrix_.minor, minor );
		}
	}
	
	MatrixRef matrix_;
	
	mixin GeneralMatrixAxpyScal;
private:
	mixin SliceIndex2dMessages;

	this( MatrixRef matrixRef ) {
		matrix_ = matrixRef;
	}
}
