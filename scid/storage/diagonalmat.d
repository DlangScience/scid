module scid.storage.diagonalmat;

import scid.matvec;
import scid.common.traits;
import scid.common.meta;
import scid.ops.common, scid.ops.expression, scid.ops.eval;
import scid.storage.cowarray;

import scid.internal.assertmessages;

import scid.blas;
import std.algorithm;

template DiagonalMatrixStorage( ElementOrArray )
		if( isScalar!(BaseElementType!ElementOrArray) ) {
	
	static if( isScalar!ElementOrArray )
		alias BasicDiagonalMatrixStorage!( CowArrayRef!(ElementOrArray), DiagonalMatrixStorageType.Root ) DiagonalMatrixStorage;
	else
		alias BasicDiagonalMatrixStorage!( ElementOrArray, DiagonalMatrixStorageType.Root ) DiagonalMatrixStorage;
}

enum DiagonalMatrixStorageType {
	Root,
	View,
	Slice
}

struct BasicDiagonalMatrixStorage( ContainerRef_, DiagonalMatrixStorageType type_ ) {
	alias BaseElementType!ContainerRef ElementType;
	alias ContainerRef_ ContainerRef;
	
	alias Vector!(DiagonalMatrixSubVectorStorage!( ContainerRef, VectorType.Column )) ColumnView;
	alias Vector!(DiagonalMatrixSubVectorStorage!( ContainerRef, VectorType.Row )) RowView;
	
	alias BasicDiagonalMatrixStorage!( ContainerRef_, DiagonalMatrixStorageType.View )
		View;
	
	private enum isRoot  = (type_ == DiagonalMatrixStorageType.Root);
	private enum isView  = (type_ == DiagonalMatrixStorageType.View);
	private enum isSlice = (type_ == DiagonalMatrixStorageType.Slice);
	
	static if( !isView )
		alias BasicDiagonalMatrixStorage!( ContainerRef_, DiagonalMatrixStorageType.Slice )
			Slice;
	else
		alias View Slice;
	
	this( A ... )( A containerArgs ) if( !is( A[ 0 ] : ElementType[][] ) && !is( A[ 0 ] : ContainerRef ) ) {
		containerRef_ = ContainerRef( containerArgs );
	}
	
	static if( !isRoot ) {
		this()( ref ContainerRef containerRef, size_t rowStart, size_t rows, size_t colStart, size_t cols ) {
			containerRef_ = containerRef;
			rows_         = rows;
			cols_         = cols;
			rowStart_     = rowStart;
			colStart_     = colStart;
		}
	}
	
	this()( ElementType[][] mat ) {
		containerRef_ = ContainerRef( mat.length );
		foreach( i ; 0 .. mat.length )
			containerRef_.indexAssign( mat[i][i], i );
	}
	
	this( this ) {
		if( isInitd_() )
			containerRef_ = ContainerRef( containerRef_.ptr );
	}
	
	void forceRefAssign( ref typeof(this) rhs ) {
		containerRef_ = rhs.containerRef_;
		static if( !isRoot ) {
			rows_     = rhs.rows_;
			cols_     = rhs.cols_;
			rowStart_ = rhs.rowStart_;
			colStart_ = rhs.colStart_;
		}
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		forceRefAssign( rhs );
		return this;
	}
	
	ElementType index( size_t i, size_t j ) const
	in {
		assert( i < rows && j < columns, boundsMsg_( i, j ) );
	} body {
		return i == j ? containerRef_.index( i ) : Zero!ElementType;
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i, size_t j )
	in {
		assert( i < rows && j < columns, boundsMsg_( i, j ) );
		assert( i == j, msgPrefix_ ~ "Assignment to zero element in diagonal matrix." );
	} body {
		containerRef_.indexAssign!op( rhs, i );
	}
	
	Slice slice( size_t rowStart, size_t colStart, size_t rowEnd, size_t colEnd ) {
		return typeof( return )( containerRef_, rowStart, rowEnd - rowStart, colStart, colEnd - colStart );
	}
	
	View view( size_t rowStart, size_t colStart, size_t rowEnd, size_t colEnd ) {
		return typeof( return )( containerRef_, rowStart, rowEnd - rowStart, colStart, colEnd - colStart );
	}
	
	RowView row( size_t i )
	in {
		assert( i < rows, sliceMsg_( i,	0, i, columns ) );
	} body {
		if( isRoot || i >= colStart_ )
			return typeof( return )( containerRef_, i + rowStart_ , i - colStart_, columns );
		else
			return typeof( return )( columns );
	}
	
	ColumnView column( size_t j )
	in {
		assert( j < columns, sliceMsg_( 0, j, rows, j ) );
	} body {
		if( isRoot || j >= rowStart_ )
			return typeof( return )( containerRef_, j + colStart_, j - rowStart_, rows );
		else
			return typeof( return )( rows );
	}
	
	void resize( size_t rows, size_t columns, void * ) {
		assert( rows == columns );
		containerRef_.resize( rows, null );
	}
	
	void resize( size_t rows, size_t columns ) {
		assert( rows == columns );
		containerRef_.resize( rows );
	}
	
	void copy( Transpose tr = Transpose.no, Source )( auto ref Source source  ) if( is( Source : typeof(this))) {
		this = source;
	}
	
	void copy( Transpose tr = Transpose.no, Source )( auto ref Source src ) if( !is( Source : typeof( this ) ) ) {
		resize( src.rows, src.columns, null );
		auto start = max( rowStart_, colStart_ );
		auto end   = min( rows_, cols_ );
		foreach( i ; start .. end ) {
			auto e = src.index( i, i );
			static if( isComplexScalar!ElementType && tr )
				e = gconj( e );
			containerRef_.indexAssign( e, i );
		}
	}
	
	void copyRight( Transpose tr = Transpose.no, Dest )( auto ref Dest dest ) {
		assert( dest.rows == rows_ && dest.columns == columns_,
		        dimMismatch_( dest.rows, dest.columns, "assignment" ) );
		auto start = max( rowStart_, colStart_ );
		auto end   = min( rows_, cols_ );
		foreach( i ; start .. end ) {
			auto e = containerRef_.index( i );
			static if( isComplexScalar!ElementType && tr )
				e = gconj( e );
			dest.indexAssign( e, i,i );
		}
	}
	
	void scale( ElementType alpha ) {
		if( isInitd_() )
			blas.scal( containerRef_.length, alpha, containerRef_.data , 1 );
	}	
	
	void scaledAddition( Transpose tr = Transpose.no, Source )( ElementType alpha, auto ref Source source ) {
		auto start = max( rowStart_, colStart_ );
		auto end   = min( rows_, cols_ );
		foreach( i ; start .. end ) {
			auto e = source.index( i, i );
			static if( isComplexScalar!ElementType && tr )
				e = gconj( e );
			containerRef_.indexAssign!"+"( e * alpha, i );
		}
	}
	
	void scaledAdditionRight( Transpose tr = Transpose.no, Source )( ElementType alpha, auto ref Source source ) {
		assert( dest.rows == rows_ && dest.columns == columns_,
		        dimMismatch_( dest.rows, dest.columns, "assignment" ) );
		auto start = max( rowStart_, colStart_ );
		auto end   = min( rows_, cols_ );
		foreach( i ; start .. end ) {
			auto e = containerRef_.index( i );
			static if( isComplexScalar!ElementType && tr )
				e = gconj( e );
			dest.indexAssign!"+"( e * alpha, i,i );
		}
	}
	
	void matrixProduct( Transpose transA = Transpose.no, Transpose transB = Transpose.no, A,B )
			( ElementType alpha, auto ref A a, auto ref B b, ElementType beta ) {
		auto start = max( rowStart_, colStart_ );
		auto end   = min( rows_, cols_ );
		auto d     = containerRef_.data[ start .. end ];
		
		auto m = transA ? a.columns : a.rows;
		auto n = transB ? b.rows : b.columns;
		
		if( !alpha ) {
			if( !beta ) resize( m, n );
			else {
				assert( rows == m && columns == n );
				blas.scal( d.length, Zero!ElementType, d.ptr, 1 );
			}
			
			return;
		}
		
		//static if( is( A : typeof( this ) ) && is( B : typeof( this ) ) ) {
		//    blas.sbmv( 'u', d.length, 0, alpha, a.cdata, 1, b.cdata, 1, beta, d.ptr, 1 );
		//} else {
			if( beta ) {	
				foreach( i, ref destElem ; d ) {
					destElem *= beta;
					destElem += alpha * rowColumnDot!( transA, transB )( a, i, b, i );
				}
			} else {
				foreach( i, ref destElem ; d ) {
					destElem = alpha * rowColumnDot!( transA, transB )( a, i, b, i );
				}
			}
		//}
	}
	
	@property {
		const(ElementType)* cdata() const { return containerRef_.cdata; }
		ElementType*        data()        { return containerRef_.data; }
		static if( isRoot ) {
			size_t size() const { return isInitd_() ? containerRef_.length() : 0; }
			alias size rows;
			alias size columns;
		} else {
			size_t rows() const { return rows_; }
			size_t columns() const { return cols_; }
		}
	}
	
private:
	mixin MatrixErrorMessages;
	
	static if( !isRoot ) {
		size_t rowStart_, rows_;
		size_t colStart_, cols_;
	} else {
		enum rowStart_ = 0;
		enum colStart_ = 0;
		alias size rows_;
		alias size cols_;
	}
	
	bool isInitd_() const {
		return containerRef_.RefCounted.isInitialized();
	}

	ContainerRef containerRef_;
}

struct DiagonalMatrixSubVectorStorage( ContainerRef_, VectorType vtype_ ) {
	alias ContainerRef_ ContainerRef;
	alias BaseElementType!ContainerRef ElementType;
	alias vtype_ vectorType;
	
	alias DiagonalMatrixSubVectorStorage!( ContainerRef, transposeVectorType!vtype_ )
		Transposed;
	
	this( ref ContainerRef containerRef, size_t realIndex, size_t fakeIndex, size_t length ) {
		containerRef_ = containerRef;
		realIndex_    = realIndex;
		fakeIndex_    = fakeIndex;
		length_       = length;
	}
	
	this( size_t length ) {
		length_ = length;
		fakeIndex_ = length;
	}
	
	void forceRefAssign( ref typeof(this) rhs ) {
		containerRef_ = rhs.containerRef_;
		realIndex_    = rhs.realIndex_;
		fakeIndex_    = rhs.fakeIndex_;
		length_       = rhs.length_;
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		forceRefAssign( rhs );
		return this;
	}
	
	ElementType index( size_t i ) const
	in {
		assert( i < length_, boundsMsg_( i ) );
	} body {
		if( i == fakeIndex_ )
			return containerRef_.index( realIndex_ );
		else
			return 0;
	}
	
	void indexAssign( string op = "" )( ElementType rhs, size_t i )
	in {
		assert( i < length_, boundsMsg_( i ) );
		assert( i == fakeIndex_, msgPrefix_ ~ "Assignment to zero element in diagonal matrix." );
	} body {
		containerRef_.indexAssign!op( rhs, realIndex_ );	
	}
	
	typeof( this ) slice( size_t start, size_t end ) {
		typeof( this ) r;
		r.containerRef_ = containerRef_;
		r.length_ = end - start;
		if( fakeIndex_ < start )
			r.fakeIndex_ = r.length_;
		else
			r.fakeIndex_ = fakeIndex_ - start;
		r.realIndex_ = realIndex_;
		
		return r;
	}
	
	alias slice view;
	
	void popFront()
	in {
		assert( !empty, emptyMsg_("popFront()") );
	} body {
		if( fakeIndex_ > 0 )
			-- fakeIndex_;
		else
			fakeIndex_ = length_;
		-- length_;
	}
	
	void popBack()
	in {
		assert( !empty, emptyMsg_("popBack()") );
	} body {
		-- length_;
	}
	
	void copy( Transpose tr = Transpose.no, Source )( auto ref Source source ) {
		assert( source.length == length_, lengthMismatch_( source.length, "assignment" ) );
		nonZeroAssign_!tr( source.index( fakeIndex_ ), realIndex_ );
	}
	
	void scaledAddition( Transpose tr = Transpose.no, Source )( ElementType alpha, auto ref Source source ) {
		assert( source.length == length_, lengthMismatch_( source.length, "scaled addition" ) );
		nonZeroScaledAddition_!tr( alpha, source.index( fakeIndex_ ) );
	}
	
	ElementType dot( Transpose tr = Transpose.no, Right )( auto ref Right right ) {
		assert( right.length == length_, lengthMismatch_( right.length, "dot product" ) );
		if( hasNonZero_() )
			return right.index( fakeIndex_ ) * nonZero_!tr();
		else
			return Zero!ElementType;
	}
	
	void copyRight( Transpose tr = Transpose.no, Dest )( auto ref Dest dest ) {
		assert( dest.length == length_, lengthMismatch_( dest.length, "assignment" ) );
		evalScaling( Zero!ElementType, dest );
		if( hasNonZero_() )
			dest.indexAssign( nonZero_!tr(), fakeIndex_ );
	}
	
	void scaledAdditionRight( Transpose tr = Transpose.no, Dest )( ElementType alpha, auto ref Dest dest ) {
		assert( dest.length == length_, lengthMismatch_( dest.length, "scaled addition" ) );
		if( hasNonZero_() )
			dest.indexAssign!"+"( alpha * nonZero_!tr(), fakeIndex_ );
	}
	
	@property {
		bool empty() const {
			return length_ == 0;
		}
		
		size_t length() const {
			return length_;
		}
	}
	
private:
	mixin ArrayErrorMessages;

	bool hasNonZero_() const {
		return fakeIndex_ < length_;
	}
		
	ElementType nonZero_( Transpose tr = Transpose.no )() const {
		auto r = containerRef_.index( realIndex_ );
		static if( isComplexScalar!ElementType && tr )
			return gconj( r );
		else
			return r;
	}
	
	void nonZeroScaledAddition_( Transpose tr = Transpose.no )( ElementType alpha, ElementType what ) {
		if( !hasNonZero_() )
			return;
		
		static if( isComplexScalar!ElementType && tr )
			what = gconj( what );
		containerRef_.indexAssign!"+"( alpha * what, realIndex_ );	
	}
	
	void nonZeroAssign_( Transpose tr = Transpose.no )( ElementType alpha, ElementType what ) {
		if( !hasNonZero_() )
			return;
		
		static if( isComplexScalar!ElementType && tr )
			what = gconj( what );
		containerRef_.indexAssign( what, realIndex_ );	
	}
		
	size_t length_, fakeIndex_, realIndex_;
	ContainerRef containerRef_;
}