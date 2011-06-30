module scid.matrix;

import scid.vectorview;
import scid.vector;
import scid.internal.cowarray;

import std.array;
import std.conv;
import std.typecons;

import std.stdio;

enum StorageOrder { RowMajor, ColumnMajor }

template otherStorageOrder( StorageOrder order ) {
	static if( order == StorageOrder.RowMajor ) 
		enum otherStorageOrder = StorageOrder.ColumnMajor;
	else
		enum otherStorageOrder = StorageOrder.RowMajor;
}

/** General matrix implementation with value semantics. Storage order (row-major / column-major)
  * can be customized, along with the Array used for storage.
  *
  * To use with BLAS functions (and other functions taking a builtin array), use the data()
  * and cdata() methods. Whenever possible, use cdata() as data() might duplicate the memory.
  *
  * MatrixView-s are not supported yet. One can obtain views of individual rows/vectors, however,
  * through slicing eg. M[0..2][0], M[][2], M[2] etc. The views returned are VectorView-s so
  * they can be used directly with BLAS.
  */

Matrix!( T, order, arrayTemplate ) matrixNoInit( T, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray )( size_t numRows, size_t numColumns ) {
	return typeof( return ).clearNoInit( numRows, numColumns );
}

struct Matrix( T, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray ) {	
	
	/** The type of elements in the vector. */
	alias T                ElementType;
	
	/** The type of the Array that's wrapped. */
	alias ArrayTemplate!T    Wrapped;
	
	/** The RefCounted array type. */
	alias RefCounted!Wrapped Array;
	
	/** Is the matrix stored in row-major order? */
	enum isRowMajor = StorageOrder.RowMajor == storageOrder_;
	
	/** Create a new matrix of a given size (MxN). Initializes the data.  */
	this( size_t sizeRows, size_t sizeColumns, T initializeWith = T.init ) {
		setSliceParams_( sizeRows, sizeColumns );
		array_.RefCounted.initialize( sizeRows * sizeColumns, initializeWith );
	}
	
	/** Create a new matrix of a given size (MxN). Do not initialize the data.  */
	this( size_t sizeRows, size_t sizeColumns, void* ) {
		setSliceParams_( sizeRows, sizeColumns );
		array_.RefCounted.initialize( sizeRows * sizeColumns, null );
	}

	/** Construct a matrix from a built-in array of rows. */
	this( T[][] matrix ) {
		this = matrix;
	}
	
	/** Construct a matrix from a given array, size and leading dimension. */
	this( Array arr, size_t numRows, size_t numCols, size_t leadingDim = 0 ) {
		setSliceParams_( numRows, numCols, leadingDim );
		array_.RefCounted.initialize( *&arr.refCountedPayload() );
	}
	
	

	/** Postblit calls opAssign on the wrapped array. */
	this( this ) {
		if( array_.RefCounted.isInitialized() )
			array_.RefCounted.initialize( *&array_.refCountedPayload() );
	}
	
	/// Create a matrix slice of another matrix. Used by SlicedProxy and op*Slice*().
	private this( ref typeof(this) rhs, size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd  ) {
		rows_    = rowEnd - rowStart;
		cols_    = colEnd - colStart;
		leading_ = rhs.leading_;

		size_t startIndex = rhs.map_( rowStart, colStart );
		size_t endIndex   = rhs.map_( rowEnd - 1, colEnd - 1 ) + 1;
		array_.RefCounted.initialize( rhs.array_.refCountedPayload()[ startIndex .. endIndex ] );
	}
	
	/** Reset the matrix with new dimensions. */
	void reset( size_t newRows, size_t newCols, T initializeWith = T.init ) {
		reset( newRows, newCols, null );
		array_.refCountedPayload()[] = initializeWith;
	}
	
	/** Reset the matrix with new dimensions, do not initialize the data. */
	void reset( size_t newRows, size_t newCols, void* ) {
		setSliceParams_( newRows, newCols );
		
		if( !array_.RefCounted.isInitialized() )
				array_.RefCounted.initialize();
		
		if( newRows * newCols != array.length )
			array_.RefCounted.initialize( newRows * newCols, null );
	}
	
	/** Return the wrapped memory block. Writable. */
	@property T[] data() {
		assert(array_.RefCounted.isInitialized());
		return array_.data;
	}
	
	/** Return the wrapped memory block. Read-only. */
	@property const(T) [] cdata() const {
		assert(array_.RefCounted.isInitialized());
		return array_.cdata;
	}
	
	/** The number of rows of the matrix. */
	@property size_t rows()    const { return rows_; }
	
	/** The number of columns of the matrix. */
	@property size_t columns() const { return cols_; }
	
	/** The leading dimension of the matrix. */
	@property size_t leading() const { return leading_; }
	
	/** The ref-counted wrapped array. */
	@property Array  array()         { return array_; }
	
	/** Assign to another matrix */
	ref typeof(this) opAssign( typeof( this ) rhs ) {
		array_   = *&rhs.array_.refCountedPayload();
		rows_    = rhs.rows_;
		cols_    = rhs.cols_;
		leading_ = rhs.leading_;
		
		return this;
	}

	/** Assign to a built-in array of rows. */
	ref typeof(this) opAssign( T[][] matrix ) {
		setSliceParams_( matrix.length, matrix[0].length );
		array_ = Array( rows_ * cols_ );
		foreach( i, row ; matrix )
			foreach( j, value ; row )
				this[ i, j ] = value;
		
		return this;
	}
	
	/** Assign to a different kind of matrix */
	void assign( S )( S rhs ) if( /* isMatrix!S */ true ) {
		reset( rhs.rows, rhs.columns );
		static if( isRowMajor ) {
			foreach( i ; 0 .. rows_ )
				this[ i ][] = rhs[ i ];
		} else {
			foreach( i ; 0 .. cols_ ) {
				this[][ i ] = rhs[][ i ];
			}
		}
	}

	/** Access to elements */
	T opIndex( size_t i, size_t j ) const
	in {
		assert( i < rows_ );
		assert( j < cols_ );
	} body {
		return array_.refCountedPayload()[ map_( i, j ) ];
	}
	
	/// ditto
	void opIndexAssign( S )( S rhs, size_t i, size_t j )
	in {
		assert( i < rows_ );
		assert( j < cols_ );
	} body {
		array_.refCountedPayload()[ map_( i, j ) ] = rhs;
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t i, size_t j )
	in {
		assert( i < rows_ );
		assert( j < cols_ );
	} body {
		mixin( "array_.refCountedPayload()[ map_( i, j ) ]" ~ op ~ "= rhs;" );
	}
	
	/** Get a view of a given row. */
	VectorView!Wrapped opIndex( size_t i )
	in {
		assert( i < rows_ );
	} body {
		auto begin  = map_( i, 0 );
		auto end    = map_( i, cols_ );
		static if( isRowMajor )
			return typeof( return )( array_, begin, end, 1 );
		else
			return typeof( return )( array_, begin, end, leading_ );
	}

	/** Slicing operators. */
	SlicerProxy opSlice() {
		return SlicerProxy( &this, 0, rows_ );
	}
	
	/// ditto
	SlicerProxy opSlice( size_t rowStart, size_t rowEnd )
	in {
		assert( rowStart <= rowEnd );
		assert( rowEnd <= rows_ );
	} body {
		return SlicerProxy( &this, rowStart, rowEnd );
	}
	
	/// ditto
	void opSliceAssign( S )( S rhs ) if( is( S == T ) ) {
		array_.refCountedPayload()[] = rhs;	
	}
	
	/** Matrix transpose. */
	@property auto t() {
		return Matrix!(T, otherStorageOrder!storageOrder_, ArrayTemplate)( array_, cols_, rows_, leading_ );
	}
	
	/** Returned by opSlice, allows for slicing the 2nd dimension and range functionality. */
	static private struct SlicerProxy {
		// TODO: Use MatrixView -> elegant solution.
		private Matrix  m;
		private size_t start, end;
		
		this( Matrix *mat, size_t start, size_t end ) {
			this.start = start;
			this.end   = end;
			
			// Dirty hack: share same array with the other matrix. Basically a matrix view.
			m.array_ = mat.array_; 
			m.setSliceParams_( mat.rows_, mat.cols_, mat.leading_ );
		}
		
		/** Slicing operators. */
		Matrix opSlice() {
			return Matrix( m, start, end, 0, m.columns );
		}
		
		/// ditto
		Matrix opSlice( size_t colStart, size_t colEnd ) {

			return Matrix( m, start, end, colStart, colEnd );
		}
		
		/** Return a view to a slice of a given column. */
		auto opIndex( size_t j ) {
			auto begin  = m.map_( start, j );
			auto end    = m.map_( end,   j );
			static if( isRowMajor )
				return VectorView!Wrapped( m.array_, begin, end, m.leading_ );
			else
				return VectorView!Wrapped( m.array_, begin, end, 1 );
		}
		
		/// ditto
		void opIndexAssign( S )( S rhs, size_t i ) {
			opIndex( i )[] = rhs;
		}

		
		/** Range of ranges interface. */
		void popFront() { ++start; }
		@property bool empty() {
			static if( isRowMajor ) return start == m.rows_;
			else                    return start == m.cols_;
		}
		
		@property VectorView!Wrapped front() {
			auto a = start * m.leading_;
			return typeof(return)( m.array_, a, a + (isRowMajor?m.cols_:m.rows_), 1 );	
		}
	}
	
	@property VectorView!Wrapped front() { assert(false); }
	
	// Provides toString().
	mixin MatrixToString;

	private size_t map_( size_t row, size_t column ) const {
		static if( isRowMajor )
			return leading_ * row + column;
		else
			return leading_ * column + row;
	}

	private void setSliceParams_( size_t r, size_t c, size_t l = 0 ) {
		static if( isRowMajor )
			leading_ = l ? l : c;
		else
			leading_ = l ? l : r;
		rows_    = r;
		cols_    = c;
	}

	private size_t leading_, rows_, cols_;
	private Array array_;
}




unittest {
	void testMat(Mat)() {
		auto A = Mat(3, 3);
		Mat  B = [ [1,2], [3,4] ];
		
		assert( A.rows == 3 && A.columns == 3, "Size constructor failed." );

		assert( B[0,0] == 1 && B[0,1] == 2 &&
				B[1,0] == 3 && B[1,1] == 4, "Array of arrays ctor equality failed." );

		int k = 0;
		foreach( i ; 0..3 )
			foreach( j ; 0..3 ) {
				assert( A[i,j] == 0, "Matrix initialization failed." );
				A[i,j] = k++;
			}

		assert( A[ 0, 0 ] == 0 && A[ 0, 1 ] == 1 && A[ 0, 2 ] == 2 &&
				A[ 1, 0 ] == 3 && A[ 1, 1 ] == 4 && A[ 1, 2 ] == 5 &&
				A[ 2, 0 ] == 6 && A[ 2, 1 ] == 7 && A[ 2, 2 ] == 8, "Set failed"  );

		auto C = A[ 1 .. 3 ][ 1 .. 3 ];
		assert( C[0,0] == A[1,1] && C[0,1] == A[1,2] &&
				C[1,0] == A[2,1] && C[1,1] == A[2,2], "Slicing failed" );

		C[0, 0] = 42;
		assert( A[ 1, 1 ] == 4 && C[ 0, 0 ] == 42, "Slicing duplication failed" );
		
		C = A;
		foreach( i ; 0..3 )
			foreach( j ; 0..3 )
				assert( C[i, j] == A[i, j], "Assignment failed." );
		
		C[0, 0] = 42;
		assert( A[ 1, 1 ] == 4 && C[ 0, 0 ] == 42, "Assignment duplication failed" );
		
		auto D = C;
		D[0, 0] = 1;
		assert( C[ 0, 0 ] == 42 && D[ 0, 0 ] == 1, "Copy duplication failed" );
		
		assert( to!string(A[0][0..2]) == "[0, 1]", "Mixed slicing 1 failed." );
		assert( to!string(A[1][0..2]) == "[3, 4]", "Mixed slicing 2 failed." );
		assert( to!string(A[2][0..2]) == "[6, 7]", "Mixed slicing 3 failed." );
		assert( to!string(A[1..3][0]) == "[3, 6]", "Mixed slicing 4 failed." );
		assert( to!string(A[1..3][1]) == "[4, 7]", "Mixed slicing 5 failed." );
		assert( to!string(A[1..3][2]) == "[5, 8]", "Mixed slicing 6 failed." );
	}
	
	testMat!(Matrix!int)();
	testMat!(Matrix!(int,StorageOrder.RowMajor))();
}


mixin template MatrixToString() {
	/** (Somewhat) Nicely formatted matrix output. */
	string toString() const {
		if( rows == 0 || columns == 0 )
			return "[]";
		
		auto app = appender!string();
		app.put("[ ");
		
		void putRow( size_t r ) {
			app.put( to!string( this[r, 0] ) );
			for( int i = 1; i < columns ; ++ i ) {
				app.put( '\t' );
				app.put( to!string( this[r, i] ) );
			}
		};
		
		auto last = rows - 1;
		for( auto i = 0; i < last ; ++ i ) {
			putRow( i );
			app.put( ";\n  " );
		}
		
		putRow( last );
		app.put(" ]");

		return app.data;
	}
}