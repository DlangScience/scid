module scid.matrix;

import scid.internal.assertmessages;
import scid.ops.common;
import scid.storage.generalmat;
import scid.storage.generalmatview;
import scid.storage.triangular;
import scid.storage.symmetric;
import scid.storage.diagonalmat;
import scid.storage.cowmatrix;
import scid.storage.external;
import scid.common.traits;

import std.algorithm, std.range, std.conv;

enum MatrixTriangle {
	Upper, Lower
}

enum StorageOrder {
	RowMajor, ColumnMajor
}

enum MatrixStorageType {
	Virtual, General,
	Triangular, Symmetric, Hermitian,
	GeneralBand, SymmetricBand, HermitianBand,
	Diagonal
};

template storageOrderOf( M ) {
	static if( is( typeof(M.storageOrder) ) )
		enum storageOrderOf = M.storageOrder;
	else static if( is( M E : RefCounted!(E,a), uint a ) )
		enum storageOrderOf = storageOrderOf!E;
	else {
		static assert( false, M.stringof ~ " hasn't got a storage order. ");
	}
}

template storageTypeOf( M ) {
	static if( is( typeof(M.storageType) ) )
		enum storageTypeOf = M.storageType;
	else static if( is( M E : RefCounted!(E,a), uint a ) )
		enum storageTypeOf = storageTypeOf!E;
	else {
		static assert( false, M.stringof ~ " hasn't got a storage type. ");
	}
}

template Matrix( ElementOrStorage, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isScalar!(BaseElementType!ElementOrStorage) ) {
	
	static if( isScalar!ElementOrStorage )
		alias BasicMatrix!( GeneralMatrixStorage!(ElementOrStorage,order_) ) Matrix;
	else
		alias BasicMatrix!( ElementOrStorage )                               Matrix;
}

template MatrixView( ElementOrStorageType, StorageOrder order_ = StorageOrder.ColumnMajor ) {
	alias Matrix!( ElementOrStorageType, order_ ).View MatrixView;
}

template TriangularMatrix( ElementOrArray, MatrixTriangle tri_ = MatrixTriangle.Upper, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isScalar!(BaseElementType!ElementOrArray) ) {
	
	alias BasicMatrix!( TriangularStorage!(ElementOrArray,tri_,order_) ) TriangularMatrix;
}

template SymmetricMatrix( ElementOrArray, MatrixTriangle tri_ = MatrixTriangle.Upper, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isScalar!(BaseElementType!ElementOrArray) ) {
	
	alias BasicMatrix!( SymmetricStorage!(ElementOrArray,tri_,order_) ) SymmetricMatrix;
}

template DiagonalMatrix( ElementOrArray )
		if( isScalar!(BaseElementType!ElementOrArray) ) {
	
	alias BasicMatrix!( DiagonalMatrixStorage!(ElementOrArray ) ) DiagonalMatrix;
}

template ExternalMatrixView( ElementOrContainer, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isScalar!( BaseElementType!ElementOrContainer ) ) {
	
	static if( isScalar!ElementOrContainer ) {
		alias BasicMatrix!(
			BasicGeneralMatrixViewStorage!(
				ExternalMatrix!( ElementOrContainer, order_,
					CowMatrixRef!( ElementOrContainer, order_ )
				)
			)
		) ExternalMatrixView;
	} else {
		alias BasicMatrix!(
			BasicGeneralMatrixViewStorage!(
				ExternalMatrix!(
					BaseElementType!ElementOrContainer,
					storageOrderOf!ElementOrContainer,
					ElementOrContainer
				)
			)
		) ExternalMatrixView;
	}
}

struct BasicMatrix( Storage_ ) {
	/** The type of the matrix elements. The same as the element type of the storage */
	alias BaseElementType!Storage_ ElementType;
	
	/** The wrapped storage type. */
	alias Storage_ Storage;
	
	/** The storage order (StorageOrder.RowMajor or StorageOrder.ColumnMajor). */
	static if( is( typeof(storageOrderOf!Storage) ) )
		alias storageOrderOf!Storage   storageOrder;
	else
		alias StorageOrder.ColumnMajor storageOrder;
	
	/** The type of the vector view of a matrix row (i.e. the return type of row(i)). */
	alias typeof(Storage.init.row(0)) RowView;
	
	/** The type of the vector view of a matrix column (i.e. the return type of column(j)). */
	alias typeof(Storage.init.column(0)) ColumnView;
	
	/** The type of a view of the storage (i.e. the return type of storage.view). */
	alias typeof(Storage.init.view(0,0,0,0)) StorageView;
	
	/** The type of a slice of the storage (i.e. the return type of storage.slice). */
	alias typeof(Storage.init.slice(0,0,0,0)) StorageSlice;
	
	/** The type of a matrix view (i.e. the return type of view( rowStart, colStart, rowEnd, colEnd )). */
	alias BasicMatrix!( StorageView ) View;
	
	/** The type of a matrix slice (i.e. the return type of slice( rowStart, colStart, rowEnd, colEnd )). */
	alias BasicMatrix!( StorageSlice ) Slice;
	
	/** Inherit any special functionality from the storage type. */
	alias storage this;
	
	/** Whether the matrix is row major or not. Shorthand for checking the storageOrder member. */
	enum isRowMajor = ( storageOrder == StorageOrder.RowMajor );
	
	/** Whether the storage can be resized. */
	enum isResizable = is( typeof( Storage.init.resize(0,0) ) );
	
	/** The type obtained by transposing the matrix. If the storage doesn't define one then it's just typeof(this). */
	static if( is( Storage.Transposed ) )
		alias BasicMatrix!( Storage.Transposed ) Transposed;
	else
		alias typeof( this ) Transposed;
	
	/** The type of a temporary that can store the result of an expression that evaluates to typeof(this). If the storage doesn't define one then it's just typeof(this). */
	static if( is( Storage.Temporary ) )
		alias BasicMatrix!( Storage.Temporary ) Temporary;
	else
		alias typeof( this ) Temporary;
	
	/** The type returned by front() and back(). */
	static if( isRowMajor )
		alias RowView MajorView;
	else
		alias ColumnView MajorView;
	
	/** Promotions for matrices:
	      - Promotes vectors to vectors of the promotion of the storages.
	      - Promotes matrices to matrices of the promotion of the storages.
	      - Promotes scalars to matrices of the promotion of the storage and the scalar type.
	*/
	template Promote( T ) {
		static if( isScalar!T )	
			alias BasicMatrix!( Promotion!(Storage,T) ) Promote;
		else static if( isMatrix!T ) {
			alias BasicMatrix!( Promotion!(Storage,T.Storage) ) Promote;
		}
	}
	
	/** Create a new matrix. Forwards the arguments to the storage type. */
	this( A... )( A args ) if( A.length > 0 && !is( A[0] : Storage ) && !isMatrix!(A[0]) && !isExpression!(A[0]) ) {
		storage = Storage( args );
	}
	
	/** Create a new matrix as a copy of a matrix with the same or a different storage type. */
	this( A )( BasicMatrix!A other ) {
		static if( is( A : Storage ) ) storage = other.storage;
		else                           storage.copy( other.storage );
	}
	
	/** Create a new matrix that stores the result of an expression. */
	this( Expr )( auto ref Expr expr ) if( isExpression!Expr ) {
		this[] = expr;
	}
	
	/** Create a new matrix from a given storage. */
	this()( Storage stor ) {
		storage = stor;
	}
	
	/** Element access. Forwarded to the storage.index method. */
	ElementType opIndex( size_t i, size_t j ) const {
		return storage.index( i, j );
	}
	
	/** Element access. Forwarded to the storage.indexAssign method. */
	void opIndexAssign( ElementType rhs, size_t i, size_t j ) {
		storage.indexAssign( rhs, i, j );
	}
	
	/** Element access. Forwarded to the storage.indexAssign!op method. */
	void opIndexOpAssign( string op )( ElementType rhs, size_t i, size_t j ) {
		storage.indexAssign!op( rhs, i, j );	
	}
	
	/** Get a view of the matrix. Forwarded to the storage.view method. */
	View view( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( storage.view( rowStart, rowEnd, colStart, colEnd ) );
	}
	
	/** Matrix slicing. Forwarded to the storage.slice method. */
	Slice slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( storage.slice( rowStart, rowEnd, colStart, colEnd ) );
	}
	
	/** Get a view of the i'th row. Forwarded to the storage.row method. */
	RowView row( size_t i ) {
		return storage.row( i );
	}
	
	/** Get a view of the j'th column. Forwarded to the storage.column method. */
	ColumnView column( size_t j ) {
		return storage.column( j );	
	}
	
	/** Get a view of a slice of the i'th row. Forwarded to the storage.rowSlice method if available.
	    Otherwise forward to storage.row(i)[ start .. end ].
	*/
	auto rowSlice( size_t i, size_t start, size_t end ) {
		static if( is(typeof(storage.rowSlice( i, start, end ))) )
			return storage.rowSlice( i, start, end );
		else
			return storage.row(i)[ start .. end ];
	}
	
	/** Get a view of a slice of the j'th column. Forwarded to the storage.columnSlice method if available.
	    Otherwise forward to storage.column(i)[ start .. end ].
	*/
	auto columnSlice( size_t i, size_t start, size_t end ) {
		static if( is(typeof(storage.columnSlice( i, start, end ))) )
			return storage.columnSlice( i, start, end );
		else
			return storage.column(i)[ start .. end ];
	}
	
	/** Same as row(i). This is to allow a consistent slicing syntax. */
	RowView opIndex( size_t i ) {
		return row( i );
	}
	
	/** Start a slice operation. Select all rows. */
	SliceProxy_ opSlice() {
		return typeof(return)( this, 0, this.rows );
	}
	
	/** Start a slice operation. Select rows from rowStart to rowEnd (not included). */
	SliceProxy_ opSlice( size_t rowStart, size_t rowEnd ) {
		return typeof(return)( this, rowStart, rowEnd );
	}
	
	/** Assign the result of a matrix operation to this one. opAssign is not used due to a bug. */
	void opSliceAssign(Rhs)( auto ref Rhs rhs ) {
		static if( is( Rhs E : E[][] ) && isConvertible!(E, ElementType) )
			evalCopy( BasicMatrix( rhs ), this );
		else
			evalCopy( rhs, this );	
	}
	
	/** Resize the matrix and leave the memory uninitialized. If not resizeable simply check that the dimensions are
	    correct.
	*/
	void resize( size_t newRows, size_t newColumns, void* ) {
		static if( isResizable ) {
			storage.resize( newRows, newColumns, null );
		} else {
			assert( rows == newRows && columns == newColumns,
				dimMismatch_( newRows, newColumns ) );
		}
	}
	
	/** Resize the matrix and set all the elements to zero. If not resizeable, check that the dimensions are correct
	    and just set the elements to zero.
	*/
	void resize( size_t newRows, size_t newColumns ) {
		static if( isResizable ) {
			storage.resize( newRows, newColumns );
		} else {
			this.resize( newRows, newColumns, null );
			evalScaling( Zero!ElementType, this );
		}
	}
	
	/// ditto
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if( op == "+" || op == "-" || op == "*" || op == "/" ) {
		static if( op == "+" ) {
			evalScaledAddition( One!ElementType, rhs, this );
		} else static if( op == "-" ) {
			evalScaledAddition( MinusOne!ElementType, rhs, this );
		} else static if( op == "*" ) {
			static if( closureOf!Rhs == Closure.Matrix )
				this[] = this * rhs;
			else {
				static assert( isConvertible!(Rhs, ElementType), "Matrix multiplication with invalid type " ~ Rhs.stringof );
				evalScaling( to!ElementType(rhs), this );
			}
		} else /* static if ( op == "/" ) */ {
			static assert( isConvertible!(Rhs, ElementType), "Matrix division with invalid type " ~ Rhs.stringof );
			evalScaling( One!ElementType / to!ElementType(rhs), this );
		}
	}
	
	/** If provided, forward popFront() to the storage type. */
	static if( is( typeof( Storage.init.popFront() ) ) )
		void popFront() { storage.popFront(); }
	
	/** If provided, forward popBack() to the storage type. */
	static if( is( typeof( Storage.init.popBack() ) ) )
		void popBack()  { storage.popBack(); }
			
	@property {		
		/** Return the first row for row-major matrices or the first column for column-major matrices. */
		MajorView front()
		in {
			assert( !empty(), msgPrefix_ ~ "front on empty." );
		} body {
			static if( isRowMajor )
				return storage.row( 0 );
			else
				return storage.column( 0 );
		}
		
		/** Return the last row for row-major matrices or the last column for column-major matrices. */
		MajorView back()
		in {
			assert( !empty(), msgPrefix_ ~ "back on empty." );
		} body {
			static if( isRowMajor )
				return storage.row( storage.rows - 1 );
			else
				return storage.column( storage.columns - 1 );
		}
		
		/** Return the number of rows of the matrix. */
		size_t rows() const {
			static if( is(typeof(Storage.init.rows) ) )
				return storage.rows;
			else static if( is(typeof(Storage.init.size) ) )
				return storage.size;
			else
				static assert( false, "Matrix storage " ~ Storage.stringof ~ " doesn't define define dimensions." );
		}
		
		/** Return the number of columns of the matrix. */
		size_t columns() const {
			static if( is(typeof(Storage.init.columns) ) )
				return storage.columns;
			else static if( is(typeof(Storage.init.size) ) )
				return storage.size;
			else
				static assert( false, "Matrix storage " ~ Storage.stringof ~ " doesn't define define dimensions." );
		}
		
		/** Return the major dimension of the matrix. */
		size_t major() const {
			static if( isRowMajor )
				return storage.rows;
			else
				return storage.columns;
		}
		
		/** Return the minor dimension of the matrix. */
		size_t minor() const {
			static if( isRowMajor )
				return storage.columns;
			else
				return storage.rows;
		}
		
		/** For square matrices return the size of the matrix. */
		static if( is(typeof(Storage.init.size)) ) {
			size_t size() const {
				return storage.size;
			}
		}
		
		/** Has the matrix zero dimensions. */
		bool empty() const {
			return major == 0;
		}
		
		/** The number of rows or columns of the matrix depending on the storage order. */
		size_t length() const {
			return this.major;
		}
		
		/** (Somewhat) Nicely formatted matrix output. */
		string pretty() const {
			if( this.rows == 0 || this.columns == 0 )
				return "[]";
		
			auto app = appender!string();
			app.put("[ ");
		
			void putRow( size_t r ) {
				app.put( to!string( this[r, 0] ) );
				for( int i = 1; i < this.columns ; ++ i ) {
					app.put( '\t' );
					app.put( to!string( this[r, i] ) );
				}
			};
		
			auto last = this.rows - 1;
			for( auto i = 0; i < last ; ++ i ) {
				putRow( i );
				app.put( ";\n  " );
			}
		
			putRow( last );
			app.put(" ]");

			return app.data;
		}
	}
	
	/** In line matrix output. */
	string toString() const {
		if( this.rows == 0 || this.columns == 0 )
			return "[]";
		
		auto app = appender!string();
		app.put('[');
		
		void putRow( size_t r ) {
			app.put( '[' );
			app.put( to!string( this[r, 0] ) );
			for( int i = 1; i < this.columns ; ++ i ) {
				app.put( ", " );
				app.put( to!string( this[r, i] ) );
			}
			app.put( ']' );
		};
		
		auto last = this.rows - 1;
		for( auto i = 0; i < last ; ++ i ) {
			putRow( i );
			app.put( ", " );
		}
		
		putRow( last );
		app.put("]");

		return app.data;
	}
	
	/** The wrapped storage object. */
	Storage storage;
	
	/** Include operator overloads for matrices. */
	mixin Operand!( Closure.Matrix );
	
	/** Proxy type used for slicing operations. Returned after the first slicing is performed. */
	static struct SliceProxy_ {
		alias BasicMatrix MatrixType;
		alias Storage     StorageType;
		
		BasicMatrix m_;
		size_t  start_, end_;
		
		this( BasicMatrix mat, size_t start, size_t end ) {
			m_.storage.forceRefAssign( mat.storage );
			start_ = start; end_ = end;
		}
		
		Slice      opSlice()                       { return m_.slice( start_, end_, 0, m_.columns ); }
		Slice      opSlice( size_t j1, size_t j2 ) { return m_.slice( start_, end_, j1, j2 ); }
		ColumnView opIndex( size_t j )             { return m_.columnSlice( j, start_, end_ ); }
		
		void opSliceAssign( Rhs )( auto ref Rhs rhs ) {
			m_.view( start_, end_, 0, m_.columns )[] = rhs;
		}
		
		void opSliceAssign( Rhs )( auto ref Rhs rhs, size_t j1, size_t j2 ) {
			m_.view( start_, end_, j1, j2 )[] = rhs;
		}
		
		void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) {
			m_.view( start_, end_, 0, m_.columns ).opSliceOpAssign!op( rhs );
		}
		
		void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t j1, size_t j2 ) {
			m_.view( start_, end_, j1, j2 ).opSliceOpAssign!op( rhs );
		}
		
		void opIndexAssign( Rhs )( auto ref Rhs rhs, size_t j ) {
			this[ j ].opSliceAssign( rhs );
		}
		
		void opIndexOpAssign( string op, Rhs )( auto ref Rhs rhs, size_t j ) {
			this[ j ].opSliceOpAssign!op( rhs );
		}
		
		void popFront()
		in {
			assert( !empty,	"popFront on empty" );
		} body {
			++ start_;	
		}
		
		@property {
			bool empty() const { return start_ == m_.major; }
			
			MajorView  front() {
				static if( isRowMajor ) return m_.row( start_ );
				else                    return m_.column( start_ );
			}
			
		}	
	}
	
private:
	mixin MatrixErrorMessages;
}

/** Matrix inversion. */
Expression!( "inv", M ) inv( M )( auto ref M matrix ) if( closureOf!M == Closure.Matrix ) {
	return typeof( return )( matrix );	
}


version( unittest ) {
	import std.stdio;
	import scid.storage.cowmatrix;
	import scid.vector;
}


//------------------------------------- Tests for General Matrices -------------------------------------//
unittest {
	alias Matrix!( double, StorageOrder.RowMajor    ) RowMat;
	alias Matrix!( double, StorageOrder.ColumnMajor ) ColMat;
	
	static assert( is( Matrix!double : ColMat ) );
	static assert( !is( Matrix!int ) );
	static assert( is( Matrix!( GeneralMatrixStorage!(CowMatrixRef!double) ) : ColMat ) );

	// empty matrix;
	{
		RowMat a;
		ColMat b;
		
		// assert( a.rows == 0 && a.columns == 0 );
		// assert( b.rows == 0 && b.columns == 0 );
		// ops on it are undefined //
	}
	
	// matrix ctors
	{
		auto a = ColMat(3,6);
		auto b = RowMat(5,4, null);
		auto c = ColMat(3, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] );
		auto d = RowMat(3, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] );
		
		assert( a.rows == 3 && a.columns == 6 );
		assert( b.rows == 5 && b.columns == 4 );
		assert( c.rows == 3 && c.columns == 3 );
		assert( d.rows == 3 && d.columns == 3 );
		
		foreach( i ; 0 .. 3 )
			foreach( j ; 0 .. 6 )
				assert( a[ i, j ] == 0.0 );
		
		writeln( c.pretty );
		double k = 1.0;
		foreach( i ; 0 .. 3 )
			foreach( j ; 0 .. 3 ) {
				assert( c[ j, i ] == k );
				assert( d[ i, j ] == k );
				k ++;
			}
	}
	
	// copy-on-write
	void cowTest( Mat )() {
		auto a = Mat( 3, 3 );
		auto b = a;
		Mat c; c = b;
		b[ 0, 0 ] =  1.0;
		c[ 0, 0 ] += 2.0;
		assert( a[ 0, 0 ] == 0.0 && b[ 0, 0 ] == 1.0 && c[ 0, 0 ] == 2.0 );
	}
	
	cowTest!RowMat();
	cowTest!ColMat();
	
	// row/column indexing
	void rowColTest( Mat )() {
		auto a = Mat( 3, 5 );
		
		double k = 0.0;
		foreach( i ; 0 .. a.rows )
			foreach( j ; 0 .. a.columns )
				a[ i, j ] = k ++;
		
		foreach( i ; 0 .. a.rows ) {
			auto v1 = a.row( i );
			auto v2 = a[ i ][];
			
			static assert( v1.vectorType == VectorType.Row );
			static assert( v2.vectorType == VectorType.Row );
			assert( v1.length == a.columns && v2.length == a.columns );
			
			foreach( j ; 0 .. a.columns )
				assert( v1[ j ] == v2[ j ] && v1[ j ] == i * a.columns + j );
		}
		
		foreach( j ; 0 .. a.columns ) {
			auto v1 = a.column( j );
			auto v2 = a[][ j ];
			
			static assert( v1.vectorType == VectorType.Column );
			static assert( v2.vectorType == VectorType.Column );
			assert( v1.length == a.rows && v2.length == a.rows );
			
			foreach( i ; 0 .. a.rows )
				assert( v1[ i ] == v2[ i ] && v1[ i ] == i * a.columns + j );
		}
		
		foreach( i ; 0 .. a.rows ) {
			auto v = a[ i ][];
			v[] *= 2.0;
		}
		
		k = 0.0;
		foreach( i ; 0 .. a.rows )
			foreach( j ; 0 .. a.columns ) {
				assert( a[ i, j ] == k );
				k += 2.0;
			}
	}
	
	rowColTest!RowMat();
	rowColTest!ColMat();
	
	// slices & views tests
	void sliceTest( Mat )() {
		auto a = Mat( 3, 5 );
		
		double k = 0.0;
		foreach( i ; 0 .. a.rows )
			foreach( j ; 0 .. a.columns )
				a[ i, j ] = ++ k;
		
		
		auto p = a[0..2][3..5];
		assert( p.rows == 2 && p.columns == 2 );
		assert( p[0, 0] == 4.0 && p[0, 1] == 5.0 &&
			    p[1, 0] == 9.0 && p[1, 1] == 10.0 );
		
		p[ 0, 0 ] = 42.0;
		assert( a[ 0, 3 ] == 4.0 && p[ 0, 0 ] == 42.0 );
		
		auto q = a.view( 1, 3, 2, 5 ); 
		auto r = q[ 0 .. 2 ][ 0 .. 3 ];
		assert( q.rows == 2 && q.columns == 3 );
		foreach( i ; 0 .. q.rows )
			foreach( j ; 0 .. q.columns )
				assert( q[ i, j ] == r[  i, j ] );
		
		q[ 0, 0 ] = 42.0;
		assert( a[ 1, 2 ] == 42.0 && q[ 0, 0 ] == 42.0 && r[ 0, 0 ] == 42.0 );
	}
	
	sliceTest!RowMat();
	sliceTest!ColMat();
	
	// range of ranges interface
	void rangeTest(Mat)() {
		auto a = Mat( 3, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] );
		auto m = [ [1.0, 2.0], [3.0, 4.0], [5.0, 6.0] ];
		
		uint k = 0;
		foreach( vec ; a ) {
		    assert( m[ k ][ 0 ] == vec[ 0 ] && m[ k ][ 1 ] == vec[ 1 ] );
		    ++k;
		}
		
		k = 0;
		foreach( vec ; a[0..2][] ) {
		    assert( m[ k ][ 0 ] == vec[ 0 ] && m[ k ][ 1 ] == vec[ 1 ] );
		    ++k;
		}
	}
	
	rangeTest!RowMat();
	rangeTest!ColMat();
}

//------------------------------------ Tests for Triangular Matrices -----------------------------------//
unittest {
	alias TriangularMatrix!( double, MatrixTriangle.Upper, StorageOrder.RowMajor    ) RowUpMat;
	alias TriangularMatrix!( double, MatrixTriangle.Upper, StorageOrder.ColumnMajor ) ColUpMat;
	alias TriangularMatrix!( double, MatrixTriangle.Lower, StorageOrder.RowMajor    ) RowLoMat;
	alias TriangularMatrix!( double, MatrixTriangle.Lower, StorageOrder.ColumnMajor ) ColLoMat;
	
	static assert( is( TriangularMatrix!double : ColUpMat ) );
	static assert( !is( TriangularMatrix!int ) );

	// empty matrix;
	{
		RowUpMat a;
		ColUpMat b;
		
		// assert( a.rows == 0 && a.columns == 0 );
		// assert( b.rows == 0 && b.columns == 0 );
		// ops on it are undefined //
	}
	
	// matrix ctors
	{
		auto cu = ColUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto ru = RowUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto cl = ColLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto rl = RowLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		
		assert( cu.size == 3 && ru.size == 3 && cl.size == 3 && rl.size == 3 );
		
		assert( cu[ 0, 0 ] == 1.0 && cu[ 0, 1 ] == 2.0 && cu[ 0, 2 ] == 4.0 &&
			    cu[ 1, 0 ] == 0.0 && cu[ 1, 1 ] == 3.0 && cu[ 1, 2 ] == 5.0 &&
			    cu[ 2, 0 ] == 0.0 && cu[ 2, 1 ] == 0.0 && cu[ 2, 2 ] == 6.0 );
		
		assert( ru[ 0, 0 ] == 1.0 && ru[ 0, 1 ] == 2.0 && ru[ 0, 2 ] == 3.0 &&
			    ru[ 1, 0 ] == 0.0 && ru[ 1, 1 ] == 4.0 && ru[ 1, 2 ] == 5.0 &&
			    ru[ 2, 0 ] == 0.0 && ru[ 2, 1 ] == 0.0 && ru[ 2, 2 ] == 6.0 );
		
		assert( cl[ 0, 0 ] == 1.0 && cl[ 0, 1 ] == 0.0 && cl[ 0, 2 ] == 0.0 &&
			    cl[ 1, 0 ] == 2.0 && cl[ 1, 1 ] == 4.0 && cl[ 1, 2 ] == 0.0 &&
			    cl[ 2, 0 ] == 3.0 && cl[ 2, 1 ] == 5.0 && cl[ 2, 2 ] == 6.0 );
		
		assert( rl[ 0, 0 ] == 1.0 && rl[ 0, 1 ] == 0.0 && rl[ 0, 2 ] == 0.0 &&
			    rl[ 1, 0 ] == 2.0 && rl[ 1, 1 ] == 3.0 && rl[ 1, 2 ] == 0.0 &&
			    rl[ 2, 0 ] == 4.0 && rl[ 2, 1 ] == 5.0 && rl[ 2, 2 ] == 6.0 );
		
	}
	
	// copy-on-write
	void cowTest( Mat )() {
		auto a = Mat( 3 );
		auto b = a;
		Mat c; c = b;
		b[ 0, 0 ] =  1.0;
		c[ 0, 0 ] += 2.0;
		assert( a[ 0, 0 ] == 0.0 && b[ 0, 0 ] == 1.0 && c[ 0, 0 ] == 2.0 );
	}
	
	cowTest!RowUpMat();
	cowTest!ColUpMat();
	cowTest!RowLoMat();
	cowTest!ColLoMat();
	
	// row/column indexing
	{
		auto cu = ColUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto ru = RowUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto cl = ColLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto rl = RowLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		
		{
			auto r = cu.row( 2 );
			r[] *= 2.0;
		}
		
		assert( cu[ 0, 0 ] == 1.0 && cu[ 0, 1 ] == 2.0 && cu[ 0, 2 ] == 4.0 &&
			    cu[ 1, 0 ] == 0.0 && cu[ 1, 1 ] == 3.0 && cu[ 1, 2 ] == 5.0 &&
			    cu[ 2, 0 ] == 0.0 && cu[ 2, 1 ] == 0.0 && cu[ 2, 2 ] == 12.0 );

		{
			auto r = ru.column( 1 );
			r[] *= 2.0;
		}
		
		assert( ru[ 0, 0 ] == 1.0 && ru[ 0, 1 ] == 4.0 && ru[ 0, 2 ] == 3.0 &&
			    ru[ 1, 0 ] == 0.0 && ru[ 1, 1 ] == 8.0 && ru[ 1, 2 ] == 5.0 &&
			    ru[ 2, 0 ] == 0.0 && ru[ 2, 1 ] == 0.0 && ru[ 2, 2 ] == 6.0 );
		
		cl[][1] = [9.0,8.0,7.0];
		assert( cl[ 0, 0 ] == 1.0 && cl[ 0, 1 ] == 0.0 && cl[ 0, 2 ] == 0.0 &&
			    cl[ 1, 0 ] == 2.0 && cl[ 1, 1 ] == 8.0 && cl[ 1, 2 ] == 0.0 &&
			    cl[ 2, 0 ] == 3.0 && cl[ 2, 1 ] == 7.0 && cl[ 2, 2 ] == 6.0 );
		
		rl[1][] = [9.0,8.0,7.0].t;
		assert( rl[ 0, 0 ] == 1.0 && rl[ 0, 1 ] == 0.0 && rl[ 0, 2 ] == 0.0 &&
			    rl[ 1, 0 ] == 9.0 && rl[ 1, 1 ] == 8.0 && rl[ 1, 2 ] == 0.0 &&
			    rl[ 2, 0 ] == 4.0 && rl[ 2, 1 ] == 5.0 && rl[ 2, 2 ] == 6.0 );
		
	}
	
	// slices & views tests
	{
		auto a = ColUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		
		auto p = a[0..2][1..3];
		assert( p.rows == 2 && p.columns == 2 );
		assert( p[0, 0] == 2.0 && p[0, 1] == 4.0 &&
			    p[1, 0] == 3.0 && p[1, 1] == 5.0 );
		
		p[ 0, 0 ] = 42.0;
		assert( a[ 0, 1 ] == 2.0 && p[ 0, 0 ] == 42.0 );
		
		auto q = a.view( 1, 3, 0, 3 ); 
		auto r = q[ 0 .. 2 ][ 0 .. 3 ];

		assert( q.rows == 2 && q.columns == 3 );
		foreach( i ; 0 .. q.rows )
			foreach( j ; 0 .. q.columns )
				assert( q[ i, j ] == r[  i, j ] );
		
		q[ 0, 1 ] = 42.0;
		assert( a[ 1, 1 ] == 42.0 && q[ 0, 1 ] == 42.0 && r[ 0, 1 ] == 42.0 );
	}
	
	// range of ranges interface
	{
		auto a = RowLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] );
		auto m = [ [1.0, 0.0, 0.0], [2.0, 3.0, 0.0], [4.0, 5.0, 6.0] ];
		
		uint k = 0;
		foreach( vec ; a ) {
			assert( m[ k ][ 0 ] == vec[ 0 ] && m[ k ][ 1 ] == vec[ 1 ] );
			++ k;
		}
		
		k = 0;
		foreach( vec ; a[0..2][] ) {
			assert( m[ k ][ 0 ] == vec[ 0 ] && m[ k ][ 1 ] == vec[ 1 ] );
			++ k;
		}
	}
}

//------------------------------------- Tests for Symmetric Matrices -----------------------------------//
unittest {
	alias SymmetricMatrix!( double, MatrixTriangle.Upper, StorageOrder.RowMajor    ) RowUpMat;
	alias SymmetricMatrix!( double, MatrixTriangle.Upper, StorageOrder.ColumnMajor ) ColUpMat;
	alias SymmetricMatrix!( double, MatrixTriangle.Lower, StorageOrder.RowMajor    ) RowLoMat;
	alias SymmetricMatrix!( double, MatrixTriangle.Lower, StorageOrder.ColumnMajor ) ColLoMat;
	
	static assert( is( SymmetricMatrix!double : ColUpMat ) );
	static assert( !is( SymmetricMatrix!int ) );

	// empty matrix;
	{
		RowUpMat a;
		ColUpMat b;
		
		// assert( a.rows == 0 && a.columns == 0 );
		// assert( b.rows == 0 && b.columns == 0 );
		// ops on it are undefined //
	}
	
	// matrix ctors
	{
		auto cu = ColUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto ru = RowUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto cl = ColLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto rl = RowLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		
		assert( cu.size == 3 && ru.size == 3 && cl.size == 3 && rl.size == 3 );
		
		assert( cu[ 0, 0 ] == 1.0 && cu[ 0, 1 ] == 2.0 && cu[ 0, 2 ] == 4.0 &&
			    cu[ 1, 0 ] == 2.0 && cu[ 1, 1 ] == 3.0 && cu[ 1, 2 ] == 5.0 &&
			    cu[ 2, 0 ] == 4.0 && cu[ 2, 1 ] == 5.0 && cu[ 2, 2 ] == 6.0 );
		
		assert( ru[ 0, 0 ] == 1.0 && ru[ 0, 1 ] == 2.0 && ru[ 0, 2 ] == 3.0 &&
			    ru[ 1, 0 ] == 2.0 && ru[ 1, 1 ] == 4.0 && ru[ 1, 2 ] == 5.0 &&
			    ru[ 2, 0 ] == 3.0 && ru[ 2, 1 ] == 5.0 && ru[ 2, 2 ] == 6.0 );
		
		assert( cl[ 0, 0 ] == 1.0 && cl[ 0, 1 ] == 2.0 && cl[ 0, 2 ] == 3.0 &&
			    cl[ 1, 0 ] == 2.0 && cl[ 1, 1 ] == 4.0 && cl[ 1, 2 ] == 5.0 &&
			    cl[ 2, 0 ] == 3.0 && cl[ 2, 1 ] == 5.0 && cl[ 2, 2 ] == 6.0 );
		
		assert( rl[ 0, 0 ] == 1.0 && rl[ 0, 1 ] == 2.0 && rl[ 0, 2 ] == 4.0 &&
			    rl[ 1, 0 ] == 2.0 && rl[ 1, 1 ] == 3.0 && rl[ 1, 2 ] == 5.0 &&
			    rl[ 2, 0 ] == 4.0 && rl[ 2, 1 ] == 5.0 && rl[ 2, 2 ] == 6.0 );
		
	}
	
	// copy-on-write
	void cowTest( Mat )() {
		auto a = Mat( 3 );
		auto b = a;
		Mat c; c = b;
		b[ 0, 0 ] =  1.0;
		c[ 0, 0 ] += 2.0;
		assert( a[ 0, 0 ] == 0.0 && b[ 0, 0 ] == 1.0 && c[ 0, 0 ] == 2.0 );
	}
	
	cowTest!RowUpMat();
	cowTest!ColUpMat();
	cowTest!RowLoMat();
	cowTest!ColLoMat();
	
	// row/column indexing
	{
		auto cu = ColUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto ru = RowUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto cl = ColLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		auto rl = RowLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		
		{
			auto r = cu.row( 2 );
			r[] *= 2.0;
		}
		
		assert( cu[ 0, 0 ] == 1.0 && cu[ 0, 1 ] == 2.0 && cu[ 0, 2 ] == 4.0 &&
			    cu[ 1, 0 ] == 2.0 && cu[ 1, 1 ] == 3.0 && cu[ 1, 2 ] == 5.0 &&
			    cu[ 2, 0 ] == 4.0 && cu[ 2, 1 ] == 5.0 && cu[ 2, 2 ] == 12.0 );

		{
			auto r = ru.column( 1 );
			r[] *= 2.0;
		}
		
		assert( ru[ 0, 0 ] == 1.0 && ru[ 0, 1 ] == 4.0 && ru[ 0, 2 ] == 3.0 &&
			    ru[ 1, 0 ] == 4.0 && ru[ 1, 1 ] == 8.0 && ru[ 1, 2 ] == 5.0 &&
			    ru[ 2, 0 ] == 3.0 && ru[ 2, 1 ] == 5.0 && ru[ 2, 2 ] == 6.0 );
		
		cl[][1] = [9.0,8.0,7.0];
		assert( cl[ 0, 0 ] == 1.0 && cl[ 0, 1 ] == 2.0 && cl[ 0, 2 ] == 3.0 &&
			    cl[ 1, 0 ] == 2.0 && cl[ 1, 1 ] == 8.0 && cl[ 1, 2 ] == 7.0 &&
			    cl[ 2, 0 ] == 3.0 && cl[ 2, 1 ] == 7.0 && cl[ 2, 2 ] == 6.0 );
		
		rl[1][] = [9.0,8.0,7.0].t;
		assert( rl[ 0, 0 ] == 1.0 && rl[ 0, 1 ] == 9.0 && rl[ 0, 2 ] == 4.0 &&
			    rl[ 1, 0 ] == 9.0 && rl[ 1, 1 ] == 8.0 && rl[ 1, 2 ] == 5.0 &&
			    rl[ 2, 0 ] == 4.0 && rl[ 2, 1 ] == 5.0 && rl[ 2, 2 ] == 6.0 );
		
	}
	
	// slices & views tests
	{
		auto a = ColUpMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0 ] );
		
		auto p = a[0..2][1..3];
		assert( p.rows == 2 && p.columns == 2 );
		assert( p[0, 0] == 2.0 && p[0, 1] == 4.0 &&
			    p[1, 0] == 3.0 && p[1, 1] == 5.0 );
		
		p[ 0, 0 ] = 42.0;
		assert( a[ 0, 1 ] == 2.0 && p[ 0, 0 ] == 42.0 );
		
		auto q = a.view( 1, 3, 0, 3 ); 
		auto r = q[ 0 .. 2 ][ 0 .. 3 ];

		assert( q.rows == 2 && q.columns == 3 );
		foreach( i ; 0 .. q.rows )
			foreach( j ; 0 .. q.columns )
				assert( q[ i, j ] == r[  i, j ] );
		
		q[ 0, 1 ] = 42.0;
		assert( a[ 1, 1 ] == 42.0 && q[ 0, 1 ] == 42.0 && r[ 0, 1 ] == 42.0 );
	}
	
	// range of ranges interface
	{
		auto a = RowLoMat( [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] );
		auto m = [ [1.0, 2.0, 4.0], [2.0, 3.0, 5.0], [4.0, 5.0, 6.0] ];
		
		uint k = 0;
		foreach( vec ; a ) {
			assert( m[ k ][ 0 ] == vec[ 0 ] && m[ k ][ 1 ] == vec[ 1 ] );
			++ k;
		}
		
		k = 0;
		foreach( vec ; a[0..2][] ) {
			assert( m[ k ][ 0 ] == vec[ 0 ] && m[ k ][ 1 ] == vec[ 1 ] );
			++ k;
		}
	}
}
