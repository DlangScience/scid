module scid.matrix;

import scid.internal.assertmessages;
import scid.internal.expression;
import scid.internal.hlblas;
import scid.storage.generalmat;
import scid.storage.triangular;
import scid.storage.symmetric;
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
		if( isFortranType!(BaseElementType!ElementOrStorage) ) {
	
	static if( isFortranType!ElementOrStorage )
		alias BasicMatrix!( GeneralMatrixStorage!(ElementOrStorage,order_) ) Matrix;
	else
		alias BasicMatrix!( ElementOrStorage )                               Matrix;
}

template MatrixView( ElementOrStorageType, StorageOrder order_ = StorageOrder.ColumnMajor ) {
	alias Matrix!( ElementOrStorageType, order_ ).View MatrixView;
}

template TriangularMatrix( ElementOrArray, MatrixTriangle tri_ = MatrixTriangle.Upper, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	alias BasicMatrix!( TriangularStorage!(ElementOrArray,tri_,order_) ) TriangularMatrix;
}

template SymmetricMatrix( ElementOrArray, MatrixTriangle tri_ = MatrixTriangle.Upper, StorageOrder order_ = StorageOrder.ColumnMajor )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	alias BasicMatrix!( SymmetricStorage!(ElementOrArray,tri_,order_) ) SymmetricMatrix;
}


template isMatrix( T ) {
	static if( is( typeof( T.init[0,0]          ) ) &&
			   is( typeof( T.init[0..1][0..1]   ) ) &&
			   is( typeof( T.init.storage       ) ) &&
			   is( typeof( T.init.view(0,0,0,0) ) ) )
		enum isMatrix = true;
	else	   
		enum isMatrix = false;
}

struct BasicMatrix( Storage_ ) {
	alias BaseElementType!Storage_      ElementType;
	alias Storage_                      Storage;
	alias BasicMatrix!( Storage.Slice ) Slice;
	alias Storage.storageOrder          storageOrder;
	alias Storage.RowView               RowView;
	alias Storage.ColumnView            ColumnView;
	alias Storage.DiagonalView          DiagonalView;
	alias BasicMatrix!( Storage.View )  View;
	alias storage                       this;
	
	
	enum isRowMajor = ( storageOrder == StorageOrder.RowMajor );
	
	static if( isRowMajor ) {
		alias RowView    MajorView;
		alias ColumnView MinorView;
	} else {
		alias RowView    MinorView;
		alias ColumnView MajorView;
	}
	
	this( A... )( A args ) if( A.length > 0 && !is( A[0] : Storage ) && !isMatrix!(A[0]) && !isExpression!(A[0]) ) {
		storage = Storage( args );
	}
	
	this( A )( BasicMatrix!(A, vectorType) other ) {
		static if( is( A : Storage ) ) storage = other.storage;
		else                           storage.copy( other.storage );
	}
	
	this( Expr )( auto ref Expr expr ) if( isExpression!Expr ) {
		this[] = expr;
	}
	
	this( this ) {}
	
	this()( Storage stor ) {
		move( stor, storage );
	}
	
	ElementType opIndex( size_t i, size_t j ) const {
		return storage.index( i, j );
	}
	
	void opIndexAssign( ElementType rhs, size_t i, size_t j ) {
		storage.indexAssign( rhs, i, j );
	}
	
	void opIndexOpAssign( string op )( ElementType rhs, size_t i, size_t j ) {
		storage.indexAssign!op( rhs, i, j );	
	}
	
	View view( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( storage.view( rowStart, rowEnd, colStart, colEnd ) );
	}
	
	Slice slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( storage.slice( rowStart, rowEnd, colStart, colEnd ) );
	}
	
	RowView row( size_t i ) {
		return storage.row( i );
	}
	
	ColumnView column( size_t j ) {
		return storage.column( j );	
	}
	
	/*
	DiagonalView diagonal( size_t k ) {
		return storage.diagonal( k );
	}
	*/
	
	RowView opIndex( size_t i ) {
		return row( i );
	}
	
	SliceProxy_ opSlice() {
		return typeof(return)( this, 0, this.rows );
	}
	
	SliceProxy_ opSlice( size_t rowStart, size_t rowEnd ) {
		return typeof(return)( this, rowStart, rowEnd );
	}
	
	void opSliceAssign(Rhs)( auto ref Rhs rhs ) {
		hlCopy( rhs, this );	
	}
	
	void opSliceOpAssign( string op, Rhs )( auto ref Rhs rhs ) if( op == "+" || op == "*" || op == "/" ) {
		static if( op == "+" )      { hlAxpy( One!ElementType, rhs, this ); }
		else static if( op == "*" ) {
			static if( closureOf!Rhs == Closure.Matrix ) this[] = this * rhs;
			else                                         hlScal( rhs, this );
		} else static if ( op == "/" ) {
			hlScal( One!ElementType	/ rhs, this );
		}
	}
			
	@property {
		MajorView front() { return storage.front; }
		MajorView back()  { return storage.back; }
		
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
	
	Storage storage;
	
	mixin Literal!( Closure.Matrix );
	
	static struct SliceProxy_ {
		alias BasicMatrix MatrixType;
		alias Storage     StorageType;
		
		BasicMatrix m_;
		size_t  start_, end_;
		
		this( ref BasicMatrix mat, size_t start, size_t end ) {
			m_.forceRefAssign( mat );
			start_ = start; end_ = end;
		}
		
		Slice      opSlice()                       { return m_.slice( start_, end_, 0, m_.columns ); }
		Slice      opSlice( size_t j1, size_t j2 ) { return m_.slice( start_, end_, j1, j2 ); }
		ColumnView opIndex( size_t j )             { return m_.columnSlice( j, start_, end_ ); }
		
		void opSliceAssign( Rhs )( auto ref Rhs rhs ) {
			hlCopy( rhs, m_.view( start_, end_, 0, m_.columns ) );
		}
		
		void opSliceAssign( Rhs )( auto ref Rhs rhs, size_t j1, size_t j2 ) {
			hlCopy( rhs, m_.view( start_, end_, j1, j2 ) );
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
	
	
}

version( unittest ) {
	import std.stdio;
	import scid.storage.cowmatrix;
	import scid.vector;
}


/************************************** Tests for General Matrices **************************************/
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
		auto a = ColMat(3,6, 1.0);
		auto b = RowMat(5,4, null);
		auto c = ColMat(3, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] );
		auto d = RowMat(3, [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0] );
		
		assert( a.rows == 3 && a.columns == 6 );
		assert( b.rows == 5 && b.columns == 4 );
		assert( c.rows == 3 && c.columns == 3 );
		assert( d.rows == 3 && d.columns == 3 );
		
		foreach( i ; 0 .. 3 )
			foreach( j ; 0 .. 6 )
				assert( a[ i, j ] == 1.0 );
		
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
		auto a = Mat( 3, 3, 0.0 );
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
			v *= 2.0;
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

/************************************** Tests for Triangular Matrices **************************************/
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
		auto a = Mat( 3, 0.0 );
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
			r *= 2.0;
		}
		
		assert( cu[ 0, 0 ] == 1.0 && cu[ 0, 1 ] == 2.0 && cu[ 0, 2 ] == 4.0 &&
			    cu[ 1, 0 ] == 0.0 && cu[ 1, 1 ] == 3.0 && cu[ 1, 2 ] == 5.0 &&
			    cu[ 2, 0 ] == 0.0 && cu[ 2, 1 ] == 0.0 && cu[ 2, 2 ] == 12.0 );

		{
			auto r = ru.column( 1 );
			r *= 2.0;
		}
		
		assert( ru[ 0, 0 ] == 1.0 && ru[ 0, 1 ] == 4.0 && ru[ 0, 2 ] == 3.0 &&
			    ru[ 1, 0 ] == 0.0 && ru[ 1, 1 ] == 8.0 && ru[ 1, 2 ] == 5.0 &&
			    ru[ 2, 0 ] == 0.0 && ru[ 2, 1 ] == 0.0 && ru[ 2, 2 ] == 6.0 );
		
		cl[][1][] = [9.0,8.0,7.0];
		assert( cl[ 0, 0 ] == 1.0 && cl[ 0, 1 ] == 0.0 && cl[ 0, 2 ] == 0.0 &&
			    cl[ 1, 0 ] == 2.0 && cl[ 1, 1 ] == 8.0 && cl[ 1, 2 ] == 0.0 &&
			    cl[ 2, 0 ] == 3.0 && cl[ 2, 1 ] == 7.0 && cl[ 2, 2 ] == 6.0 );
		
		rl[1][][] = [9.0,8.0,7.0];
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


/************************************** Tests for Symmetric Matrices **************************************/
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
		auto a = Mat( 3, 0.0 );
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
			r *= 2.0;
		}
		
		assert( cu[ 0, 0 ] == 1.0 && cu[ 0, 1 ] == 2.0 && cu[ 0, 2 ] == 4.0 &&
			    cu[ 1, 0 ] == 2.0 && cu[ 1, 1 ] == 3.0 && cu[ 1, 2 ] == 5.0 &&
			    cu[ 2, 0 ] == 4.0 && cu[ 2, 1 ] == 5.0 && cu[ 2, 2 ] == 12.0 );

		{
			auto r = ru.column( 1 );
			r *= 2.0;
		}
		
		assert( ru[ 0, 0 ] == 1.0 && ru[ 0, 1 ] == 4.0 && ru[ 0, 2 ] == 3.0 &&
			    ru[ 1, 0 ] == 4.0 && ru[ 1, 1 ] == 8.0 && ru[ 1, 2 ] == 5.0 &&
			    ru[ 2, 0 ] == 3.0 && ru[ 2, 1 ] == 5.0 && ru[ 2, 2 ] == 6.0 );
		
		cl[][1][] = [9.0,8.0,7.0];
		assert( cl[ 0, 0 ] == 1.0 && cl[ 0, 1 ] == 2.0 && cl[ 0, 2 ] == 3.0 &&
			    cl[ 1, 0 ] == 2.0 && cl[ 1, 1 ] == 8.0 && cl[ 1, 2 ] == 7.0 &&
			    cl[ 2, 0 ] == 3.0 && cl[ 2, 1 ] == 7.0 && cl[ 2, 2 ] == 6.0 );
		
		rl[1][][] = [9.0,8.0,7.0];
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