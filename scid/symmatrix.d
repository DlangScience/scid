module scid.symmatrix;

import scid.matrix;
import scid.trimatrix;
import std.typecons;

import scid.internal.cowarray;
import scid.internal.vecadapters;
import scid.internal.arraydata;

import scid.common.traits;

import std.conv, std.complex, std.array, std.math;


/** Symmetric and hermitian matrix with packed storage implementation. For helpful interfacing both the stored
  * triangle (upper/lower) and the storage order (row major / column major) can be specified. The storage can
  * also be set to a custom container. If T is complex, the matrix is automatically hermitian.
  *
  * Unlike general matrices, triangular & symmetric matrices offer only 'virtual' views to their rows and
  * columns that need to be copied to Vector-s before BLAS usage.
  *
  * To use with BLAS functions (and other functions taking a builtin array), use the data() and cdata()
  * Whenever possible, use cdata() as data() might cause unnecessary memory duplication.
  */
struct SymmetricMatrix( T, MatrixTriangle triangle_ = MatrixTriangle.Upper, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray ) {
	/** The type of elements in the vector. */
	alias T                ElementType;
	
	/** The type of the Array that's wrapped. */
	alias ArrayTemplate!T    Wrapped;
	
	/** The RefCounted array type. */
	alias RefCounted!Wrapped Array;
	
	/** Is the matrix stored in row-major order? */
	enum isRowMajor  = StorageOrder.RowMajor == storageOrder_;
	
	/** Is the upper triangle stored? */
	enum isUpper     = MatrixTriangle.Upper  == triangle_;
	
	/** Is the matrix hermitian? */
	enum isHermitian = isComplex!T;

	/** Create a new matrix of a given size (NxN). */
	this( size_t size, T initializeWith = T.init ) {
		array_.RefCounted.initialize( size * ( size + 1 ) / 2, initializeWith );
		size_ = size;
	}
	
	/** Create a new matrix of a given size (NxN). Do not initialize. */
	this( size_t size, void* ) {
		array_.RefCounted.initialize( size * ( size + 1 ) / 2, null );
		size_ = size;
	}

	/** Create a matrix given a builtin array of arrays.
	  * The type of arrays is the same as in TriangularMatrix-s.
	  */
	this( T[][] matrix ) {
		this = matrix;
	}

	/** Postblit ctor calls opAssign on the wrapped array. */
	this( this ) {
		// Workaround for bug 6199.
		if( array_.RefCounted.isInitialized() )
			array_.RefCounted.initialize( *&array_.refCountedPayload() );
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
	
	/** Return the ref-counted array. */
	@property Array  array()         { return array_; }
	
	/** Return the number of rows. Defined for compatibility with generic algorithms. */
	@property size_t rows() const    { return size_; }
	
	/** Return the number of columns. Defined for compatibility with generic algorithms. */
	@property size_t columns() const { return size_; }
	
	/** Return the number of rows & columns. */
	@property size_t size() const    { return size_; }

	/** Assignment to another matrix. Forwards to the wrapped array. */
	ref typeof(this) opAssign( typeof( this ) rhs ) {
		// Workaround for bug 6199.
		array_   = *&rhs.array_.refCountedPayload();
		size_    = rhs.size_;
		return this;
	}
	
	/** Assignment to a built-in array of arrays. See TriangularMatrix. */
	ref typeof(this) opAssign( T[][] matrix ) {
		array_.RefCounted.initialize( matrix.length * ( matrix.length + 1 ) / 2 );
		size_ = matrix.length;

		foreach( i, row ; matrix ) {
			static if( isUpper ) {
				assert( row.length == size_ - i, "Invalid array for symmetric matrix." );	
				foreach( j, value ; row )
					this[ i, j + i ] = value;
			} else {
				assert( row.length == i + 1, "Invalid array for symmetric matrix." );
				foreach( j, value ; row )
					this[ i, j ] = value;
			}
		}

		return this;
	}
	
	/** Access to individual elements. */
	T opIndex( size_t row, size_t column ) const {
		return index_( row, column );
	}
	
	/// ditto
	void opIndexAssign( S )( S rhs, size_t row, size_t column ) if( is(S : T ) ) { 
		indexAssign_( rhs, row, column );
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t row, size_t column ) { indexOpAssign_!op( rhs, row, column ); }
	
	/** Access to individual rows. */
	RowView opIndex( size_t row )
	in   {
		assert( row < size_ );
	} body {
		return RowView( array_, size_, 0, size_, row );
	}

	/** Slicing ops. Only vector slices are allowed e.g A[0..2][3] or A[1][]. */
	SlicedProxy opSlice( size_t a, size_t b ) {
		return typeof( return )( &this, a, b );
	}

	/// ditto
	SlicedProxy opSlice() {
		return opSlice( 0, size_ );	
	}
	
	/// ditto
	void opSliceAssign( S )( S rhs ) if( is( S == T ) ) {
		array_.refCountedPayload()[] = rhs;
	}

	/** Returned by opSlice, provides access to individual columns */
	struct SlicedProxy {
		private SymmetricMatrix* m_;
		private size_t            start_, end_;
		
		/** Access to individual columns. */
		ColumnView opIndex( size_t j ) {
			return ColumnView( m_.array_, m_.size, start_, end_, j );
		}
		
		/// ditto
		void opIndexAssign( S )( S rhs, size_t i ) {
			opIndex( i )[] = rhs;
		}
	}
	
	/** An abstract view of a given row. Provides Vector-like functionality. */
	struct RowView {
		mixin SubVector!(VectorType.Row);
		
		/** Access to elements. */
		T opIndex( size_t i ) const                          { return index_( fixed_, start_ + i ); }
		/// ditto
		void opIndexAssign( T rhs, size_t i )                { return indexAssign_( rhs, fixed_, start_ + i ); }
		/// ditto
		void opIndexOpAssign( string op )( T rhs, size_t i ) { return indexOpAssign_!op( rhs, fixed_, start_ + i ); }
	}
	
	/** An abstract view of a given column. Provides Vector-like functionality. */
	struct ColumnView {
		mixin SubVector!(VectorType.Column);
		
		/** Access to elements. */
		T opIndex( size_t i ) const                          { return index_( start_ + i, fixed_ ); }
		
		/// ditto
		void opIndexAssign( T rhs, size_t i )                { return indexAssign_( rhs, start_ + i, fixed_ ); }
		
		/// ditto
		void opIndexOpAssign( string op )( T rhs, size_t i ) { return indexOpAssign_!op( rhs, start_ + i, fixed_ ); }
	}
	
	static if( isHermitian ) {
		/** Generic conjugate - works on both cdouble & Complex!double. */
		static const (C) genConj( C )( C z ) {
			static if( __traits( compiles, z.conj ) )
				return z.conj;
			else
				return conj(z);
			
		}
	}
	
	private mixin template IndexMapper() {
		private size_t mapHelper_( bool colUpper )( size_t i, size_t j ) const {
			static if( colUpper ) return i + j * (j + 1) / 2;
			else                  return i + ( ( size_ + size_ - j - 1 ) * j ) / 2;
		}
		
		private size_t map_( size_t i, size_t j ) const {
			static if( isRowMajor )
				return mapHelper_!( !isUpper )( j, i );
			else
				return mapHelper_!( isUpper )( i, j );
		}

		private bool needSwap_( size_t i, size_t j ) const {
			static if( isUpper ) {
				return i > j;
			} else {
				return i < j;
			}
		}
		
		private T index_( size_t row, size_t column ) const
		in {
			assert( row < size_ );
			assert( column < size_ );
		} body {
			if( needSwap_( row, column ) ) {
				static if( isHermitian ) {
					return genConj(array_.refCountedPayload()[ map_( column, row ) ]);
				} else {
					return array_.refCountedPayload()[ map_( column, row ) ];
				}
			} else {
				return array_.refCountedPayload()[ map_( row, column ) ];
			}
		}	

		private void indexAssign_( T rhs, size_t row, size_t column )
		in {
			assert( row < size_ );
			assert( column < size_ );
		} body {
			if( needSwap_( row, column ) ) {
				static if( isHermitian ) {
					array_.refCountedPayload()[ map_( column, row ) ] = genConj(rhs);
				} else {
					array_.refCountedPayload()[ map_( column, row ) ] = rhs;
				}
			} else {
				array_.refCountedPayload()[ map_( row, column ) ] = rhs;
			}
		}

		private void indexOpAssign_( string op )( T rhs, size_t row, size_t column )
		in {
			assert( row < size_ );
			assert( column < size_ );
		} body {
			if( needSwap_( row, column ) ) {
				static if( isHermitian ) {
					mixin( "array_.refCountedPayload()[ map_( column, row ) ]" ~ op ~ "= genConj(rhs);" );
				} else {
					mixin( "array_.refCountedPayload()[ map_( column, row ) ]" ~ op ~ "= rhs;" );
				}
			} else {
				mixin( "array_.refCountedPayload()[ map_( row, column ) ]" ~ op ~ "= rhs;" );
			}
		
		}
	}
	
	private mixin template SubVector( VectorType vtype ) {
		mixin IndexMapper;
		mixin MatrixVectorView!( T, Array, vtype );
		
		
		private this( Array array, size_t size, size_t start, size_t end, size_t fixed ) {
			array_ = array;
			size_  = size;
			setView_( start, end, fixed );
		}
		
		/** Assignment to a different sub vector. */
		ref typeof( this ) opAssign( typeof( this ) rhs ) {
			array_ = rhs.array_;
			size_  = rhs.size_;
			setView_( rhs );
			return this;
		}
			
		private size_t size_;
	}
	
	// Provide toString()
	mixin MatrixToString;
	
	// The index mapping is wrapped in a mixin.
	mixin IndexMapper;
	
	private size_t size_;
	private Array  array_;
}

import std.stdio;

unittest {
	void testNonHerm( StorageOrder order )() {
		alias SymmetricMatrix!(int, MatrixTriangle.Upper, order ) UpMat;
		alias SymmetricMatrix!(int, MatrixTriangle.Lower, order ) LoMat;
		
		auto x = UpMat( 3 );
		auto y = LoMat( 3 );
		
		assert( x.size == 3, "Upper symmetric array ctor size failed." );
		assert( y.size == 3, "Lower symmetric array ctor size failed." );
		
		x[ 0,1 ] = 42;
		y[ 1,0 ] = 42;
		assert( x[0, 1] == 42 && y[1, 0] == 42, "Set failed" );
		assert( x[1, 0] == 42 && y[1, 0] == 42, "Symmetry failed");
		
		auto z = x;
		z[ 0, 1 ] = 1;
		assert( z[ 0, 1 ] == 1 && x[ 0, 1 ] == 42, "Copy duplication failed." );
		assert( z[ 1, 0 ] == 1, "Copy symmetry failed." );
		
		
		auto a = UpMat([ [1,2,3], [4,5], [6] ]);
		auto b = LoMat([ [1], [2,3], [4,5,6] ]);
		
		assert( a.size == 3, "Upper symmetric array ctor size failed." );
		assert( b.size == 3, "Lower symmetric array ctor size failed." );
		
		assert( a[0, 0] == 1 && a[0, 1] == 2 && a[0, 2] == 3 &&
			    a[1, 0] == 2 && a[1, 1] == 4 && a[1, 2] == 5 &&
			    a[2, 0] == 3 && a[2, 1] == 5 && a[2, 2] == 6,
			    "Upper symmetric array ctor equality failed" );
		
		assert( b[0, 0] == 1 && b[0, 1] == 2 && b[0, 2] == 4 &&
			    b[1, 0] == 2 && b[1, 1] == 3 && b[1, 2] == 5 &&
			    b[2, 0] == 4 && b[2, 1] == 5 && b[2, 2] == 6,
			   "Lower symmetric array ctor equality failed" );
		
		// Check row/column views
		assert( to!string( a[][ 0 ] ) == "[1, 2, 3]" &&
		        to!string( a[][ 1 ] ) == "[2, 4, 5]" &&
		        to!string( a[][ 2 ] ) == "[3, 5, 6]",
			   "Column extraction failed on upper symmetric matrix" );
		
		assert( to!string( b[][ 0 ] ) == "[1, 2, 4]" &&
		        to!string( b[][ 1 ] ) == "[2, 3, 5]" &&
		        to!string( b[][ 2 ] ) == "[4, 5, 6]",
			   "Column extraction failed on lower symmetric matrix" );
		
		assert( to!string( a[ 0 ][] ) == "[1, 2, 3]" &&
		        to!string( a[ 1 ][] ) == "[2, 4, 5]" &&
		        to!string( a[ 2 ][] ) == "[3, 5, 6]",
			   "Row extraction failed on upper symmetric matrix" );
		
		assert( to!string( b[ 0 ][] ) == "[1, 2, 4]" &&
		        to!string( b[ 1 ][] ) == "[2, 3, 5]" &&
		        to!string( b[ 2 ][] ) == "[4, 5, 6]",
			   "Row extraction failed on lower symmetric matrix" );
		
	}
	
	void testHerm( StorageOrder order )() {
		alias Complex!double Z;
		alias SymmetricMatrix!(Z, MatrixTriangle.Upper, order ) UpMat;
		alias SymmetricMatrix!(Z, MatrixTriangle.Lower, order ) LoMat;
		
		UpMat x = UpMat( 3 );
		LoMat y = LoMat( 3 );
		
		x[0, 1] = Z(42,42);
		auto z = x;
		z[ 0, 1 ] = Z(1,1);
		assert( z[ 0, 1 ] == Z(1,1) && x[ 0, 1 ] == Z(42,42), "Copy duplication failed." );
		assert( z[ 1, 0 ] == Z(1,-1), "Copy symmetry failed." );
		
		auto k = Z(0.0,0.0);
		foreach( i ; 0 .. 3 )
			foreach( j ; 0 .. (i+1) ) {
				k += Z(1.0,1.0);
				x[i, j]     = k;
				y[2-i, 2-j] = k;
			}
	
		k = Z(0.0,0.0);
		foreach( i ; 0 .. 3 )
			foreach( j ; 0 .. (i+1) ) {
				k += Z(1.0,1.0);
				assert( x[j, i] == (i==j ? k : UpMat.genConj(k) )     && x[i,j] == k,     "Upper Hermitian set symmetry failed." );
				assert( y[2-j, 2-i] == (i==j ? k : UpMat.genConj(k) ) && y[2-i,2-j] == k, "Lower Hermitian set symmetry failed." );
			}
		
		// Check row/column views
		assert( to!string( x[][ 0 ] ) == "[1+1i, 2+2i, 4+4i]" &&
		        to!string( x[][ 1 ] ) == "[2-2i, 3+3i, 5+5i]" &&
		        to!string( x[][ 2 ] ) == "[4-4i, 5-5i, 6+6i]",
			   "Column extraction failed on upper hermitian matrix" );
		
		
		assert( to!string( y[][ 0 ] ) == "[6+6i, 5-5i, 4-4i]" &&
		        to!string( y[][ 1 ] ) == "[5+5i, 3+3i, 2-2i]" &&
		        to!string( y[][ 2 ] ) == "[4+4i, 2+2i, 1+1i]",
			   "Column extraction failed on lower hermitian matrix" );		
		
		assert( to!string( x[ 0 ][] ) == "[1+1i, 2-2i, 4-4i]" &&
		        to!string( x[ 1 ][] ) == "[2+2i, 3+3i, 5-5i]" &&
		        to!string( x[ 2 ][] ) == "[4+4i, 5+5i, 6+6i]",
			   "Row extraction failed on upper hermitian matrix" );
		
		assert( to!string( y[ 0 ][] ) == "[6+6i, 5+5i, 4+4i]" &&
		        to!string( y[ 1 ][] ) == "[5-5i, 3+3i, 2+2i]" &&
		        to!string( y[ 2 ][] ) == "[4-4i, 2-2i, 1+1i]",
			   "Row extraction failed on lower hermitian matrix" );
	}
	
	testNonHerm!( StorageOrder.ColumnMajor )();
	testNonHerm!( StorageOrder.RowMajor    )();
	testHerm!   ( StorageOrder.ColumnMajor )();
	testHerm!   ( StorageOrder.RowMajor    )();
	
}