module scid.trimatrix;

import scid.matrix;
import scid.internal.cowarray;
import scid.internal.vecadapters;
import std.array;
import std.conv;
import std.typecons;

enum MatrixTriangle {
	Upper, Lower	
}

/** Triangular matrix implementation with packed storage. Can be parametrized with the type of triangle it represents (upper/lower),
  * storage order (row-major/column-major) and the container used for storage.
  *
  * Unlike general matrices, triangular & symmetric matrices offer only 'virtual' views to their rows and
  * columns that need to be copied to Vector-s before BLAS usage.
  *
  * To use with BLAS functions (and other functions taking a builtin array), use the data() and cdata()
  * Whenever possible, use cdata() as data() might cause unnecessary memory duplication.
  */
struct TriangularMatrix( T, MatrixTriangle triangle_, StorageOrder storageOrder_ = StorageOrder.ColumnMajor, alias ArrayTemplate = CowArray ) {
	/** The type of elements in the vector. */
	alias T                ElementType;
	
	/** The type of the Array that's wrapped. */
	alias ArrayTemplate!T    Wrapped;
	
	/** The RefCounted array type. */
	alias RefCounted!Wrapped Array;
	
	/** Is the matrix stored in row-major order? */
	enum isRowMajor = StorageOrder.RowMajor == storageOrder_;
	
	/** Is the upper triangle stored? */
	enum isUpper    = MatrixTriangle.Upper  == triangle_;
	
	/** Create a new matrix of a given size (NxN). */
	this( size_t size ) {
		array_.RefCounted.initialize(  size * ( size + 1 ) / 2  );
		size_ = size;
	}

	/** Create a matrix given a builtin array of arrays.
	  * The sizes need to be in decreasing order for upper triangular and increasing for lower.
	  */
	this( T[][] matrix ) {
		this = matrix;
	}

	//** Postblit ctor calls opAssign on the wrapped array. */
	this( this ) {
		// Workaround for bug 6199.
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
	
	/** Return the number of rows. Defined for compatibility with generic algorithms. */
	@property size_t columns() const { return size_; }
	
	/** Return the number of rows & columns. */
	@property size_t size() const    { return size_; }
	
	
	/** Assignment to another matrix. Forwards to the wrapped array. */
	ref typeof(this) opAssign( typeof( this ) rhs ) {
		array_   = *&rhs.array_.refCountedPayload();
		size_    = rhs.size_;
		
		return this;	
	}
	
	/** Assignment to a built-in array of arrays. See constructor. */
	ref typeof(this) opAssign( T[][] matrix ) {
		array_.RefCounted.initialize( matrix.length * ( matrix.length + 1 ) / 2 );
		size_ = matrix.length;
		
		foreach( i, row ; matrix ) {
			static if( isUpper ) {
				assert( row.length == size_ - i, "Invalid array for triangular matrix." );	
				foreach( j, value ; row )
					this[ i, j + i ] = value;
			} else {
				assert( row.length == i + 1, "Invalid array for triangular matrix." );
				foreach( j, value ; row )
					this[ i, j ] = value;
			}
		}

		return this;
	}
	
	/** Element access. */
	T opIndex( size_t row, size_t column ) const
	in {
		assert( row < size_ );
		assert( column < size_ );
	} body {
		static if( isUpper ) {
			if( row > column )
				return 0;
		} else {
			if( row < column )
				return 0;
		}

		return array_.refCountedPayload()[ map_( row, column ) ];
	}
	
	/// ditto
	void opIndexAssign( T rhs, size_t row, size_t column )
	in {
		assert( row < size_ );
		assert( column < size_ );
	} body {
		static if( isUpper )
			assert( row <= column, "Modification of zero element in triangle matrix." );
		else
			assert( row >= column, "Modification of zero element in triangle matrix." );

		array_.refCountedPayload()[ map_( row, column ) ] = rhs;
	}
	
	/// ditto
	void opIndexOpAssign( string op )( T rhs, size_t row, size_t column )
	in {
		assert( row < size_ );
		assert( column < size_ );
	} body {
		static if( isUpper )
			assert( row <= column, "Modification of zero element in triangle matrix." );
		else
			assert( row >= column, "Modification of zero element in triangle matrix." );

		mixin( "array_.refCountedPayload()[ map_( row, column ) ]" ~ op ~ "= rhs;" );
	}
	
	/** Access to individual rows. */
	RowView opIndex( size_t row )
	in   { assert( row < size_ ); }
	body { return RowView( array_, size_, 0, size_, row ); }
	
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
		private TriangularMatrix* m_;
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
		
		private bool bypass_( size_t i ) const {
			static if( isUpper ) {
				if( fixed_ > i ) return true;
			} else {
				if( fixed_ < i ) return true;
			}
			
			return false;
		}
		
	}
	
	/** An abstract view of a given column. Provides Vector-like functionality. */
	struct ColumnView {
		mixin SubVector!(VectorType.Column);
		
		private bool bypass_( size_t i ) const {
			static if( isUpper ) {
				if( i > fixed_ ) return true;
			} else {
				if( i < fixed_ ) return true;
			}
			
			return false;
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
		
		/** Element access. */
		T opIndex( size_t i ) const {
			i += start_;
			
			if( bypass_( i ) )
				return 0;
			
			return array_.refCountedPayload()[ premap_( fixed_, i ) ];
		}
		
		/// ditto
		void opIndexAssign( T rhs, size_t i ) {
			i += start_;
			assert( !bypass_( i ), "Trying to modify read-only matrix element." );
			array_.refCountedPayload()[ premap_( fixed_, i ) ] = rhs;
		}
		
		/// ditto
		void opIndexOpAssign( string op )( T rhs, size_t i ) {
			i += start_;
			assert( !bypass_( i ), "Trying to modify read-only matrix element." );
			mixin("array_.refCountedPayload()[ premap_( fixed_, i ) ] " ~ op ~ "= rhs;");
		}
			
		private size_t size_;
	}
	
	mixin IndexMapper;
	mixin MatrixToString;
	
	private size_t size_;
	private Array  array_;
}


unittest {
	
	void testMat( LoMat, UpMat )() {
		auto x = UpMat( 3 );
		auto y = LoMat( 3 );
		
		assert( x.size == 3, "Upper triangle array ctor size failed." );
		assert( y.size == 3, "Lower triangle array ctor size failed." );
		
		
		x[ 0,1 ] = 42;
		y[ 1,0 ] = 42;
		assert( x[0, 1] == 42 && y[1, 0] == 42, "Set failed" );
		
		auto z = x;
		z[ 0, 1 ] = 1;
		assert( z[ 0, 1 ] == 1 && x[ 0, 1 ] == 42, "Copy duplication failed." );
		
		
		auto a = UpMat([ [1,2,3], [4,5], [6] ]);
		auto b = LoMat([ [1], [2,3], [4,5,6] ]);
		
		assert( a.size == 3, "Upper triangle array ctor size failed." );
		assert( b.size == 3, "Lower triangle array ctor size failed." );

		assert( a[0, 0] == 1 && a[0, 1] == 2 && a[0, 2] == 3 &&
			    a[1, 0] == 0 && a[1, 1] == 4 && a[1, 2] == 5 &&
			    a[2, 0] == 0 && a[2, 1] == 0 && a[2, 2] == 6,
			    "Upper triangle array ctor equality failed" );
		
		assert( b[0, 0] == 1 && b[0, 1] == 0 && b[0, 2] == 0 &&
			    b[1, 0] == 2 && b[1, 1] == 3 && b[1, 2] == 0 &&
			    b[2, 0] == 4 && b[2, 1] == 5 && b[2, 2] == 6,
			   "Lower triangle array ctor equality failed" );
		
		// Check row/column views
		assert( to!string( a[][ 0 ] ) == "[1, 0, 0]" &&
		        to!string( a[][ 1 ] ) == "[2, 4, 0]" &&
		        to!string( a[][ 2 ] ) == "[3, 5, 6]",
			   "Column extraction failed on upper triangular matrix" );
		
		assert( to!string( b[][ 0 ] ) == "[1, 2, 4]" &&
		        to!string( b[][ 1 ] ) == "[0, 3, 5]" &&
		        to!string( b[][ 2 ] ) == "[0, 0, 6]",
			   "Column extraction failed on lower triangular matrix" );
		
		assert( to!string( a[ 0 ][] ) == "[1, 2, 3]" &&
		        to!string( a[ 1 ][] ) == "[0, 4, 5]" &&
		        to!string( a[ 2 ][] ) == "[0, 0, 6]",
			   "Row extraction failed on upper triangular matrix" );
		
		assert( to!string( b[ 0 ][] ) == "[1, 0, 0]" &&
		        to!string( b[ 1 ][] ) == "[2, 3, 0]" &&
		        to!string( b[ 2 ][] ) == "[4, 5, 6]",
			   "Row extraction failed on lower triangular matrix" );
		
		
	}
	
	testMat!( TriangularMatrix!(int,MatrixTriangle.Lower), TriangularMatrix!(int,MatrixTriangle.Upper) )();
	testMat!( TriangularMatrix!(int,MatrixTriangle.Lower, StorageOrder.RowMajor), TriangularMatrix!(int,MatrixTriangle.Upper, StorageOrder.RowMajor) )();
}