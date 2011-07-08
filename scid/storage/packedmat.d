module scid.storage.packedmat;

import scid.storage.packedsubmat;
import scid.storage.packedsubvec;
import scid.internal.assertmessages;
import scid.common.traits;
import scid.matrix;
import scid.vector;
import std.algorithm;

struct PackedStorage( MatrixRef_ ) {
	alias MatrixRef_                                                         MatrixRef;
	alias BaseElementType!MatrixRef                                          ElementType;
	alias storageOrderOf!MatrixRef                                           storageOrder;
	alias PackedSubMatrixStorage!( MatrixRef, SubMatrixType.Slice)           Slice;
	alias PackedSubMatrixStorage!( MatrixRef, SubMatrixType.View)            View;
	alias Vector!( PackedSubVectorStorage!( MatrixRef, VectorType.Row ) )    RowView;
	alias Vector!( PackedSubVectorStorage!( MatrixRef, VectorType.Column ) ) ColumnView;
	alias ColumnView                                                         DiagonalView;
	alias matrix_                                                            this;
	
	enum isRowMajor   = (storageOrder == StorageOrder.RowMajor );
	
	this(A...)( A args ) {
		matrix_ = MatrixRef( args );
	}
	
	this( this ) {
		matrix_ = MatrixRef( matrix_.ptr );	
	}
	
	void forceRefAssign( ref typeof(this) rhs ) {
		matrix_ = rhs.matrix_;
	}
		
	ref typeof(this) opAssign( typeof(this) rhs ) {
		move( rhs.matrix_, matrix_ );
		return this;
	}
	
	RowView row( size_t i )
	in {
		assert( i < this.rows, sliceMsg_( i,0, i,this.columns ) );
	} body {
		return RowView( matrix_, i, 0, this.columns );
	}
	
	ColumnView column( size_t j )
	in {
		assert( j < this.columns, sliceMsg_( 0,j, this.rows,j ) );
	} body {
		return ColumnView( matrix_, j, 0, this.rows );
	}
	
	RowView rowSlice( size_t i, size_t start, size_t end )
	in {
		assert( i < rows && start < end && end <= columns, sliceMsg_( i, start, i, end ) );
	} body {
		return RowView( matrix_, i, start, end - start );
	}
	
	ColumnView columnSlice( size_t j, size_t start, size_t end )
	in {
		assert( j < columns && start < end && end <= columns, sliceMsg_( start,j, end,j ) );
	} body {
		return ColumnView( matrix_, j, start, end - start );
	}
	
	View view( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( matrix_, rowStart, rowEnd - rowStart, colStart, colEnd - colStart );
	}
	
	Slice slice( size_t rowStart, size_t rowEnd, size_t colStart, size_t colEnd ) {
		return typeof( return )( matrix_, rowStart, rowEnd - rowStart, colStart, colEnd - colStart );
	}
	
	@property {
		auto front() {
			static if( isRowMajor ) return row( 0 );
			else                    return column( 0 );
		}
		
		auto back() {
			static if( isRowMajor ) return row( this.rows - 1 );
			else                    return column( this.columns - 1 );
		}
		
		// workaround for alias this problems {
		auto size()    const { return matrix_.size; }
		
		alias size rows;
		alias size columns;
		// }
	}
	
	MatrixRef matrix_;
	
private:
	mixin SliceIndex2dMessages;
}
