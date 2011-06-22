module scid.vector;

import scid.internal.arraystorage;
import scid.internal.arrayviewstorage;

struct BaseVector( Storage_ ) {
	alias Storage_                  Storage;
	alias Storage.ElementType       ElementType;
	alias BaseVector!(Storage.View) View;
	
	@property
	Storage storage() { return stor_; }
	
	this( A... )( A args ) {
		stor_ = Storage( args );
	}
	
	this( this ) {
		Storage oldStor_ = stor_;
		stor_ = Storage();
		stor_.assign( oldStor_ );
	}
	
	ref BaseVector opAssign( ref BaseVector rhs ) {
		stor_.assign( rhs.stor_ );
		return this;
	}
	
	ref BaseVector opAssign( ref BaseVector!(Storage.View) rhs ) {
		stor_.assign( rhs.stor_ );
		return this;
	}
	
	ElementType opIndex( size_t i ) {
		return stor_.at( i );
	}
	
	void opIndexAssign( T )( T rhs, size_t i ) {
		stor_.set( rhs, i );
	}
	
	void opIndexOpAssign( string op, T )( T rhs, size_t i ) {
		stor_.set( mixin("stor_.at(i)" ~ op ~ "rhs"), i );
	}
	
	View opSlice() {
		return View( stor_.slice( 0, stor_.length ) );
	}
	
	View opSlice( size_t i, size_t j ) {
		return View( stor_.slice( i, j ) );
	}
	
	void opSliceAssign( T )( T rhs, size_t i, size_t j ) {
		stor_.sliceAssign( rhs, i, j );
	}
	
	void opSliceAssign( T )( T rhs ){
		stor_.sliceAssign( rhs, 0, stor_.length );
	}
	
	@property
	size_t length() { return stor_.length; }
	
	string toString() {
		return stor_.toString();
	}
	
	/** InputRange method. Return the first element in the array. */
	@property ElementType front() const { return stor_.front; }
	
	/** InputRange method. Return the last element in the array. */
	@property ElementType back() const  { return stor_.back; }
	
	/** InputRange method. Return whether the array is empty. */
	@property bool        empty() const { return stor_.empty; }
	
	/** InputRange method. Remove the first element. */
	void                  popFront()    { stor_.popFront(); }
	
	/** InputRange method. Remove the last element. */
	void                  popBack()     { stor_.popBack(); }
	
	
	private Storage stor_;
	
}

template Vector( T )     { alias BaseVector!( ArrayStorage!T )     Vector;     }
template VectorView( T ) { alias BaseVector!( ArrayViewStorage!T ) VectorView; }