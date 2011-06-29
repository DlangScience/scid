module scid.internal.vecadapters;

public import std.typecons : RefCounted;
public import scid.vector  : VectorType;

mixin template MatrixVectorView( T, ArrayT_, VectorType vtype ) {
	alias T ElementType;
	enum isRowVector = (vtype == VectorType.Row);
	
	typeof( this ) opSlice() {
		return opSlice( 0, length );
	}
	
	typeof( this ) opSlice( size_t newStart, size_t newEnd ) {
		auto ret = this;
		ret.start_ = start_ + newStart;
		ret.end_   = start_ + newEnd;
		
		return ret;
	}
	
	typeof( this ) opSliceAssign( S )( S rhs ) {
		opSliceAssign( rhs, 0, length );
	}
	
	void opSliceAssign( S )( S rhs, size_t a, size_t b ) if( is( S == T ) ) {
		foreach( i ; a .. b )
			this[ i ] = rhs;
	}
	
	void opSliceAssign( S )( S rhs, size_t a, size_t b ) if( isInputRange!S ) {
		static if( __traits( compiles, rhs.length ) ) {
			assert( rhs.length == b - a );	
		}
		
		foreach( i ; a .. b ) {
			assert( !S.empty(), "Range underrun in slice assignment." );
			this[ i ] = S.front();
			S.popFront();
		}
	}
	
	@property bool		  empty()  const { return end_ == start_; }
	@property size_t	  length() const { return end_ -  start_; }
	@property ElementType front()  const { return this[ 0 ]; }
	@property ElementType back()   const { return this[ length - 1 ]; }
	
	void popFront()
	in   { assert( !empty(), "popFront() on empty " ~ typeof( this ).stringof ); }
	body { ++ start_; }
	
	void popBack()
	in   { assert( !empty(), "popBack() on empty " ~ typeof( this ).stringof ); }
	body { -- end_; }
	
	private size_t premap_( size_t f, size_t v ) const {
		static if( isRowVector ) return map_( f, v );
		else                     return map_( v, f );
	}
	
	private void setView_( ref typeof( this ) rhs ) {
		start_ = rhs.start_;
		end_   = rhs.end_;
		fixed_ = rhs.fixed_;
	}
	
	private void setView_( size_t s, size_t e, size_t f ) {
		start_ = s;
		end_   = e;
		fixed_ = f;
	}
	
	private ArrayT_ array_;
	private size_t  start_, end_, fixed_;
}

