module scid.vector;

import scid.storage.array;
import scid.storage.arrayview;
import scid.common.traits;

import std.traits, std.range, std.algorithm, std.conv;

enum VectorType {
	Row, Column
}

template Vector( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrStorage) ) {
	
	static if( isFortranType!ElementOrStorage )
		alias BasicVector!( ArrayStorage!ElementOrStorage, vectorType ) Vector;
	else
		alias BasicVector!( ElementOrStorage, vectorType )              Vector;
}

template VectorView( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrStorage) ) {
	
	alias BasicVector!( ArrayViewStorage!ElementOrStorage, vectorType) VectorView;
}

template StridedVectorView( ElementOrStorage, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrStorage) ) {
	
	alias BasicVector!( StridedArrayViewStorage!ElementOrStorage, vectorType).View StridedVectorView;
}

template isVector( T ) {
	static if( is( typeof( T.init[0]          ) ) &&
			   is( typeof( T.init[0..1]       ) ) &&
			   is( typeof( T.init.storage     ) ) &&
			   is( typeof( T.init.view(0,0)   ) ) &&
			   is( typeof( T.init.view(0,0,1) ) ) &&
			   isInputRange!T )
		enum isVector = true;
	else
		enum isVector = false;
}

template signOfOp( string op, T ) {
	static if( op == "+" )
		enum T signOfOp = 1;
	else static if( op == "-" )
		enum T signOfOp = -1;
}

struct BasicVector( Storage_, VectorType vectorType_ = VectorType.Column ) {
	alias BaseElementType!Storage_                        ElementType;
	alias Storage_                                        Storage;
	alias vectorType_                                     vectorType;
	alias BasicVector!( Storage.View, vectorType )        View;
	alias BasicVector!( Storage.StridedView, vectorType ) StridedView;
	
	template isCompatible( V ) {
		static if( is( V : ElementType[] ) || (isVector!V && is(V.ElementType == ElementType)) )
			enum isCompatible = true;
		else
			enum isCompatible = false;
	}
	
	this( A... )( A args ) if( A.length > 0 && !is( A[0] : Storage ) && !isVector!(A[0]) ) {
		stor_ = Storage(args);
		//emplace!Storage(&stor_, args );
	}
	
	this( A )( BasicVector!(A, vectorType) other ) {
		static if( is( A : Storage ) ) stor_ = other.stor_;
		else                           stor_.copy( other.stor_ );
	}
	
	this()( Storage stor ) {
		stor_ = stor;
		//move( stor, stor_ );
	}
	
	ElementType opIndex( size_t i ) const {
		return stor_.index( i );
	}
	
	void opIndexAssign( ElementType rhs, size_t i ) {
		stor_.indexAssign( rhs, i );
	}
	
	void opIndexOpAssign( string op )( ElementType rhs, size_t i ) {
		stor_.indexAssign!op( rhs, i );
	}
	
	ref typeof(this) opAssign( typeof(this) rhs ) {
		move( rhs.stor_, stor_ );
		return this;
	}
	
	ref typeof(this) opAssign( ElementType[] rhs ) {
		stor_.copyLeft( rhs );
		return this;
	}
	
	ref typeof(this) opOpAssign( string op, Rhs )( Rhs rhs ) if( (op == "+" || op == "-") && isCompatible!Rhs ) {
		stor_.axpyLeft( getStorage(rhs), signOfOp!(op,ElementType) );
		return this;
	}
	
	ref typeof(this) opOpAssign( string op )( ElementType rhs ) if(op == "*" || op == "/" ) {
		static if( op == "/" ) rhs = 1/rhs;
		stor_.scalLeft( rhs );
		return this;
	}
	
	typeof( this ) opSlice() {
		return typeof(this)( stor_ );
	}
	
	typeof( this ) opSlice( size_t start, size_t end ) {
		return typeof(this)( stor_.slice( start, end ) );	
	}

	void opSliceAssign( Rhs )( Rhs rhs ) if( isCompatible!Rhs ) {
		stor_.copyLeft( getStorage(rhs) );
	}
	
	void opSliceAssign( Rhs )( Rhs rhs, size_t start, size_t end ) if( isCompatible!Rhs ) {
		stor_.slicedCopyLeft( getStorage(rhs), start, end );
	}
	
	void opSliceOpAssign( string op, Rhs )( Rhs rhs, size_t start, size_t end ) if( (op == "+" || op == "-") && isCompatible!Rhs ) {
		stor_.slicedAxpyLeft( getStorage(rhs), signOfOp!(op,ElementType), start, end );
	}
	
	void opSliceOpAssign( string op )( ElementType rhs, size_t start, size_t end ) if(op == "*" || op == "/" ) {
		static if( op == "/" ) rhs = 1/rhs;
		stor_.slicedScalLeft( rhs, start, end );
	}
	
	View view( size_t start, size_t end ) {
		return typeof( return )( stor_.view( start, end ) );
	}
	
	StridedView view( size_t start, size_t end, size_t stride ) {
		return typeof( return )( stor_.view( start, end, stride ) );	
	}
	
	void popFront() { stor_.popFront(); }
	void popBack()  { stor_.popBack(); }
	
	@property {
		ref Storage          storage()       { return stor_; }
		ref const(Storage)   storage() const { return stor_; }
		
		ElementType[]        data()          { return stor_.data; }
		const(ElementType[]) cdata()   const { return stor_.cdata; }
		
		bool                 empty()   const { return stor_.empty; }
		ElementType          front()   const { return stor_.front; }
		ElementType          back()    const { return stor_.back; }
		size_t               length()  const { return stor_.length; }
	}
	
private:
	ref auto getStorage( V )( ref V rhs ) if( isCompatible!V ) {
		static if( is( V : ElementType[] ) )
			return rhs;
		else
			return rhs.storage;
	}

	Storage stor_;
}

unittest {
	// TODO: Write tests for Vector.
}