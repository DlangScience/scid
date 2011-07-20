module scid.storage.array;

import scid.internal.assertmessages;
import scid.internal.expression;
import scid.internal.hlblas;
import scid.storage.cowarray;
import scid.storage.arrayview;
import scid.common.traits;
import std.traits, std.range, std.algorithm;
import scid.vector;

template ArrayStorage( ElementOrArray, VectorType vectorType = VectorType.Column )
		if( isFortranType!(BaseElementType!ElementOrArray) ) {
	
	static if( isFortranType!ElementOrArray )
		alias BasicArrayStorage!( CowArrayRef!ElementOrArray, vectorType ) ArrayStorage;
	else
		alias BasicArrayStorage!( ElementOrArray, vectorType )             ArrayStorage;
}

struct BasicArrayStorage( ArrayRef_, VectorType vectorType_ ) {
	alias ArrayRef_                                                       ArrayRef;
	alias BaseElementType!ArrayRef                                        ElementType;	
	alias vectorType_                                                     vectorType;
	alias ArrayViewStorage!( ArrayRef, vectorType )                       View;
	alias StridedArrayViewStorage!( ArrayRef, vectorType )                StridedView;
	alias BasicArrayStorage!( ArrayRef_, transposeVectorType!vectorType ) Transposed;
	alias array_                                                          this;
	
	enum  stride     = 1;
	enum  firstIndex = 0;
	
	this( size_t newLength ) {
		array_.RefCounted.ensureInitialized();
		array_ = ArrayRef( newLength );
	}
	
	this( size_t newLength, void* ) {
		array_.RefCounted.ensureInitialized();
		array_ = ArrayRef( newLength, null );
	}
	
	this( ElementType[] initializer ) {
		array_.RefCounted.ensureInitialized();
		array_ = ArrayRef( initializer );
	}
	
	this( ArrayRef arrayRef ) {
		array_.RefCounted.ensureInitialized();
		array_ = arrayRef;
	}
	
	this( this ) {
		array_.RefCounted.ensureInitialized();
		array_ = ArrayRef( array_.ptr );
	}
	
	void forceRefAssign( ref typeof(this) rhs ) {
		array_.RefCounted.ensureInitialized();
		array_ = rhs.array;	
	}
		
	ref typeof(this) opAssign( typeof(this) rhs ) {
		array_.RefCounted.ensureInitialized();
		array_ = ArrayRef( rhs.array_.ptr );
		return this;
	}
	
	typeof( this ) slice( size_t start, size_t end ) {
		return typeof(this)( ArrayRef( array_.slice( start, end ).ptr )  );
	}
	
	View view( size_t start, size_t end ) {
		return typeof( return )( array_, start, end - start );
	}
	
	StridedView view( size_t start, size_t end, size_t stride ) {
		size_t len = (end - start);
		len = len / stride + (len % stride != 0);
		return typeof( return )( array_, start, len, stride );
	}
	
	@property {
		ArrayRef     array()       { return array_; }
	}
	
	// Workaround: alias this doesn't make this a range //
	@property {
		bool                empty()  const { return array_.empty; }
		ElementType         front()  const { return array_.front; }
		ElementType         back()   const { return array_.back; }
		size_t              length() const { return array_.RefCounted.isInitialized() ? array_.length : 0; }
		const(ElementType)* cdata()  const { return array_.RefCounted.isInitialized() ? array_.cdata : null; }
		ElementType*        data()         { return array_.RefCounted.isInitialized() ? array_.data : null; }
	}
	
	void popFront() { array_.popFront(); }
	void popBack()  { array_.popBack(); }
	// Workaround ends. //
	
	// vector ops
	void resizeOrClear( size_t rlength ) {
		if( array_.RefCounted.isInitialized() ) 
			array_.resizeOrClear( rlength );
		else
			array_ = ArrayRef( rlength );
	}
	
	void resizeOrClear( size_t rlength, void* ) {
		if( array_.RefCounted.isInitialized() ) 
			array_.resizeOrClear( rlength, null );
		else
			array_ = ArrayRef( rlength, null );
	}
	
	void copy( Transpose tr, Source )( ref Source rhs ) if( isStridedStorage!(Source,ElementType) ) {
		static if( is( Source : ArrayStorage!( ArrayRef, transposeVectorType!(vectorType,tr) ) ) ) {
			array_.RefCounted.ensureInitialized();
			this.array_ = ArrayRef( rhs.array_.ptr );	
		} else static if( is( Source : ArrayViewStorage!( ArrayRef, transposeVectorType!(vectorType,tr) ) ) ) {
			array_.RefCounted.ensureInitialized();
			this.array_ = rhs.array.slice( rhs.firstIndex, rhs.firstIndex + rhs.length );	
		} else {
			
			hlStridedCopy( rhs, this );
		}
	}
	
	mixin StridedAxpyScalDot;
	
	ArrayRef array_;
private:
	mixin SliceIndexMessages;

	
}
unittest {
	// TODO: Write tests for ArrayStorage.
}

