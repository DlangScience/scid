module scid.common.storagetraits;

public import scid.common.traits;
import scid.storages, std.range;
import scid.matvec;

template isBasicContainer( T ) {
	enum isBasicContainer =
		isFortranType!(BaseElementType!T) &&
		is( typeof(T.init.data) == BaseElementType!T* ) &&
		is( typeof(T.init.cdata) == const(BaseElementType!T)* ) &&
		is( typeof(T.init.ptr) ) &&
		is( typeof((){ T container; container = container; }) );
}

template isBasicStorage( T ) {
	static if( !is(T.ContainerRef ) )
		enum isBasicStorage = false;
	else
		enum isBasicStorage =
			isBasicContainer!(T.ContainerRef) &&
			is( BaseElementType!(T.ContainerRef) : BaseElementType!T );
}

template hasArrayStorageIndexing( T ) {
	enum hasArrayStorageIndexing =
		is( typeof( T.init.index(0) ) : BaseElementType!T ) &&
		is( typeof( T.init.indexAssign( BaseElementType!T.init, 0 ) )  ) &&
		is( typeof( T.init.indexAssign!"+"( BaseElementType!T.init, 0 ) )  ) &&
		is( typeof( T.init.indexAssign!"*"( BaseElementType!T.init, 0 ) )  );
}

template hasMatrixStorageIndexing( T ) {
	enum hasMatrixStorageIndexing =
		is( typeof( T.init.index(0,0) ) : BaseElementType!T ) &&
		is( typeof( T.init.indexAssign( BaseElementType!T.init, 0, 0 ) )  ) &&
		is( typeof( T.init.indexAssign!"+"( BaseElementType!T.init, 0, 0 ) )  ) &&
		is( typeof( T.init.indexAssign!"*"( BaseElementType!T.init, 0, 0 ) )  );
}

template hasMatrixStorageResizing( T ) {
	enum hasMatrixStorageResizing =
		is( typeof(T.init.resize( 10, 10 )) ) &&
		is( typeof(T.init.resize( 10, 10, null )) );
}

template hasArrayStorageResizing( T ) {
	enum hasArrayStorageResizing =
		is( typeof(T.init.resize( 10 )) ) &&
		is( typeof(T.init.resize( 10, null )) );
}

template hasMatrixDimensions( T ) {
	enum hasMatrixDimensions =
		is( typeof(T.init.rows) ) &&
		is( typeof(T.init.columns) ) &&
		is( typeof(T.init.minor) ) &&
		is( typeof(T.init.major) );
}

template isVectorStorage( T ) {
	enum isVectorStorage =
		isBasicStorage!T &&
		hasArrayStorageIndexing!T &&
		hasArrayStorageResizing!T &&
		is( typeof((){
			auto s = T.init.slice( 0, 10 );
			auto v = T.init.view( 0, 10 );
			auto w = T.init.view( 0, 10, 1 );
		}()) );
}

template isMatrixStorage( T ) {
	enum isMatrixStorage =
		//isBasicStorage!T &&
		hasMatrixStorageIndexing!T &&
		hasMatrixStorageResizing!T &&
		hasMatrixDimensions!T &&
		is( T.Transposed ) &&
		is( typeof( storageOrderOf!T ) ) &&
		is( typeof((){
			auto s = T.init.slice(0,0,0,0);
			auto v = T.init.view(0,0,0,0);
		}) );
}

template isArrayContainer( T ) {
	enum isArrayContainer =
		isBasicContainer!T &&
		hasArrayStorageIndexing!T &&
		hasArrayStorageResizing!T &&
		isInputRange!T &&
		hasLength!T;
}

template isMatrixContainer( T ) {
	enum isMatrixContainer =
		isFortranType!(BaseElementType!T) &&
		hasMatrixDimensions!T;
}
		
version( unittest ) {
	static assert( isArrayContainer!(CowArray!double) );
	static assert( isMatrixContainer!(CowMatrix!double) );
	static assert( isVectorStorage!(ArrayStorage!double) );
	static assert( isVectorStorage!(ArrayViewStorage!double) );
	static assert( isVectorStorage!(PackedSubVectorStorage!(SymmetricArrayAdapter!(CowArrayRef!double, MatrixTriangle.Upper, StorageOrder.ColumnMajor), VectorType.Column)) );
	static assert( isVectorStorage!(PackedSubVectorStorage!(TriangularArrayAdapter!(CowArrayRef!double, MatrixTriangle.Upper, StorageOrder.ColumnMajor), VectorType.Column)) );
	static assert( isMatrixStorage!(GeneralMatrixStorage!double) );
	static assert( isMatrixStorage!(SymmetricStorage!double) );
	static assert( isMatrixStorage!(TriangularStorage!double) );
	static assert( isMatrixStorage!(GeneralMatrixViewStorage!double) );
}