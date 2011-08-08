//module scid.storage.triangularview;
//
//import scid.common.storagetraits;
//import scid.common.meta;
//import scid.matrix;
//
//enum UnitDiagonal { yes, no }
//
//enum MatrixStorageAdapterType {
//    Symmetric,
//    Triangular,
//    Diagonal
//}
//
//struct MatrixStorageAdapter(
//        Storage_,
//        MatrixStorageAdapterType adapterType_,
//        MatrixTriangle tri_ = MatrixTriangle.Upper,
//        UnitDiagonal diag_ = UnitDiagonal.no
//) {
//    alias Storage_ Storage;	
//    alias BaseElementType!Storage ElementType;
//    alias storageOrderOf!Storage storageOrder;
//    alias tri_ triangle;
//    
//    enum isRowMajor      = ( storageOrder == StorageOrder.RowMajor );
//    enum isUpper         = ( triangle == MatrixTriangle.Upper );
//    enum isSymmetric     = ( adapterType_ == MatrixStorageAdapterType.Symmetric );
//    enum isHermitian     = isSymmetric && isComplexScalar!ElementType;
//    enum hasUnitDiagonal = !isSymmetric && ( diag_ = UnitDiagonal.yes );
//    enum adapterType     = adapterType_;
//    
//    this()( Storage otherStorage ) {
//        storage = otherStorage;
//    }
//    
//    this( A ... )( A args ) if( A.length > 0 && !is( A[ 0 ] : Storage ) ) {
//        storage = Storage( args );
//    }
//    
//    void forceRefAssign( ref typeof( this ) rhs ) {
//        storage.forceRefAssign( rhs.storage );
//    }
//    
//    ref typeof( this ) opAssign( typeof( this ) rhs ) {
//        move( rhs.storage, storage );
//    }
//    
//    ElementType index( size_t i, size_t j ) const
//    in {
//        assert( i < this.rows && j < this.columns, boundsMsg_( i, j ) );
//    } body {
//
//        static if( hasUnitDiagonal )
//            if( i == j )
//                return One!ElementType;
//        
//        static if( isUpper ) {
//            if( i <= j )
//                return storage.index( i, j );
//        } else {
//            if( j <= i )
//                return storage.index( i, j );
//        }
//        
//        static if( isHermitian )
//            return blas.xconj( storage.index( j, i ) );
//        else static if( isSymmetric )
//            return storage.index( j, i );
//        else 
//            return Zero!ElementType;
//    }
//    
//    void indexAssign( string op = "" )( ElementType rhs, size_t i, size_t j )
//    in {
//        assert( i < rows_ && j < cols_, boundsMsg_( i, j ) );
//    } body {
//        static if( isUpper )
//            bool noswap = i <= j;
//        else
//            bool noswap = j <= i;
//        
//        static if( isSymmetric ) {		
//            static if( isHermitian ) {
//                if( noswap )
//                    rhs = blas.xconj( rhs );
//            }
//            
//            if( noswap )
//                storage.indexAssign!op( rhs, i, j );
//            else
//                storage.indexAssign!op( rhs, j, i );
//            
//        } else {
//            assert( noswap && (!hasUnitDiagonal || i!=j),
//                   format("Assignment to read-only element in triangular matrix view (%d,%d).", i, j ) );
//            
//            storage.indexAssign!op( rhs, i, j );
//        }
//    }
//    
//    RowView row( size_t i )
//    in {
//        assert( i < rows, sliceMsg_(i,0,i,columns) );
//    } body {
//        return subVector_( storage.row( i ), i );
//    }
//    
//    ColumnView column( size_t j )
//    in {
//        assert( j < columns, sliceMsg_(0,j,rows,j) );
//    } body {
//        return subVector_( storage.column( j ), j );
//    }
//    
//    typeof( this ) slice( size_t rowStart, size_t colStart, size_t rowEnd, size_t colEnd )
//    in {
//        assert( rowStart < rowEnd && rowEnd <= this.rows && colStart < colEnd && colEnd <= this.columns,
//               sliceMsg_( rowStart, colStart, rowEnd, colEnd ) );
//        assert( rowStart == colStart && rowEnd == colEnd,
//            format( "%sNon-square slice of triangular/symmetric view (%d,%d)-(%d,%d)",
//                msgPrefix_, rowStart, colStart, rowEnd, colEnd )
//        );
//    } body {
//        return typeof( this )( stor.slice( rowStart, rowStart, rowEnd, rowEnd ) );
//    }
//    
//    typeof( this ) view( size_t rowStart, size_t colStart, size_t rowEnd, size_t colEnd )
//    in {
//        assert( rowStart < rowEnd && rowEnd <= this.rows && colStart < colEnd && colEnd <= this.columns,
//               sliceMsg_( rowStart, colStart, rowEnd, colEnd ) );
//        assert( rowStart == colStart && rowEnd == colEnd,
//            format( "%sNon-square view of triangular/symmetric view (%d,%d)-(%d,%d)",
//                msgPrefix_, rowStart, colStart, rowEnd, colEnd )
//        );
//    } body {
//        return typeof( this )( stor.slice( rowStart, rowStart, rowEnd, rowEnd ) );
//    }
//    
//    void copy()           { static assert( false, "Use fallbackCopy." );           }
//    void scale()          { static assert( false, "Use fallbackScale." );          }
//    void scaledAddition() { static assert( false, "Use fallbackScaledAddition." ); }
//    
//    void invert() {
//        enum uplo = isUpper ? 'U' : 'L';
//        enum diag = hasUnitDiagonal ? 'U' : 'N';
//        int info = 0;
//        size_t n = this.rows;
//        static if( isSymmetric ) {
//            auto alloc = newRegionAllocator();
//            auto ipiv  = alloc.uninitializedArray!( int[] )( this.rows );
//            auto work  = alloc.uninitializedArray!( ElementType[] )( this.rows );
//            // TODO: copy me.
//            static if( isHermitian ) {
//                lapack.hetrf!( uplo, n, this.data, this.leading, ipiv, work, work.length, info );
//                lapack.hetri!( uplo, n, this.data, this.leading, ipiv, work, info );
//            } else {
//                lapack.sytrf!( uplo, n, this.data, this.leading, ipiv, work, work.length, info );
//                lapack.sytri!( uplo, n, this.data, this.leading, ipiv, work, info );
//            }
//        } else {
//            lapack.trtri!( uplo, diag )( n, this.data, this.leading, info );
//        }
//        
//        assert( info == 0,
//               format( "invert: Inversion of singular triangular view matrix." ) );
//    }
//    
//    void solveRight( Transpose transM, Side side_, Dest )( auto ref Dest dest )
//            if( isStridedVectorStorage!Dest || isGeneralMatrixStorage!Dest ) {
//        
//        enum vectorRhs = isStridedVectorStorage!Dest;
//        enum uplo      = isUpper ? 'U' : 'L';
//        enum side      = side_ == Side.Left ? 'L' : 'R';
//        size_t n       = this.rows;
//        
//        static if( isSymmetric ) {
//            static assert( false, "Symmetric Views not fully implemented yet." );
//        } else {
//            enum diag = hasUnitDiagonal ? 'U' : 'N';
//            static if( isVectorRhs ) {
//                static if( isComplex!ElementType ) {
//                    static if( side == 'R' && transM ) {
//                        // annoying special case v * A.H = ((A.H).T * x.T).T = (conj(A) * x)
//                        // to solve this we'll call trsm if stride == 1, or do conj(A*conj(x)) otherwise.
//                        if( dest.stride == 1 ) {
//                            blas.trsm!( 'R', uplo, 'C', diag )( dest.rows, dest.columns, this.cdata, this.leading, dest.data, dest.leading );
//                        } else {
//                            blas.xcopyc( dest.length, dest.data, dest.stride );
//                            blas.trsv!( uplo, trans, diag )( this.cdata, this.leading, dest.data, dest.stride );
//                            blas.xcopyc( dest.length, dest.data, dest.stride );
//                        }
//                    } else {
//                        static if( side == 'R' )
//                            enum trans = 'T'; // we know transM == false
//                        else
//                            enum trans = transM ? 'C' : 'N';
//                        
//                        blas.trsv!( uplo, trans, diag )( this.cdata, this.leading, dest.data, dest.stride );
//                    }
//                } else {
//                    enum trans = side == 'R' ^ transM ? 'T' : 'N';
//                    blas.trsv!( uplo, trans, diag )( this.cdata, this.leading, dest.data, dest.stride );
//                }
//            } else {
//                // TODO: Treat different storageOrders.
//                static assert( storageOrderOf!Dest == storageOrder,
//                    typeof(this).stringof ~ "Different storage orders in solveRight not supported yet." );
//                enum trans = transM ? (isComplex!ElementType ? 'C' : 'T') : 'N';
//                blas.trsm!( side, uplo, trans, diag )(
//                    dest.rows, dest.columns, One!T, this.cdata, this.leading, dest.data, dest.leading
//                );
//            }
//        }
//    }
//    
//    @property {
//        // size_t              size()    const { return isInitd_() ? storage.rows    : 0;    }
//        // size_t              leading() const { return isInitd_() ? storage.leading : 0;    }
//        // ElementType*        data()          { return isInitd_() ? storage.data    : null; }
//        // const(ElementType)* cdata()   const { return isInitd_() ? storage.cdata   : null; }
//    }
//    
//    private import scid.storage.generalmat;
//    template Promote( Other ) {
//        static if( is( Other : typeof(this) ) ) {
//            alias typeof( this ) Promote;
//        } else {
//            alias Promotion!( BasicGeneralMatrixStorage!(MatrixTypeOf!Storage), Other ) Promote;
//        }
//    }
//    
//    static struct SubVector( VectorType vtype ) {
//        alias typeof( Storage.init.row( 0 ) )    StorageRow;
//        alias typeof( Storage.init.column( 0 ) ) StorageColumn;
//        
//        static assert( is( StorageRow.ContainerRef : StorageColumn.ContainerRef ) );
//        alias StorageRow.ContainerRef ContainerRef;
//        alias BaseElementType!ContainerRef ElementType;
//        alias vtype vectorType;
//        alias SubVector!( transposeVectorType!vectorType ) Transposed;
//        
//        enum isRow = vectorType == VectorType.Row;
//        
//        static if( isRow ) {
//            alias StorageRow StorageA;
//            static if( isSymmetric )
//                alias StorageColumn StorageB;
//        } else {
//            alias StorageColumn StorageA;
//            static if( isSymmetric )
//                alias StorageRow StorageB;
//        }
//        
//        ElementType index( size_t i ) const {
//            static if( isRow ^ isUpper ) {
//                if( i < storageA_.length )
//                    return storageA_.index( i );
//                
//                static if( isSymmetric ) {
//                    i -= storageA_.length;
//                    static if( isHermitian )
//                        return xconj( storageB_.index( i - storageA_.length ) );
//                    else
//                        return storageB_.index( i - storageA_.length );
//                } else
//                    return Zero!ElementType;
//            } else {
//                    
//            }
//            
//            
//        }
//        
//        void indexAssign( string op = "" )( ElementType rhs, size_t i ) {
//            static if( isSymmetric ) {
//                if( i < storageA_.length )
//                    return storageA_.index( i );
//                
//                if( isHermitian ) {
//                    i -= 	
//                }
//                
//            }
//        }
//        
//        
//        StorageA storageA_;
//        static if( isSymmetric )
//            StorageB storageB_;
//        
//    }
//    
//    Storage storage;
//}