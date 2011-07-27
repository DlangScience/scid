/**
$(D RegionAllocator) is a memory allocator based on a thread-local segmented
stack.  A segmented stack is similar to a regular stack in that memory is
allocated and freed in last in, first out order.  When memory is requested from
a segmented stack, it first checks whether enough space is available in the
current segment, and if so increments the stack pointer and returns.  If not,
a new segment is allocated.  When memory is freed, the stack pointer is
decremented.  If the last segment is no longer in use, it may be returned to
where it was allocated from or retained for future use.

$(D RegionAllocator) has the following advantages compared to allocation on the
call stack:

1.  Pointers to memory allocated on the $(D RegionAllocator) stack are still
    valid when the function they were allocated from returns, unless the
    last copy of the $(D RegionAllocator) object they were allocated from
    goes out of scope.  Functions can be written to create and return data
    structures on the $(D RegionAllocator) stack.

2.  Since it is a segmented stack, large allocations can be performed with no
    danger of stack overflow errors.

3.  It guarantees 16-byte alignment of all allocated objects.

It has thet following disadvantages:

1.  The extra checks due to stack segmentation make allocation slightly
    slower and prevent simultanelous deallocation of several objects from
    being a single decrement of the stack pointer.

2.  Memory allocated on RegionAllocator is accessed through an extra layer of
    pointer indirection compared to memory on the call stack.

It has the following advantages compared to heap allocation:

1.  Both allocation and deallocation are extremely fast.  Most allocations
    consist of verifying enough space is available, incrementing a pointer and
    a performing a few cheap bookkeeping operations.  Most deallocations
    consist decrementing a pointer and performing a few cheap bookkeeping
    operations.

2.  The segmented stack is thread-local, so synchronization is only needed
    when a segment needs to be allocated or freed.

3.  Fragmentation is not an issue when allocating memory on the
    $(D RegionAllocator) stack, though it can be an issue when trying to allocate
    a new segment.

It has the following disadvantages compared to heap allocation:

1.  The requirement that memory be freed in last in, first out order.

2.  No automatic garbage collection.

Note:

The first segment of the $(D RegionAllocator) stack is allocated lazily, so
no space is allocated in any thread that does not use $(D RegionAllocator).

Synopsis:
---
void fun() {
    // Create a new RegionAllocator.
    auto alloc = newRegionAllocator();

    // Allocate a temporary array on the RegionAllocator stack.
    auto arr = alloc.uninitializedArray!(double[][][])(8, 6, 7);
    assert(arr.length == 8);
    assert(arr[0].length == 6);
    assert(arr[0][0].length == 7);

    // When alloc goes out of scope, all memory allocated by it is freed.
    // If alloc has been copied, then the memory is freed when all copies
    // have gone out of scope.
}
---

Author:  David Simcha
Copyright:  Copyright (c) 2008-2011, David Simcha.
License:    $(WEB boost.org/LICENSE_1_0.txt, Boost License 1.0)
*/

module scid.internal.regionallocator;

import std.traits, core.memory, std.range, core.exception, std.conv,
    std.algorithm, std.typetuple, std.exception;

static import core.stdc.stdlib;

// This is just for convenience/code readability/saving typing.
private enum ptrSize = (void*).sizeof;

// This was accidentally assumed in a few places and I'm too lazy to fix it
// until I see proof that it needs to be fixed.
static assert(bool.sizeof == 1);

// Memory allocation routines.  These wrap allocate(), free() and realloc(),
// and guarantee alignment.
private enum size_t alignBytes = 16;

/**
This struct provides an interface to the $(D RegionAllocator) functionality
and enforces scoped deletion.  A new instance is created using the
$(D newRegionAllocator) function.  Each instance has reference semantics
in that any copy will allocate from the same memory.  When the last copy
of an instance goes out of scope, all memory allocated via that instance
is freed.  Only the most recently created $(D RegionAllocator) instance
still in existence may be used to allocate and free memory at any given time.
This is checked via assertions in debug mode, but for performance reasons
it is not checked in release mode.

Examples:
---
void foo() {
    auto alloc = newRegionAllocator();
    auto ptr1 = bar(alloc);
    auto ptr2 = alloc.allocate(42);

    // The last instance of the region allocator used to allocate ptr1
    // and ptr2 is going out of scope here.  The memory pointed to by
    // both ptr1 and ptr2 will be freed.
}

void* bar(RegionAllocator alloc) {
    auto ret = alloc.allocate(42);

    auto alloc2 = newRegionAllocator();
    auto ptr3 = alloc2.allocate(42);

    // ptr3 was allocated using alloc2, which is going out of scope.
    // Its memory will therefore be freed.  ret was allocated using alloc.
    // An instance of this RegionAllocator is still alive in foo() after
    // bar() executes.  Therefore, ret will not be freed on returning and
    // is still valid after bar() returns.

    return ret;
}

void* dontDoThis(RegionAllocator alloc) {
    auto alloc2 = newRegionAllocator();

    // Error:  Allocating from a RegionAllocator instance other than the
    // most recently created one that's still alive.
    auto ptr4 = alloc.allocate(42);
}
---
*/
struct RegionAllocator {
private:
    static struct Stack(T) {  // Simple, fast stack w/o error checking.
        private size_t capacity;
        private size_t index;
        private T* data;
        private enum sz = T.sizeof;

        private static size_t max(size_t lhs, size_t rhs) pure {
            return (rhs > lhs) ? rhs : lhs;
        }

        void push(T elem) {
            if (capacity == index) {
                capacity = max(16, capacity * 2);
                data = cast(T*) core.stdc.stdlib.realloc(data, capacity * sz);
            }
            data[index++] = elem;
        }

        T pop() {
            index--;
            auto ret = data[index];
            data[index] = T.init;  // Prevent false ptrs.
            return ret;
        }

        void destroy() {
            if(data) {
                core.stdc.stdlib.free(data);
                data = null;
            }
        }
    }

    static struct Block {
        size_t used = 0;
        void* space = null;
    }

    final class State {
        size_t used;
        void* space;
        size_t totalAllocs;
        void*[] lastAlloc;
        uint nblocks;
        uint nfree;
        size_t regionIndex = size_t.max;

        // inUse holds info for all blocks except the one currently being
        // allocated from.  freelist holds space ptrs for all free blocks.
        Stack!(Block) inUse;
        Stack!(void*) freelist;

        // Add an element to lastAlloc, checking length first.
        void putLast(void* last) {
            if (totalAllocs == lastAlloc.length)
                doubleSize(lastAlloc);
            lastAlloc[totalAllocs] = cast(void*) last;
            totalAllocs++;
        }

        // Hacky use of the same array to store frame indices, reference
        // counts and previous pointers.
        void putLast(size_t num) {
            return putLast(cast(void*) num);
        }

        void destroy() {
            if(space) {
                alignedFree(space);
                space = null;
            }

            if(lastAlloc) {
                core.stdc.stdlib.free(lastAlloc.ptr);
                lastAlloc = null;
            }

            while(inUse.index > 0) {
                auto toFree = inUse.pop();
                alignedFree(toFree.space);
            }

            inUse.destroy();

            while(freelist.index > 0) {
                auto toFree = freelist.pop();
                alignedFree(toFree);
            }

            freelist.destroy();
        }

        ~this() {
            destroy();
        }
    }

    static ~this() {
        if(state) {
            state.destroy();
            state = null;
        }
    }

    enum size_t blockSize = 4 * 1024 * 1024;
    enum size_t nBookKeep = blockSize / .alignBytes;
    static State state;

    static void doubleSize(ref void*[] lastAlloc) {
        size_t newSize = lastAlloc.length * 2;
        void** ptr = cast(void**) core.stdc.stdlib.realloc(
            lastAlloc.ptr, newSize * ptrSize);
        lastAlloc = ptr[0..newSize];
    }

    static State stateInit() {
        State stateRef;
        try { stateRef = new State; } catch { outOfMemory(); }

        with(stateRef) {
            space = alignedallocate(blockSize);

            // We don't need 16-byte alignment for the bookkeeping array.
            lastAlloc = (cast(void**) core.stdc.stdlib.malloc(nBookKeep))
                        [0..nBookKeep / ptrSize];
            nblocks++;
        }

        state = stateRef;
        return stateRef;
    }

    // CTFE function, for static assertions.  Can't use bsr/bsf b/c it has
    // to be usable at compile time.
    static bool isPowerOf2(size_t num) pure {
        int nBitsSet;
        foreach(shift; 0..size_t.sizeof * 8) {
            size_t mask = (cast(size_t) 1) << shift;
            nBitsSet += ((num & mask) > 0);
        }

        return nBitsSet == 1;
    }

    static State getState() {
        State stateRef = state;
        return (stateRef is null) ? stateInit : stateRef;
    }

    static void regionInit(State stateRef) {
        with(stateRef) {
            putLast(regionIndex);
            putLast(1);
            regionIndex = totalAllocs;
        }
    }

    static void regionFree(State stateRef) {
        with(stateRef) {
            while (totalAllocs > regionIndex) {
                free(stateRef);
            }
            totalAllocs -= 2;  // Reference count, frame index.
            regionIndex = cast(size_t) lastAlloc[totalAllocs];
        }
    }

    void incrementRefCount() {
        auto stateRef = getState();
        stateRef.lastAlloc[correctRegionIndex - 1]++;
    }

    void decrementRefCount() {
        auto stateRef = getState();
        stateRef.lastAlloc[correctRegionIndex - 1]--;
        if(cast(size_t) stateRef.lastAlloc[correctRegionIndex - 1] == 0) {
            enforce(stateRef.regionIndex == correctRegionIndex,
                "Cannot free RegionAlloc regions in non-last in first out " ~
                "order.  Did you return a RegionAllocator from a function?"
            );
            regionFree(stateRef);
        }
    }

    static void* allocate()(size_t nBytes, State stateRef) {
        nBytes = allocSize(nBytes);
        with(stateRef) {
            void* ret;
            if (blockSize - used >= nBytes) {
                ret = space + used;
                used += nBytes;
            } else if (nBytes > blockSize) {
                ret = alignedallocate(nBytes);
            } else if (nfree > 0) {
                inUse.push(Block(used, space));
                space = freelist.pop;
                used = nBytes;
                nfree--;
                nblocks++;
                ret = space;
            } else { // Allocate more space.
                inUse.push(Block(used, space));
                space = alignedallocate(blockSize);
                nblocks++;
                used = nBytes;
                ret = space;
            }
            putLast(ret);
            return ret;
        }
    }

    static void free()(State stateRef) {
        with(stateRef) {
            void* lastPos = lastAlloc[--totalAllocs];

            // Handle large blocks.
            if (lastPos > space + blockSize || lastPos < space) {
                alignedFree(lastPos);
                return;
            }

            used = (cast(size_t) lastPos) - (cast(size_t) space);
            if (nblocks > 1 && used == 0) {
                freelist.push(space);
                Block newHead = inUse.pop;
                space = newHead.space;
                used = newHead.used;
                nblocks--;
                nfree++;

                if (nfree >= nblocks * 2) {
                    foreach(i; 0..nfree / 2) {
                        alignedFree(freelist.pop);
                        nfree--;
                    }
                }
            }
        }
    }

    // The region index that should be current anytime this instance is
    // being used.  This is checked for in debug mode to prevent any
    // but the last created RegionAllocator instance from being used.
    size_t correctRegionIndex = size_t.max;

    bool initialized() @property const pure nothrow @safe {
        return correctRegionIndex < size_t.max;
    }

    void ensureInitialized() {
        if(!initialized) {
            auto stateRef = getState();
            regionInit(stateRef);
            correctRegionIndex = stateRef.regionIndex;
        }
    }

public:

    this(this) {
        incrementRefCount();
    }

    ~this() {
        if(initialized) decrementRefCount();
    }

    // Called on initialization of a RegionAllocator.
    static void regionInit() {
        regionInit(getState());
    }

    // Called on destruction of the last instance of a RegionAllocator.
    static void regionFree() {
        regionFree(getState());
    }

    /**
    Assignment of one RegionAllocator to another is disabled.
    */
    void opAssign(RegionAllocator rhs) @disable {
        // BUGS:  DMD Bug 6330:  @disable doesn't do anything.
        assert(0);
    }

    /**
    Obtains a copy of the last created $(D RegionAllocator) instance.
    */
    RegionAllocator lastCreatedInstance() @property {
        auto stateRef = getState();
        if(stateRef.regionIndex == size_t.max) {
            regionInit();
        }

        auto ret = RegionAllocator(stateRef.regionIndex);
        ret.incrementRefCount();
        return ret;
    }

    /**
    Allocates $(D nBytes) bytes on the $(D RegionAllocator) stack.  The last
    block allocated in the current thread can be freed by calling
    $(D RegionAllocator.free) or $(D RegionAllocator.freeLast) or will be
    automatically freed when the last copy of this $(D RegionAllocator)
    instance goes out of scope.  The memory returned by this function is not
    scanned for pointers by the garbage collector unless $(D GC.addRange) is
    called.
    */
    void* allocate()(size_t nBytes) {
        // This is a template to allow assertions to be enabled or disabled
        // based on how the call site is compiled.

        ensureInitialized();
        auto stateRef = getState();
        assert(stateRef.regionIndex == this.correctRegionIndex,
            "Cannot allocate memory from a RegionAllocator that is not " ~
            "currently at the top of the stack." ~ text(stateRef.regionIndex,
            '\t', this.correctRegionIndex)
        );
        return allocate(nBytes, stateRef);
    }

    /**
    Frees the last block of memory allocated in the current thread by
    $(D RegionAllocator).
    */
    void freeLast()() {
        auto stateRef = getState();
        assert(stateRef.regionIndex == this.correctRegionIndex,
            "Cannot allocate memory from a RegionAllocator that is not " ~
            "currently at the top of the stack, or memory that has not been " ~
            "allocated with this instance."
        );
        free(stateRef);
    }

    /**
    Checks that $(D ptr) is a pointer to the last chunk allocated and then
    frees it.  Throws an exception if the last allocated-first freed
    requirement is violated.
    */
    void free()(void* ptr) {
        auto stateRef = getState();
        void* lastPos = stateRef.lastAlloc[stateRef.totalAllocs - 1];
        enforce(ptr is lastPos);
        freeLast();
    }

    /**Allocates an array of type $(D T) on the $(D RegionAllocator) stack.
    The returned array is not scanned for pointers by the garbage collector
    unless $(D GC.addRange) is called.   $(D T) may be a multidimensional array.
    In this case sizes may be specified for any number of dimensions from 1 to
    the number in $(D T).

    Examples:
    ---
    auto alloc = newRegionAllocator();
    double[] arr = alloc.newArray!(double[])(100);
    assert(arr.length == 100);

    double[][] matrix = alloc.newArray!(double[][])(42, 31);
    assert(matrix.length == 42);
    assert(matrix[0].length == 31);
    ---
    */
    auto newArray(T, I...)(I sizes)
    if(allSatisfy!(isIntegral, I)) {

        static void initialize(R)(R toInitialize) {
            static if(isArray!(ElementType!R)) {
                foreach(elem; toInitialize) {
                    initialize(elem);
                }
            } else {
                toInitialize[] = ElementType!(R).init;
            }
        }

        auto ret = uninitializedArray!(T, I)(sizes);
        initialize(ret);
        return ret;
    }

    /**
    Same as $(D uninitializedArray), except skips initialization of elements for
    performance reasons.
    */
    auto uninitializedArray(T, I...)(I sizes)
    if(allSatisfy!(isIntegral, I)) {
        static assert(sizes.length >= 1,
            "Cannot allocate an array without the size of at least the first " ~
            " dimension.");
        static assert(sizes.length <= nDimensions!T,
            to!string(sizes.length) ~ " dimensions specified for a " ~
            to!string(nDimensions!T) ~ " dimensional array.");

        alias typeof(T.init[0]) E;

        auto ptr = cast(E*) allocate(sizes[0] * E.sizeof);
        auto ret = ptr[0..sizes[0]];

        static if(sizes.length > 1) {
            foreach(ref elem; ret) {
                elem = uninitializedArray!(E)(sizes[1..$]);
            }
        }

        return ret;
    }

    /**
    Dummy implementation for compatibility with allocator interface.
    Always returns false.
    */
    static bool resize(void* p, size_t newSize) {
        return false;
    }

    /**
    Returns the number of bytes to which an allocation of size s is
    guaranteed to be aligned.
    */
    static size_t alignBytes(size_t nBytes) {
        return .alignBytes;
    }

    /**
    Returns the number of bytes used to satisfy an allocation request
    of $(D nBytes).  Will return a value greater than or equal to
    $(D nBytes) to account for alignment overhead.
    */
    static size_t allocSize(size_t nBytes) pure nothrow {
        static assert(isPowerOf2(.alignBytes));
        return (nBytes + (.alignBytes - 1)) & (~(.alignBytes - 1));
    }

    /**
    False because memory allocated by this allocator is not automatically
    reclaimed by the garbage collector.
    */
    enum isAutomatic = false;

    /**
    True because, when the last last copy of a $(RegionAllocator) instance
    goes out of scope, the memory it references is automatically freed.
    */
    enum isScoped = true;

    /**
    If memory is freed via $(D free()) instead of $(D freeLast()) then
    the pointer is checked for validity.
    */
    enum freeIsChecked = true;

    /**
    Returns the maximum number of bytes that may be allocated in the
    current stack segment.
    */
    static size_t segmentSlack() @property {
        return blockSize - getState().used;
    }

    /**
    Copies $(D range) to an array.  The array will be located on the
    $(D RegionAllocator) stack and not scanned for pointers by the garbage
    collector under either of the following conditions:

    1.  $(D std.traits.hasIndirections!(ElementType!R)) is false, or

    2.  $(D R) is a builtin array.  In this case $(D range) maintains pointers
        to all elements at least until $(D tempArray) returns, preventing the
        elements from being freed by the garbage collector.  A similar assumption
        cannot be made for ranges other than builtin arrays.

    If neither condition is met, the array is returned on the C heap
    and $(D GC.addRange) is called.  In either case, $(D RegionAllocator.free)
    or $(D RegionAllocator.regionFree) will free the array as if it had been
    allocated on the $(D RegionAllocator) stack.

    Rationale:  The most common reason to call $(D tempArray) on an array is to
                modify its contents inside a function without affecting the
                caller's view.  In this case $(D range) is not modified and
                prevents the elements from being freed by the garbage
                collector.

    Examples:
    ---
    auto alloc = newRegionAllocator();
    auto arr = alloc.tempArray(iota(5));
    assert(arr == [0, 1, 2, 3, 4]);
    ---
     */
    Unqual!(ElementType!(R))[] tempArray(R)(R range)
    if(isInputRange!(R) && (isArray!(R) || !hasIndirections!(ElementType!(R)))) {
        alias ElementType!(R) E;
        alias Unqual!(E) U;
        static if(hasLength!(R)) {
            U[] ret = uninitializedArray!(U[])(range.length);
            copy(range, ret);
            return ret;
        } else {
            auto state = RegionAllocator.getState();
            auto startPtr = allocate(0);
            size_t bytesCopied = 0;

            while(!range.empty) {  // Make sure range interface is being used.
                auto elem = range.front;
                if(state.used + U.sizeof <= RegionAllocator.blockSize) {
                    range.popFront;
                    *(cast(U*) (startPtr + bytesCopied)) = elem;
                    bytesCopied += U.sizeof;
                    state.used += U.sizeof;
                } else {
                    if(bytesCopied + U.sizeof >= RegionAllocator.blockSize / 2) {
                        // Then just heap-allocate.
                        U[] result = (cast(U*) alignedallocate(bytesCopied * 2))
                            [0..bytesCopied / U.sizeof * 2];

                        immutable elemsCopied = bytesCopied / U.sizeof;
                        result[0..elemsCopied] = (cast(U*) startPtr)
                            [0..elemsCopied];
                        finishCopy(result, range, elemsCopied);
                        freeLast();
                        state.putLast(result.ptr);
                        return result;
                    } else {
                        U[] oldData = (cast(U*) startPtr)
                            [0..bytesCopied / U.sizeof];
                        state.used -= bytesCopied;
                        state.totalAllocs--;
                        U[] uninitializedArray = uninitializedArray!(U[])
                            (bytesCopied / U.sizeof + 1);
                        uninitializedArray[0..oldData.length] = oldData[];
                        startPtr = state.space;
                        uninitializedArray[$ - 1] = elem;
                        bytesCopied += U.sizeof;
                        range.popFront;
                    }
                }
            }
            auto rem = bytesCopied % .alignBytes;
            if(rem != 0) {
                auto toAdd = 16 - rem;
                if(state.used + toAdd < RegionAllocator.blockSize) {
                    state.used += toAdd;
                } else {
                    state.used = RegionAllocator.blockSize;
                }
            }
            return (cast(U*) startPtr)[0..bytesCopied / U.sizeof];
        }
    }

    // Ditto but not worth its own ddoc.
    Unqual!(ElementType!(R))[] tempArray(R)(R range)
    if(isInputRange!(R) && !(isArray!(R) || !hasIndirections!(ElementType!(R)))) {
        // Initial guess of how much space to allocate.  It's relatively large b/c
        // the object will be short lived, so speed is more important than space
        // efficiency.
        enum initialGuess = 128;

        alias Unqual!(ElementType!R) E;
        auto arr = (cast(E*) alignedallocate(E.sizeof * initialGuess, true))
            [0..initialGuess];

        finishCopy(arr, range, 0);
        RegionAllocator.getState().putLast(arr.ptr);
        return arr;
    }
}

/**
Creates a new $(D RegionAllocator) region.  $(D RegionAllocator) has reference
semantics in that, when copied, all copies of the same instance allocate from
the same region of memory.  When the last copy of an instance goes out of scope,
all memory allocated by any instance is freed.  See the $(D RegionAllocator)
documentation for examples.
*/
RegionAllocator newRegionAllocator() {
    return RegionAllocator();
}

// Finishes copying a range to a C heap allocated array.  Assumes the first
// half of the input array is stuff already copied and the second half is
// free space.
private void finishCopy(T, U)(ref T[] result, U range, size_t alreadyCopied) {
    void doRealloc() {
        auto newPtr = cast(T*) alignedRealloc(
            result.ptr, result.length * T.sizeof * 2, result.length * T.sizeof
        );
        result = newPtr[0..result.length * 2];
    }

    auto index = alreadyCopied;
    foreach(elem; range) {
        if(index == result.length) doRealloc();
        result[index++] = elem;
    }

    result = result[0..index];
}

unittest {
    RegionAllocator alloc;
    auto arr = alloc.tempArray(iota(5));
    assert(arr == [0, 1, 2, 3, 4]);

    // Create quick and dirty finite but lengthless range.
    static struct Count {
        uint num;
        uint upTo;
        @property size_t front() {
            return num;
        }
        void popFront() {
            num++;
        }
        @property bool empty() {
            return num >= upTo;
        }
    }

    alloc.allocate(1024 * 1024 * 3);
    Count count;
    count.upTo = 1024 * 1025;
    auto asArray = alloc.tempArray(count);
    foreach(i, elem; asArray) {
        assert(i == elem, to!(string)(i) ~ "\t" ~ to!(string)(elem));
    }
    assert(asArray.length == 1024 * 1025);
    alloc.freeLast();
    alloc.freeLast();
    while(RegionAllocator.getState().freelist.index > 0) {
        alignedFree(RegionAllocator.getState().freelist.pop);
    }
}

unittest {
    RegionAllocator alloc;
    double[] arr = alloc.uninitializedArray!(double[])(100);
    assert(arr.length == 100);

    double[][] matrix = alloc.uninitializedArray!(double[][])(42, 31);
    assert(matrix.length == 42);
    assert(matrix[0].length == 31);

    double[][] mat2 = alloc.newArray!(double[][])(3, 1);
    assert(mat2.length == 3);
    assert(mat2[0].length == 1);

    import std.math;
    assert(isNaN(mat2[0][0]));
}

unittest {
    /* Not a particularly good unittest in that it depends on knowing the
     * internals of RegionAllocator, but it's the best I could come up w/.  This
     * is really more of a stress test/sanity check than a normal unittest.*/

    // Make sure state is completely reset.
    if(RegionAllocator.state) RegionAllocator.state.destroy();
    RegionAllocator.state = null;

     // First test to make sure a large number of allocations does what it's
     // supposed to in terms of reallocing lastAlloc[], etc.
     enum nIter =  RegionAllocator.blockSize * 5 / alignBytes;

    {
         RegionAllocator alloc;
         foreach(i; 0..nIter) {
             alloc.allocate(alignBytes);
         }
         assert(RegionAllocator.getState().nblocks == 5,
            to!string(RegionAllocator.getState().nblocks));
         assert(RegionAllocator.getState().nfree == 0);
         foreach(i; 0..nIter) {
            alloc.freeLast();
        }
        assert(RegionAllocator.getState().nblocks == 1);
        assert(RegionAllocator.getState().nfree == 2);

        // Make sure logic for freeing excess blocks works.  If it doesn't this
        // test will run out of memory.
        enum allocSize = RegionAllocator.blockSize / 2;
        foreach(i; 0..50) {
            foreach(j; 0..50) {
                alloc.allocate(allocSize);
            }
            foreach(j; 0..50) {
                alloc.freeLast();
            }
        }

        // Make sure data is stored properly.
        foreach(i; 0..10) {
            alloc.allocate(allocSize);
        }
        foreach(i; 0..5) {
            alloc.freeLast();
        }
        void* space = RegionAllocator.state.space;
        size_t used = RegionAllocator.state.used;

        {
            RegionAllocator alloc2;

            // This array of arrays should not be scanned by the GC because
            // otherwise bugs caused the not having the GC scan certain internal
            // things in RegionAllocator that it should would not be exposed.
            uint[][] arrays = (cast(uint[]*) GC.malloc((uint[]).sizeof * 10,
                               GC.BlkAttr.NO_SCAN))[0..10];
            foreach(i; 0..10) {
                uint[] data = alloc2.uninitializedArray!(uint[])(250_000);
                foreach(j, ref e; data) {
                    e = cast(uint) (j * (i + 1));  // Arbitrary values that can be read back later.
                }
                arrays[i] = data;
            }

            // Make stuff get overwrriten if blocks are getting GC'd when
            // they're not supposed to.
            GC.minimize;  // Free up all excess pools.
            uint[][] foo;
            foreach(i; 0..40) {
                foo ~= new uint[1_048_576];
            }
            foo = null;

            for(size_t i = 9; i != size_t.max; i--) {
                foreach(j, e; arrays[i]) {
                    assert(e == j * (i + 1));
                }
            }
        }

        assert(space == RegionAllocator.state.space,
            text(space, '\t', RegionAllocator.state.space));
        assert(used == RegionAllocator.state.used);
        while(RegionAllocator.state.nblocks > 1 || RegionAllocator.state.used > 0) {
            alloc.freeLast();
        }
    }

    // Test that everything is really getting destroyed properly when
    // destroy() is called.  If not then this test will run out of memory.
    foreach(i; 0..1000) {
        RegionAllocator.state.destroy();
        RegionAllocator.state = null;

        RegionAllocator alloc;
        foreach(j; 0..1_000) {
            auto ptr = alloc.allocate(20_000);
            assert((cast(size_t) ptr) % alignBytes == 0);
        }

        foreach(j; 0..500) {
            alloc.freeLast();
        }
    }
}

private  void outOfMemory()  {
    throw new OutOfMemoryError("Out of memory in RegionAllocator.");
}

private void* alignedallocate(size_t size, bool shouldAddRange = false) {
    // We need (alignBytes - 1) extra bytes to guarantee alignment, 1 byte
    // to store the shouldAddRange flag, and ptrSize bytes to store
    // the pointer to the beginning of the block.
    void* toFree = core.stdc.stdlib.malloc(
        alignBytes + ptrSize + size
    );

    if(toFree is null) outOfMemory();

    // Add the offset for the flag and the base pointer.
    auto intPtr = cast(size_t) toFree + ptrSize + 1;

    // Align it.
    intPtr = (intPtr + alignBytes - 1) & (~(alignBytes - 1));
    auto ret = cast(void**) intPtr;

    // Store base pointer.
    (cast(void**) ret)[-1] = toFree;

    // Store flag.
    (cast(bool*) ret)[-1 - ptrSize] = shouldAddRange;

    if(shouldAddRange) {
        GC.addRange(ret, size);
    }

    return ret;
}

private void alignedFree(void* ptr) {
    // If it was allocated with alignedallocate() then the pointer to the
    // beginning is at ptr[-1].
    auto addedRange = (cast(bool*) ptr)[-1 - ptrSize];

    if(addedRange) {
        GC.removeRange(ptr);
    }

    core.stdc.stdlib.free( (cast(void**) ptr)[-1]);
}

// This is used by RegionAllocator, but I'm not sure enough that its interface
// isn't going to change to make it public and document it.
private void* alignedRealloc(void* ptr, size_t newLen, size_t oldLen) {
    auto storedRange = (cast(bool*) ptr)[-1 - ptrSize];
    auto newPtr = alignedallocate(newLen, storedRange);
    memcpy(newPtr, ptr, oldLen);

    alignedFree(ptr);
    return newPtr;
}
