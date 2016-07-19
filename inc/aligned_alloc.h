//Copyright (c) 2016, Christopher Weaver
//All rights reserved.
//
//Redistribution and use in source and binary forms, with or without modification,
//are permitted provided that the following conditions are met:
//
//1. Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
//2. Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
//THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
//IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
//INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
//NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
//WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
//ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//POSSIBILITY OF SUCH DAMAGE.

#ifndef NUSQUIDS_ALIGNED_ALLOC_H
#define NUSQUIDS_ALIGNED_ALLOC_H

#include <limits>

//An allocator which allocates blocks aligned to 2^n bytes for a user selected n.
//The user is responsible for picking n be at least large enough for the natural 
//alignment of the type to be stored in the allocated region!
//For simplicity, this implementation doesn't bother using facilities which may 
//already be present like posix_memalign. This could be added. 
template <typename T>
class aligned_allocator{
public:
	
	//Types

	typedef size_t size_type;
	typedef ptrdiff_t difference_type;
	typedef T* pointer;
	typedef const T* const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef T value_type;
	template <class U> struct rebind{
		typedef aligned_allocator<U> other;
	};
	typedef std::true_type propagate_on_container_move_assignment;
	
	//'structors
	
	aligned_allocator(uint8_t a) noexcept:
	align_pow(a)
	{
		assert(a<32 && "2^31 is the maximum supported alignment");
	}
	
	aligned_allocator(const aligned_allocator& other) noexcept:
	align_pow(other.align_pow){}
	
	template <class U>
	aligned_allocator(const aligned_allocator<U>& other) noexcept:
	align_pow(other.align_pow){}
	
	~aligned_allocator(){}
	
	//Helper functions
	
	pointer address(reference x) const noexcept{
		return(&x); //nothing clever required
	}
	
	const_pointer address(const_reference x) const noexcept{
		return(&x); //nothing clever required
	}
	
	///\returns The largest value N for which the call allocate(N,0) might succeed.
	size_type max_size(){
		return((std::numeric_limits<size_type>::max()-get_align())/sizeof(T));
	}
	
	template <class U, class... Args>
	void construct(U* p, Args&&... args){
		::new((void *)p) U(std::forward<Args>(args)...);
	}
	
	template <class U>
	void destroy(U* p){
		p->~U();
	}
	
	template<typename U>
	bool operator==(const aligned_allocator<U>& other) const noexcept{
		return(align_pow==other.align_pow);
	}
	
	template<typename U>
	bool operator!=(const aligned_allocator<U>& other) const noexcept{
		return(align_pow!=other.align_pow);
	}
	
	//Main interface
	
	pointer allocate(size_type n, const_pointer hint = 0){
		//can't allocate too many objects
		if(n>max_size()) 
			throw std::bad_alloc();
		
		size_type needed=n*sizeof(T);
		size_type alignment=get_align();
		needed+=alignment; //include space to do our thing
		
		//try to get memory
		void* ptr=new char[needed];
		
		//align!
		uint8_t offset=alignment-reinterpret_cast<size_type>(ptr)%alignment;
		//Note that since align is a power of 2, align-1 forms a mask with the low
		//align_pow bits set. It's negation will then preserve only the higher bits
		//of ptr, effectively flooring it to the nearest multiple of align. Adding
		//align then preserves this alignment, and ensures that at least one byte 
		//of padding is before ptr. 
		//This trick borrowed from Johan Mabille
		//https://jmabille.github.io/blog/2014/12/06/aligned-memory-allocator/
		ptr=reinterpret_cast<void*>((reinterpret_cast<size_type>(ptr) & ~(alignment-1))+alignment);
		//Stash the offset in our guaranteed byte of padding where we can find it later
		*((uint8_t*)ptr-1)=offset;
		
		return(pointer(ptr));
	}
	
	void deallocate(pointer p, size_type n){
		if(!p) //ignore NULL pointers
			return;
		
		uint8_t* ptr=(uint8_t*)p;
		//find our stashed offset value
		uint8_t offset=*(ptr-1);
		//recover the real pointer
		ptr-=offset;
		delete[] ptr;
	}
	
private:
	///Log base 2 of alignment, in bytes. 
	///Must be <32. 
	uint8_t align_pow;
	
	///Compute the actual alignment
	uint32_t get_align() const noexcept{
		return(1UL << align_pow);
	}
};

#endif //NUSQUIDS_ALIGNED_ALLOC_H
