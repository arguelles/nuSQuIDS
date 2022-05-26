//Copyright (c) 2015, Christopher Weaver
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

#ifndef MARRAY_H
#define MARRAY_H

#if __cplusplus < 201103L
#error C++11 compiler required. Update your compiler and use the flag -std=c++11
#endif

#include <array>
#include <cassert>
#include <functional>
#include <numeric>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>

#include <iosfwd>

#ifdef MARRAY_USE_LIBDIVIDE
#include "libdivide.h"
#endif

namespace nusquids{

namespace detail{
	
//This simple compressed pair is based heavily on the implementation in libc++
//It does not bother to handle:
//	the two element types being the same
//	either element type being both empty and final
//	constructing element types in any way other than copy construction

template<typename T1, typename T2, bool=std::is_empty<T1>::value, bool=std::is_empty<T2>::value>
struct compressed_pair_base;
	
template<typename T1, typename T2>
struct compressed_pair_base<T1,T2,false,false>{
private:
	T1 first_;
	T2 second_;
public:
	compressed_pair_base();
	compressed_pair_base(T1 t1, T2 t2):
	first_(std::move(t1)),second_(std::move(t2)){}
	compressed_pair_base(const compressed_pair_base& p)
	noexcept(std::is_nothrow_copy_constructible<T1>::value && std::is_nothrow_copy_constructible<T2>::value):
	first_(p.first_),second_(p.second_){}
	compressed_pair_base& operator=(const compressed_pair_base& p)
	noexcept(std::is_nothrow_copy_assignable<T1>::value && std::is_nothrow_copy_assignable<T2>::value){
		first_=p.first_;
		second_=p.second_;
		return(*this);
	}
	compressed_pair_base(compressed_pair_base&& p)
	noexcept(std::is_nothrow_move_constructible<T1>::value && std::is_nothrow_move_constructible<T2>::value):
	first_(std::move(p.first_)),second_(std::move(p.second_)){}
	compressed_pair_base& operator=(compressed_pair_base&& p)
	noexcept(std::is_nothrow_move_assignable<T1>::value && std::is_nothrow_move_assignable<T2>::value){
		first_=std::move(p.first_);
		second_=std::move(p.second_);
		return(*this);
	}
	typename std::add_lvalue_reference<T1>::type first() noexcept{ return(first_); }
	typename std::add_lvalue_reference<const T1>::type first() const noexcept{ return(first_); }
	typename std::add_lvalue_reference<T2>::type second() noexcept{ return(second_); }
	typename std::add_lvalue_reference<const T2>::type second() const noexcept{ return(second_); }
	void swap(compressed_pair_base& p) /*missing noexcept*/{
		using std::swap;
		swap(first_,p.first_);
		swap(second_,p.second_);
	}
};
	
//template<typename T1, typename T2>
//struct compressed_pair_base<T1,T2,true,false> : private T1{
//private:
//	T2 second_;
//public:
//	compressed_pair_base();
//	compressed_pair_base(T1 t1, T2 t2):
//	T1(std::move(t1)),second_(std::move(t2)){}
//	compressed_pair_base(const compressed_pair_base& p)
//	noexcept(std::is_nothrow_copy_constructible<T1>::value && std::is_nothrow_copy_constructible<T2>::value):
//	T1(p.first_),second_(p.second_){}
//	compressed_pair_base& operator=(const compressed_pair_base& p)
//	noexcept(std::is_nothrow_copy_assignable<T1>::value && std::is_nothrow_copy_assignable<T2>::value){
//		T1::operator=(p.first());
//		second_=p.second_;
//		return(*this);
//	}
//	compressed_pair_base(compressed_pair_base&& p)
//	noexcept(std::is_nothrow_move_constructible<T1>::value && std::is_nothrow_move_constructible<T2>::value):
//	T1(std::move(p.first_)),second_(std::move(p.second_)){}
//	compressed_pair_base& operator=(compressed_pair_base&& p)
//	noexcept(std::is_nothrow_move_assignable<T1>::value && std::is_nothrow_move_assignable<T2>::value){
//		T1::operator=(std::move(p.first()));
//		second_=std::move(p.second_);
//		return(*this);
//	}
//	typename std::add_lvalue_reference<T1>::type first() noexcept{ return(*this); }
//	typename std::add_lvalue_reference<const T1>::type first() const noexcept{ return(*this); }
//	typename std::add_lvalue_reference<T2>::type second() noexcept{ return(second_); }
//	typename std::add_lvalue_reference<const T2>::type second() const noexcept{ return(second_); }
//	void swap(compressed_pair_base& p) /*missing noexcept*/{
//		using std::swap;
//		swap(second_,p.second_);
//	}
//};
	
template<typename T1, typename T2>
struct compressed_pair_base<T1,T2,false,true> : private T2{
private:
	T1 first_;
public:
	compressed_pair_base();
	compressed_pair_base(T1 t1, T2 t2):
	T2(std::move(t2)),first_(std::move(t1)){}
	compressed_pair_base(const compressed_pair_base& p)
	noexcept(std::is_nothrow_copy_constructible<T1>::value && std::is_nothrow_copy_constructible<T2>::value):
	T2(p.second_),first_(p.first_){}
	compressed_pair_base& operator=(const compressed_pair_base& p)
	noexcept(std::is_nothrow_copy_assignable<T1>::value && std::is_nothrow_copy_assignable<T2>::value){
		first_=p.first_;
		T2::operator=(p.second());
		return(*this);
	}
	compressed_pair_base(compressed_pair_base&& p)
	noexcept(std::is_nothrow_move_constructible<T1>::value && std::is_nothrow_move_constructible<T2>::value):
	first_(std::move(p.first_)),T2(std::move(p.second_)){}
	compressed_pair_base& operator=(compressed_pair_base&& p)
	noexcept(std::is_nothrow_move_assignable<T1>::value && std::is_nothrow_move_assignable<T2>::value){
		first_=std::move(p.first_);
		T2::operator=(std::move(p.second()));
		return(*this);
	}
	typename std::add_lvalue_reference<T1>::type first() noexcept{ return(first_); }
	typename std::add_lvalue_reference<const T1>::type first() const noexcept{ return(first_); }
	typename std::add_lvalue_reference<T2>::type second() noexcept{ return(*this); }
	typename std::add_lvalue_reference<const T2>::type second() const noexcept{ return(*this); }
	void swap(compressed_pair_base& p) /*missing noexcept*/{
		using std::swap;
		swap(first_,p.first_);
	}
};
	
//template<typename T1, typename T2>
//struct compressed_pair_base<T1,T2,true,true> : private T1, private T2{
//public:
//	compressed_pair_base();
//	compressed_pair_base(T1 t1, T2 t2):
//	T1(std::move(t1)),T2(std::move(t2)){}
//	compressed_pair_base(const compressed_pair_base& p)
//	noexcept(std::is_nothrow_copy_constructible<T1>::value && std::is_nothrow_copy_constructible<T2>::value):
//	T1(p.first_),T2(p.second_){}
//	compressed_pair_base& operator=(const compressed_pair_base& p)
//	noexcept(std::is_nothrow_copy_assignable<T1>::value && std::is_nothrow_copy_assignable<T2>::value){
//		T1::operator=(p.first());
//		T2::operator=(p.second());
//		return(*this);
//	}
//	compressed_pair_base(compressed_pair_base&& p)
//	noexcept(std::is_nothrow_move_constructible<T1>::value && std::is_nothrow_move_constructible<T2>::value):
//	T1(std::move(p.first_)),T2(std::move(p.second_)){}
//	compressed_pair_base& operator=(compressed_pair_base&& p)
//	noexcept(std::is_nothrow_move_assignable<T1>::value && std::is_nothrow_move_assignable<T2>::value){
//		T1::operator=(std::move(p.first()));
//		T2::operator=(std::move(p.second()));
//		return(*this);
//	}
//	typename std::add_lvalue_reference<T1>::type first() noexcept{ return(*this); }
//	typename std::add_lvalue_reference<const T1>::type first() const noexcept{ return(*this); }
//	typename std::add_lvalue_reference<T2>::type second() noexcept{ return(*this); }
//	typename std::add_lvalue_reference<const T2>::type second() const noexcept{ return(*this); }
//	void swap(compressed_pair_base& p) noexcept{}
//};

template<typename T1, typename T2>
class compressed_pair : private compressed_pair_base<T1,T2>{
private:
	using base=compressed_pair_base<T1,T2>;
public:
	compressed_pair(){}
	compressed_pair(T1 t1, T2 t2):base(std::move(t1),std::move(t2)){}
	compressed_pair(const compressed_pair& p)
	noexcept(std::is_nothrow_copy_constructible<base>::value):
	base(p){}
	compressed_pair& operator=(const compressed_pair& p)
	noexcept(std::is_nothrow_copy_assignable<base>::value){
		base::operator=(p);
		return(*this);
	}
	compressed_pair(compressed_pair&& p)
	noexcept(std::is_nothrow_move_constructible<base>::value):
	base(p){}
	compressed_pair& operator=(compressed_pair&& p)
	noexcept(std::is_nothrow_move_assignable<base>::value){
		base::operator=(std::move(p));
		return(*this);
	}
	typename std::add_lvalue_reference<T1>::type first() noexcept{ return(base::first()); }
	typename std::add_lvalue_reference<const T1>::type first() const noexcept{ return(base::first()); }
	typename std::add_lvalue_reference<T2>::type second() noexcept{ return(base::second()); }
	typename std::add_lvalue_reference<const T2>::type second() const noexcept{ return(base::second()); }
	void swap(compressed_pair& p) /*missing noexcept*/{ base::swap(p); }
};
	
template<typename T1, typename T2>
void swap(compressed_pair<T1,T2>& p1, compressed_pair<T1,T2>& p2){ p1.swap(p2); }

template<typename Alloc>
Alloc select_on_container_move_construction(Alloc& a, std::false_type){ return(a); }
template<typename Alloc>
Alloc select_on_container_move_construction(Alloc& a, std::true_type){ return(std::move(a)); }
template<typename Alloc>
Alloc select_on_container_move_construction(Alloc& a){
	using namespace std;
	return(select_on_container_move_construction(a,integral_constant<bool,allocator_traits<Alloc>::propagate_on_container_move_assignment::value>()));
}

template<typename Alloc>
void select_on_container_move_assignment(Alloc& src, Alloc& trg, std::false_type) noexcept{/*do nothing*/}
template<typename Alloc>
void select_on_container_move_assignment(Alloc& src, Alloc& trg, std::true_type) noexcept{ trg=std::move(src); }
template<typename Alloc>
void select_on_container_move_assignment(Alloc& src, Alloc& trg) noexcept{
	using namespace std;
	return(select_on_container_move_assignment(src,trg,integral_constant<bool,allocator_traits<Alloc>::propagate_on_container_move_assignment::value>()));
}

template<typename T, unsigned int N, unsigned int FullRank>
struct subscript_proxy{
private:
	T* data;
	const std::array<size_t,FullRank>& extents;
	const size_t offset;
public:
	subscript_proxy(T* data, const std::array<size_t,FullRank>& extents, size_t offset):
	data(data),extents(extents),offset(offset*extents[FullRank-N]){}
	
	subscript_proxy<T,N-1,FullRank> operator[](size_t i){
		return(subscript_proxy<T,N-1,FullRank>{data,extents,offset+i});
	}
};

template<typename T,unsigned int FullRank>
struct subscript_proxy<T,1,FullRank>{
private:
	T* data;
public:
	subscript_proxy(T* data, const std::array<size_t,FullRank>& extents, size_t off):
	data(data+off*extents[FullRank-1]){}
	
	T& operator[](size_t i){
		return(*(data+i));
	}
};

template<typename T, size_t Rank>
struct subscript_traits{
	using result_type=subscript_proxy<T,Rank-1,Rank>;
	static result_type make_subscript(T* data, const std::array<size_t,Rank>& extents, size_t offset){
		return(result_type(data,extents,offset));
	}
};

template<typename T>
struct subscript_traits<T,1>{
	using result_type=T&;
	template<size_t FullRank>
	static result_type make_subscript(T* data, const std::array<size_t,FullRank>& extents, size_t offset){
		return(*(data+offset));
	}
};
	
template<typename Target, typename Source>
using match_const=typename std::conditional<std::is_const<Source>::value,
                                            typename std::add_const<Target>::type,
                                            typename std::remove_const<Target>::type>;

template<typename T, unsigned int Rank>
struct initializer_for_rank{
	using type=std::initializer_list<typename initializer_for_rank<T,Rank-1>::type>;
};
template<typename T>
struct initializer_for_rank<T,1>{
	using type=std::initializer_list<T>;
};

template<typename T, unsigned int Rank, unsigned int IRank>
struct initializer_shape_checker{
	static void check(const std::array<size_t,Rank>& extents, const typename detail::initializer_for_rank<T,IRank>::type& init, unsigned int entry){
		if(init.size()!=extents[Rank-IRank])
			throw std::logic_error("Incorrect shape for multidimensional array initializer:\nentry "+
								   std::to_string(entry)+" in dimension "+std::to_string(Rank-IRank)+
								   " has size "+std::to_string(init.size())+" but should have size "+
								   std::to_string(extents[Rank-IRank]));
		unsigned int subentry=0;
		for(const auto& i : init)
			initializer_shape_checker<T,Rank,IRank-1>::check(extents,i,subentry++);
	}
};
template<typename T, unsigned int Rank>
struct initializer_shape_checker<T,Rank,1>{
	static constexpr unsigned int IRank=1;
	static void check(const std::array<size_t,Rank>& extents, const typename detail::initializer_for_rank<T,IRank>::type& init, unsigned int entry){
		if(init.size()!=extents[Rank-IRank])
			throw std::logic_error("Incorrect shape for multidimensional array initializer:\nentry "+
								   std::to_string(entry)+" in dimension "+std::to_string(Rank-IRank)+
								   " has size "+std::to_string(init.size())+" but should have size "+
								   std::to_string(extents[Rank-IRank]));
	}
};

template<typename T, unsigned int Rank, typename Allocator, unsigned int IRank>
struct multidimensional_initializer{
	static void initialize(const typename detail::initializer_for_rank<T,IRank>::type& init, Allocator& alloc, T*& write_ptr){
		for(const auto& i : init)
			multidimensional_initializer<T,Rank,Allocator,IRank-1>::initialize(i,alloc,write_ptr);
	}
	static void assign(const typename detail::initializer_for_rank<T,IRank>::type& init, T*& write_ptr){
		for(const auto& i : init)
			multidimensional_initializer<T,Rank,Allocator,IRank-1>::assign(i,write_ptr);
	}
};
template<typename T, unsigned int Rank, typename Allocator>
struct multidimensional_initializer<T,Rank,Allocator,1>{
	static void initialize(const typename detail::initializer_for_rank<T,1>::type& init, Allocator& alloc, T*& write_ptr){
		for(const auto& i : init)
			std::allocator_traits<Allocator>::construct(alloc,write_ptr++,i);
	}
	static void assign(const typename detail::initializer_for_rank<T,1>::type& init, T*& write_ptr){
		for(const auto& i : init)
			*(write_ptr++)=i;
	}
};
	
template<typename... > struct valid_types_impl { using type = void; };
template<typename... Types> using valid_types = typename valid_types_impl<Types...>::type;
	
template<typename C, typename = void>
struct is_index_container : std::false_type{};
	
//An index container is something we can iterate over, whose iterators refer to size_ts
template<typename C>
struct is_index_container<C,
		valid_types<typename C::value_type
		            ,decltype(std::declval<C>().size())
		            ,decltype(std::declval<C>().begin())
		            ,decltype(std::declval<C>().end())
		>
	>
	: std::true_type{};
	
template<unsigned int a, unsigned int b>
struct Max{
	static constexpr unsigned int value=(a>b ? a : b);
};
	
template<typename T, typename = void>
struct is_marray : std::false_type{};
template<typename T>
struct is_marray<T,valid_types<decltype(T::marray_tag)>> : std::true_type{};

//TODO: this really needs to go somewhere else
struct no_init_tag{};
static constexpr no_init_tag no_init={};
	
} //namespace detail

template<typename T, unsigned int Rank, typename Alloc = std::allocator<T>>
class marray{
public:
	static_assert(Rank>0,"Multidimensional arrays must have at least one dimension");
	
	using value_type=T;
	using allocator_type=Alloc;
	using reference=T&;
	using const_reference=const T&;
	
	using allocator_traits=typename std::allocator_traits<allocator_type>;
	
	using pointer=typename allocator_traits::pointer;
	using const_pointer=typename allocator_traits::const_pointer;
	using size_type=size_t;
	
	enum{marray_tag};
	
private:
	template<typename DerefType, typename DerivedType, typename PointerType=typename std::add_pointer<DerefType>::type>
	struct iterator_base : public std::iterator<std::random_access_iterator_tag, DerefType, std::ptrdiff_t, PointerType>{
	private:
		using derived_type=DerivedType;
		using base_type=std::iterator<std::random_access_iterator_tag, DerefType, std::ptrdiff_t, PointerType>;
		friend class marray;
	protected:
		using array_type=typename detail::match_const<marray,DerefType>::type;
		PointerType ptr;
		iterator_base(array_type* a, size_type i):ptr(a->data+i){}
	public:
		enum{marray_iterator_tag};
		using typename base_type::value_type;
		using typename base_type::difference_type;
		using typename base_type::pointer;
		using typename base_type::reference;
		
		iterator_base():ptr(nullptr){}
		iterator_base(const iterator_base<typename std::remove_const<value_type>::type,DerivedType>& other):
		ptr(other.ptr){}
		bool operator==(const iterator_base& other) const{
			return(ptr==other.ptr);
		}
		bool operator!=(const iterator_base& other) const{
			return(ptr!=other.ptr);
		}
		reference operator*(){ return(*ptr); }
		pointer operator->(){ return(ptr); }
		const DerefType* operator->() const{ return(ptr); }
		///prefix increment
		derived_type& operator++(){
			ptr=static_cast<derived_type*>(this)->increment(1);
			return(*static_cast<derived_type*>(this));
		}
		///postfix increment
		derived_type operator++(int){
			derived_type old=*static_cast<derived_type*>(this);
			ptr=static_cast<derived_type*>(this)->increment(1);
			return(old);
		}
		///prefix decrement
		derived_type& operator--(){
			ptr=static_cast<derived_type*>(this)->increment(-1);
			return(*static_cast<derived_type*>(this));
		}
		///postfix decrement
		derived_type operator--(int){
			derived_type old=*static_cast<derived_type*>(this);
			ptr=static_cast<derived_type*>(this)->increment(-1);
			return(old);
		}
		derived_type& operator+=(difference_type n){
			ptr=static_cast<derived_type*>(this)->increment(n);
			return(*static_cast<derived_type*>(this));
		}
		derived_type& operator-=(difference_type n){
			ptr=static_cast<derived_type*>(this)->increment(-n);
			return(*static_cast<derived_type*>(this));
		}
		derived_type operator+(difference_type n) const{
			return(derived_type(static_cast<const derived_type&>(*this))+=n);
		}
		derived_type operator-(difference_type n) const{
			return(derived_type(static_cast<const derived_type&>(*this))-=n);
		}
		reference operator[](difference_type n){
			return(*(static_cast<derived_type*>(this)->increment(n)));
		}
	};
	
	void copy_assign(const marray& other, std::true_type do_replace_allocator){
		allocator_type new_allocator(other.allocator());
		value_type* new_data=allocator_traits::allocate(new_allocator,other.total_size()); //may throw
		try{
			copy_init_buffer(new_data,other.total_size(),new_allocator,other.data); //may throw
		}catch(...){
			allocator_traits::deallocate(allocator(),new_data,total_size());
			throw;
		}
		//----
		std::swap(data,new_data);
		destroy_buffer(new_data,total_size(),allocator());
		allocator()=std::move(new_allocator);
	}
	void copy_assign(const marray& other, std::false_type dont_replace_allocator){
		value_type* new_data=allocator_traits::allocate(allocator(),other.total_size()); //may throw
		try{
			copy_init_buffer(new_data,other.total_size(),allocator(),other.data); //may throw
		}catch(...){
			allocator_traits::deallocate(allocator(),new_data,total_size());
			throw;
		}
		//----
		std::swap(data,new_data);
		destroy_buffer(new_data,total_size(),allocator());
	}
	
	void move_assign(marray&& other, std::true_type do_replace_allocator){
		using std::swap;
		swap(extents,other.extents);
		swap(data,other.data);
		swap(total_size(),other.total_size());
		allocator()=std::move(other.allocator());
	}
	void move_assign(marray&& other, std::false_type dont_replace_allocator){
		using std::swap;
		if(allocator()==other.allocator()){
			swap(extents,other.extents);
			swap(data,other.data);
			swap(total_size(),other.total_size());
		}
		else{
			value_type* new_data=allocator_traits::allocate(allocator(),other.total_size()); //may throw
			try{
				copy_init_buffer(new_data,other.total_size(),allocator(),other.data); //may throw
			}catch(...){
				allocator_traits::deallocate(allocator(),new_data,total_size());
				throw;
			}
			//----
			swap(data,new_data);
			destroy_buffer(new_data,total_size(),allocator());
			std::copy_n(other.extents.begin(),Rank,extents.begin());
			total_size()=other.total_size();
		}
	}
	
	void swap_allocator(allocator_type& other_alloc, std::true_type do_replace_allocator){
		using std::swap;
		swap(allocator(),other_alloc);
	}
	void swap_allocator(allocator_type& other_alloc, std::false_type dont_replace_allocator){/*do nothing*/}
	
	template<unsigned int IRank>
	void check_initializer_shape(const typename detail::initializer_for_rank<value_type,IRank>::type& init, unsigned int entry){
		if(init.size()!=extents[Rank-IRank])
			throw std::logic_error("Incorrect shape for multidimensional array initializer:\nentry "+
								   std::to_string(entry)+" in dimension "+std::to_string(Rank-IRank)+
								   " has size "+std::to_string(init.size())+" but should have size "+
								   std::to_string(extents[Rank-IRank]));
		if(IRank>1){
			unsigned int subentry=0;
			for(const auto& i : init)
				check_initializer_shape<IRank-1>(i,subentry++);
		}
	}
	
	template<unsigned int IRank>
	void initialize_multidim(const typename detail::initializer_for_rank<value_type,IRank>::type& init, value_type*& pos){
		if(IRank!=1){
			for(const auto& i : init)
				initialize_multidim(i,pos);
		}
		else{
			for(const auto& i : init)
				allocator_traits::construct(allocator(),pos++,i);
		}
	}
	
public:
	template<typename DerefType>
	struct iterator_type : public iterator_base<DerefType,iterator_type<DerefType>>{
	private:
		using base_type=iterator_base<DerefType,iterator_type<DerefType>>;
		template<typename DummyType> struct dummy{typedef DummyType type;};
		friend typename dummy<base_type>::type;
		friend class marray;
	public:
		using base_type::base_type;
		using difference_type=typename base_type::difference_type;
		using pointer=typename base_type::pointer;
		difference_type operator-(const iterator_type& other) const{
			return(this->ptr-other.ptr);
		}
		bool operator<(const iterator_type& other) const{
			return(this->ptr<other.ptr);
		}
		bool operator>(const iterator_type& other) const{
			return(this->ptr>other.ptr);
		}
		bool operator<=(const iterator_type& other) const{
			return(this->ptr<=other.ptr);
		}
		bool operator>=(const iterator_type& other) const{
			return(this->ptr>=other.ptr);
		}
	private:
		pointer increment(difference_type n){ return(this->ptr+n); }
		pointer decrement(difference_type n){ return(this->ptr-n); }
	};
	
	template<typename DerefType>
	struct reverse_iterator_type : public iterator_base<DerefType,reverse_iterator_type<DerefType>>{
	private:
		using base_type=iterator_base<DerefType,reverse_iterator_type<DerefType>>;
		template<typename DummyType> struct dummy{typedef DummyType type;};
		friend typename dummy<base_type>::type;
		friend class marray;
	public:
		using base_type::base_type;
		using difference_type=typename base_type::difference_type;
		using pointer=typename base_type::pointer;
		difference_type operator-(const reverse_iterator_type& other) const{
			return(other.ptr-this->ptr);
		}
		bool operator<(const reverse_iterator_type& other) const{
			return(this->ptr>other.ptr);
		}
		bool operator>(const reverse_iterator_type& other) const{
			return(this->ptr<other.ptr);
		}
		bool operator<=(const reverse_iterator_type& other) const{
			return(this->ptr>=other.ptr);
		}
		bool operator>=(const reverse_iterator_type& other) const{
			return(this->ptr<=other.ptr);
		}
	private:
		pointer increment(difference_type n){ return(this->ptr-n); }
		pointer decrement(difference_type n){ return(this->ptr+n); }
	};
	
	using iterator=iterator_type<value_type>;
	using const_iterator=iterator_type<const value_type>;
	using reverse_iterator=reverse_iterator_type<value_type>;
	using const_reverse_iterator=reverse_iterator_type<const value_type>;
	
	using subscript_traits=typename detail::subscript_traits<value_type,Rank>;
	using subscript_type=typename subscript_traits::result_type;
	using const_subscript_traits=typename detail::subscript_traits<const value_type,Rank>;
	using const_subscript_type=typename const_subscript_traits::result_type;
	
	using difference_type = typename iterator::difference_type;
	
	///Construct an empty marray
	///\param alloc the allocator to be used by the array to obtain memory
	explicit marray(allocator_type alloc=allocator_type()):
	total_size_and_alloc(0,alloc),
	data(nullptr){
		std::fill(extents.begin(),extents.end(),0);
	}
	
	///Construct an marray with a given size, whose contents will be default constructed
	///\param dims the extents of the array in all of its dimensions
	///\param alloc the allocator to be used by the array to obtain memory
	explicit marray(std::initializer_list<size_type> dims, allocator_type alloc=allocator_type()):
	total_size_and_alloc(compute_total_size(dims),alloc),
	data(total_size()?allocator_traits::allocate(allocator(),total_size()):nullptr){
		std::copy(dims.begin(),dims.end(),extents.begin());
		size_type i=0;
		try{
			for(size_type end=total_size(); i!=end; i++)
				allocator_traits::construct(allocator(),data+i);
		}catch(...){
			//if something went wrong destroy all of the already constructed objects
			for(size_type j=0; j!=i; j++)
				allocator_traits::destroy(alloc,data+j);
			allocator_traits::deallocate(allocator(),data,total_size());
			throw;
		}
	}
	
	///Construct an marray with a given size, whose contents will be default constructed
	///\param dims the extents of the array in all of its dimensions
	///\param alloc the allocator to be used by the array to obtain memory
	template<typename IndexContainer,
	         typename=typename std::enable_if<detail::is_index_container<IndexContainer>::value>::type>
	explicit marray(const IndexContainer& dims, allocator_type alloc=allocator_type()):
	total_size_and_alloc(compute_total_size(dims),alloc),
	data(total_size()?allocator_traits::allocate(allocator(),total_size()):nullptr){
		std::copy(dims.begin(),dims.end(),extents.begin());
		size_type i=0;
		try{
			for(size_type end=total_size(); i!=end; i++)
				allocator_traits::construct(allocator(),data+i);
		}catch(...){
			//if something went wrong destroy all of the already constructed objects
			for(size_type j=0; j!=i; j++)
				allocator_traits::destroy(alloc,data+j);
			allocator_traits::deallocate(allocator(),data,total_size());
			throw;
		}
	}
	
	///Construct an marray with the given size, leaving the contents unitinialized.
	///Naturally this is only valid for types which do not require initialization,
	///i.e. trivial types.
	///\param dims the extents of the array in all of its dimensions
	///\param alloc the allocator to be used by the array to obtain memory
	explicit marray(std::initializer_list<size_type> dims, const detail::no_init_tag&, allocator_type alloc=allocator_type()):
	total_size_and_alloc(compute_total_size(dims),alloc),
	data(total_size()?allocator_traits::allocate(allocator(),total_size()):nullptr){
		std::copy(dims.begin(),dims.end(),extents.begin());
		static_assert(std::is_trivial<T>::value,"Only containers of trivial types can be constructed without initialization");
	}
	
	///Construct an marray with the given size, leaving the contents unitinialized.
	///Naturally this is only valid for types which do not require initialization,
	///i.e. trivial types.
	///\param dims the extents of the array in all of its dimensions
	///\param alloc the allocator to be used by the array to obtain memory
	template<typename IndexContainer,
	typename=typename std::enable_if<detail::is_index_container<IndexContainer>::value>::type>
	explicit marray(const IndexContainer& dims, const detail::no_init_tag&, allocator_type alloc=allocator_type()):
	extents(dims),
	total_size_and_alloc(compute_total_size(dims),alloc),
	data(total_size()?allocator_traits::allocate(allocator(),total_size()):nullptr){
		static_assert(std::is_trivial<T>::value,"Only containers of trivial types can be constructed without initialization");
	}
	
	///Construct an marray with a given size and initialize its contents to given values
	///\param dims the extents of the array in all of its dimensions
	///\param init a nested set of initializer lists which define the contents of the array
	///\param alloc the allocator to be used by the array to obtain memory
	marray(std::initializer_list<size_type> dims, typename detail::initializer_for_rank<value_type,Rank>::type init, allocator_type alloc=allocator_type()):
	total_size_and_alloc(compute_total_size(dims),alloc),
	data(total_size()?allocator_traits::allocate(allocator(),total_size()):nullptr){
		std::copy(dims.begin(),dims.end(),extents.begin());
		
		//ensure that the initializer has the correct shape
		detail::initializer_shape_checker<value_type,Rank,Rank>::check(extents,init,0);
		
		pointer write_pos=data;
		try{
			detail::multidimensional_initializer<value_type,Rank,allocator_type,Rank>::initialize(init,allocator(),write_pos);
		}catch(...){
			for(pointer p=data; p!=write_pos; p++)
				allocator_traits::destroy(allocator(),p);
			allocator_traits::deallocate(allocator(),data,total_size());
			throw;
		}
	}
	
	///Construct an marray with a given size initialize its contents to values
	///from another container. Since the input iterator pair is interpreted as a
	///one dimensional range this is equivalent to constructing the marray without
	///initialization and then using std::copy to copy the data into it; that is
	///the values will be stored in iteration order.
	///\param dims the extents of the array in all of its dimensions
	///\param init_begin iterator to the start of the data to be copied into the array
	///\param init_end iterator to the start of the data to be copied into the array
	///\param alloc the allocator to be used by the array to obtain memory
	///\pre std::distance(init_begin,init_end) must equal the product of the entries in dims
	template<typename ForwardIterator>
	marray(std::initializer_list<size_type> dims, ForwardIterator init_begin, ForwardIterator init_end, allocator_type alloc=allocator_type()):
	total_size_and_alloc(compute_total_size(dims),alloc),
	data(total_size()?allocator_traits::allocate(allocator(),total_size()):nullptr){
		std::copy(dims.begin(),dims.end(),extents.begin());
		
		if(std::distance(init_begin,init_end)!=total_size())
			throw std::logic_error("Iterator range has wrong size to initialize marray");
		
		pointer write_pos=data;
		try{
			for(; init_begin!=init_end; init_begin++, write_pos++)
				*write_pos = *init_begin;
		}catch(...){
			for(pointer p=data; p!=write_pos; p++)
				allocator_traits::destroy(allocator(),p);
			allocator_traits::deallocate(allocator(),data,total_size());
			throw;
		}
	}
	
	///Construct an marray by copying from an existing instance
	marray(const marray& other):
	extents(other.extents),
	total_size_and_alloc(other.total_size(),allocator_traits::select_on_container_copy_construction(other.allocator())),
	data(total_size()?allocator_traits::allocate(allocator(),total_size()):nullptr){
		try{
			copy_init_buffer(data,total_size(),allocator(),other.data);
		}catch(...){
			allocator_traits::deallocate(allocator(),data,total_size());
			throw;
		}
	}
	
	///Construct an marray by moving from an existing instance
	marray(marray&& other):
	extents(other.extents),
	total_size_and_alloc(other.total_size(),detail::select_on_container_move_construction(other.allocator())),
	data(other.data){
		std::fill_n(other.extents.begin(),Rank,0);
		other.total_size()=0;
		other.data=nullptr;
	}
	
	~marray() noexcept{
		if(data)
			destroy_buffer(data,total_size(),allocator());
	}
	
	marray& operator=(const marray& other){
		if(&other==this)
			return(*this);
		copy_assign(other,typename allocator_traits::propagate_on_container_copy_assignment());
		std::copy_n(other.extents.begin(),Rank,extents.begin());
		total_size()=other.total_size();
		return(*this);
	}
	
	marray& operator=(marray&& other)
	noexcept(allocator_traits::propagate_on_container_move_assignment::value && std::is_nothrow_move_assignable<allocator_type>::value)
	{
		if(&other==this)
			return(*this);
		move_assign(std::move(other),typename allocator_traits::propagate_on_container_move_assignment());
		return(*this);
	}
	
	///Assignment from a multidimensional initializer
	///\todo: make this exception safe
	///\todo: make this work when the initializer shape does not match the array shape
	///\todo: when the shapes do match and value_type is nothrow-copy-assignable do not reallocate
	marray& operator=(typename detail::initializer_for_rank<value_type,Rank>::type init){
		//ensure that the initializer has the same shape
		detail::initializer_shape_checker<value_type,Rank,Rank>::check(extents,init,0);
		
		pointer write_pos=data;
		detail::multidimensional_initializer<value_type,Rank,allocator_type,Rank>::assign(init,write_pos);
		return(*this);
	}
	
	subscript_type operator[](size_type i){
		return(subscript_traits::make_subscript(data,extents,i));
	}
	
	const_subscript_type operator[](size_type i) const{
		return(const_subscript_traits::make_subscript(data,extents,i));
	}
	
	
	
	const_reference operator[](std::initializer_list<size_t> indices) const{
		assert(indices.size()==Rank);
		size_t offset=0;
		//Since Rank>0 this is always a valid iterator, although it may be extents.end()
		auto extentIt=extents.begin()+1, extentEnd=extents.end();
		size_t stride=std::accumulate(extentIt,extentEnd,1,std::multiplies<size_type>());
		for(size_t index : indices){
			offset+=index*stride;
			if(extentIt!=extentEnd){ //avoid undefined behavior on last iteration
				stride/=*extentIt;
				extentIt++;
			}
		}
		return(*(data+offset));
	}
	
	reference operator[](std::initializer_list<size_t> indices){
		assert(indices.size()==Rank);
		size_t offset=0;
		//Since Rank>0 this is always a valid iterator, although it may be extents.end()
		auto extentIt=extents.begin()+1, extentEnd=extents.end();
		size_t stride=std::accumulate(extentIt,extentEnd,1,std::multiplies<size_type>());
		for(size_t index : indices){
			offset+=index*stride;
			if(extentIt!=extentEnd){ //avoid undefined behavior on last iteration
				stride/=*extentIt;
				extentIt++;
			}
		}
		return(*(data+offset));
	}
	
	template<typename Container,
	typename=typename std::enable_if<detail::is_index_container<Container>::value>::type>
	reference operator[](Container indices){
		assert(indices.size()==Rank);
		size_t offset=0;
		//Since Rank>0 this is always a valid iterator, although it may be extents.end()
		auto extentIt=extents.begin()+1, extentEnd=extents.end();
		size_t stride=std::accumulate(extentIt,extentEnd,1,std::multiplies<size_type>());
		for(size_t index : indices){
			offset+=index*stride;
			if(extentIt!=extentEnd){ //avoid undefined behavior on last iteration
				stride/=*extentIt;
				extentIt++;
			}
		}
		return(*(data+offset));
	}
	
	template<typename Container,
	typename=typename std::enable_if<detail::is_index_container<Container>::value>::type>
	const_reference operator[](Container indices) const{
		assert(indices.size()==Rank);
		size_t offset=0;
		//Since Rank>0 this is always a valid iterator, although it may be extents.end()
		auto extentIt=extents.begin()+1, extentEnd=extents.end();
		size_t stride=std::accumulate(extentIt,extentEnd,1,std::multiplies<size_type>());
		for(size_t index : indices){
			offset+=index*stride;
			if(extentIt!=extentEnd){ //avoid undefined behavior on last iteration
				stride/=*extentIt;
				extentIt++;
			}
		}
		return(*(data+offset));
	}
	
	//TODO: noexcept condition is wrong if swap(allocator_type&,allocator_type&) is found via ADL
	///Swap contents with another array
	void swap(marray& other)
	noexcept(!allocator_traits::propagate_on_container_swap::value || noexcept(std::swap(std::declval<allocator_type>(),std::declval<allocator_type>())))
	{
		using std::swap;
		for(unsigned int i=0; i<Rank; i++)
			swap(extents[i],other.extents[i]);
		swap(total_size(),other.total_size());
		swap_allocator(other.allocator(),allocator_traits::propagate_on_container_swap());
		swap(data,other.data);
	}
	
	///Change the size of the array along one dimension
	///\param dim the dimension in which to resize
	///\param size the new size of the array in the given dimension
	void resize(size_type dim, size_type size){
		assert(dim<Rank);
		if(size<extents[dim])
			erase(dim,size,extents[dim]-size);
		else if(size>extents[dim])
			insert(dim,extents[dim],size-extents[dim]);
	}
	
	template<typename Container>
	void resize(const Container& new_extents){
		if(std::distance(std::begin(new_extents),std::end(new_extents))!=Rank)
			throw std::logic_error("Incorrect number of dimensions when attempting to resize multidimensional array");
		size_type dim=0;
		for(size_type size : new_extents)
			resize(dim++,size);
	}
	
private:
	void insert(size_type dim, size_type pos, size_t num, value_type val, std::true_type move){
		const auto mult=std::multiplies<size_type>();
		const size_type iterations=std::accumulate(extents.begin(),extents.begin()+dim,1u,mult);
		const size_type remainder=extents[dim]-pos;
		const size_type chunk_size=std::accumulate(extents.begin()+dim+1,extents.end(),1u,mult);
		
		const size_type new_total_size=iterations*(extents[dim]+num)*chunk_size;
		value_type* new_data=allocator_traits::allocate(allocator(),new_total_size); //may throw
		
		value_type* write_pos=new_data;
		value_type* read_pos=data;
		size_type i=0, new_j=0;
		try{
			for(i=0; i<iterations; i++){
				//move a segment of old items
				for(size_type j=0; j!=pos*chunk_size; j++,write_pos++,read_pos++)
					allocator_traits::construct(allocator(),write_pos,std::move(*read_pos));
				//insert a segment of new items
				for(new_j=0; new_j!=num*chunk_size; new_j++,write_pos++)
					allocator_traits::construct(allocator(),write_pos,val);
				//move another segment of old items
				for(size_type j=0; j!=remainder*chunk_size; j++,write_pos++,read_pos++)
					allocator_traits::construct(allocator(),write_pos,std::move(*read_pos));
			}
		}catch(...){
			//need to move back each moved object, and destroy all objects which were constructed
			//Note that since move construction cannot throw, we can unconditionally process whole
			//stretches of moved objects; only the stretch of copy constructed objects might have
			//been broken.
			write_pos=new_data;
			read_pos=data;
			for(size_type i_=0; i_<=i; i_++){
				//put back a segment of moved items
				for(size_type j_=0; j_!=pos*chunk_size; j_++,write_pos++,read_pos++){
					*read_pos=std::move(*write_pos);
					allocator_traits::destroy(allocator(),write_pos);
				}
				//destroy a segment of newly constructed items
				for(size_type j_=0, end=(i_==i?new_j:num*chunk_size); j_<end; j_++,write_pos++)
					allocator_traits::destroy(allocator(),write_pos);
				if(i_==i)
					break;
				//put back another segment of moved items
				for(size_type j_=0; j_!=remainder*chunk_size; j_++,write_pos++,read_pos++){
					*read_pos=std::move(*write_pos);
					allocator_traits::destroy(allocator(),write_pos);
				}
			}
			allocator_traits::deallocate(allocator(),new_data,new_total_size);
			throw;
		}
		destroy_buffer(data, total_size(), allocator());
		data=new_data;
		total_size()=new_total_size;
		extents[dim]+=num;
	}
	
	void insert(size_type dim, size_type pos, size_t num, value_type val, std::false_type dont_move){
		const auto mult=std::multiplies<size_type>();
		const size_type iterations=std::accumulate(extents.begin(),extents.begin()+dim,1u,mult);
		const size_type remainder=extents[dim]-pos;
		const size_type chunk_size=std::accumulate(extents.begin()+dim+1,extents.end(),1u,mult);
		
		const size_type new_total_size=iterations*(extents[dim]+num)*chunk_size;
		value_type* new_data=allocator_traits::allocate(allocator(),new_total_size); //may throw
		
		value_type* write_pos=new_data;
		value_type* read_pos=data;
		try{
			for(size_type i=0; i<iterations; i++){
				//copy a segment of old items
				for(size_type j=0; j!=pos*chunk_size; j++,write_pos++,read_pos++)
					allocator_traits::construct(allocator(),write_pos,*read_pos);
				//insert a segment of new items
				for(size_type new_j=0; new_j!=num*chunk_size; new_j++,write_pos++)
					allocator_traits::construct(allocator(),write_pos,val);
				//copy another segment of old items
				for(size_type j=0; j!=remainder*chunk_size; j++,write_pos++,read_pos++)
					allocator_traits::construct(allocator(),write_pos,*read_pos);
			}
		}catch(...){
			//destroy all constructed objects
			for(value_type* ptr=new_data; ptr!=write_pos; ptr++)
				allocator_traits::destroy(allocator(),ptr);
			allocator_traits::deallocate(allocator(),new_data,new_total_size);
			throw;
		}
		
		destroy_buffer(data, total_size(), allocator());
		data=new_data;
		total_size()=new_total_size;
		extents[dim]+=num;
	}

public:
	///insert size entries before index pos in dimension dim
	///\param dim the dimension in which to expand the array
	///\param pos the index in the given dimension before which to insert the new entries
	///\param num the number of entries to add in the given dimension
	///\param val the initial value to assign to the added entries
	void insert(size_type dim, size_type pos, size_t num, value_type val=T()){
		assert(dim<Rank);
		assert(pos<=extents[dim]);
		
		//choose to use move operations only if we can be certain they won't throw
		constexpr bool move=std::is_nothrow_move_constructible<value_type>::value
		                    && std::is_nothrow_move_assignable<value_type>::value;
		insert(dim,pos,num,val,std::integral_constant<bool,move>());
	}

private:
	void erase(size_type dim, size_type pos, size_t num, std::true_type move){
		const auto mult=std::multiplies<size_type>();
		const size_type iterations=std::accumulate(extents.begin(),extents.begin()+dim,1u,mult);
		const size_type remainder=extents[dim]-(pos+num);
		const size_type chunk_size=std::accumulate(extents.begin()+dim+1,extents.end(),1u,mult);
		
		const size_type new_total_size=iterations*(extents[dim]-num)*chunk_size;
		value_type* new_data=allocator_traits::allocate(allocator(),new_total_size); //may throw
		
		value_type* write_pos=new_data;
		value_type* read_pos=data;
		
		for(unsigned int i=0; i<iterations; i++){
			for(size_type j=0; j!=pos*chunk_size; j++,write_pos++,read_pos++)
				allocator_traits::construct(allocator(),write_pos,std::move(*read_pos));
			read_pos+=num*chunk_size; //skip the items being erased
			for(size_type j=0; j!=remainder*chunk_size; j++,write_pos++,read_pos++)
				allocator_traits::construct(allocator(),write_pos,std::move(*read_pos));
		}
		
		destroy_buffer(data, total_size(), allocator());
		data=new_data;
		total_size()=new_total_size;
		extents[dim]-=num;
	}
	
	void erase(size_type dim, size_type pos, size_t num, std::false_type dont_move){
		const auto mult=std::multiplies<size_type>();
		const size_type iterations=std::accumulate(extents.begin(),extents.begin()+dim,1u,mult);
		const size_type remainder=extents[dim]-(pos+num);
		const size_type chunk_size=std::accumulate(extents.begin()+dim+1,extents.end(),1u,mult);
		
		const size_type new_total_size=iterations*(extents[dim]-num)*chunk_size;
		value_type* new_data=allocator_traits::allocate(allocator(),new_total_size); //may throw
		
		value_type* write_pos=new_data;
		value_type* read_pos=data;
		
		try{
			for(unsigned int i=0; i<iterations; i++){
				for(size_type j=0; j!=pos*chunk_size; j++,write_pos++,read_pos++)
					allocator_traits::construct(allocator(),write_pos,*read_pos);
				read_pos+=num*chunk_size; //skip the items being erased
				for(size_type j=0; j!=remainder*chunk_size; j++,write_pos++,read_pos++)
					allocator_traits::construct(allocator(),write_pos,*read_pos);
			}
		}catch(...){
			//just need to destroy the objects which were successfully copied
			for(value_type* ptr=new_data; ptr!=write_pos; ptr++)
				allocator_traits::destroy(allocator(),ptr);
			allocator_traits::deallocate(allocator(),new_data,new_total_size);
			throw;
		}
		
		destroy_buffer(data, total_size(), allocator());
		data=new_data;
		total_size()=new_total_size;
		extents[dim]-=num;
	}
	
public:
	///remove size entries beginning index pos in dimension dim
	///\param dim the dimension in which to contract the array
	///\param pos the index in the given dimension at which to begin removing entries
	///\param num the number of entries to remove in the given dimension
	void erase(size_type dim, size_type pos, size_t num){
		assert(dim<Rank);
		assert(pos<extents[dim]);
		assert(pos+num<=extents[dim]);
		
		const auto mult=std::multiplies<size_type>();
		const size_type iterations=std::accumulate(extents.begin(),extents.begin()+dim,1u,mult);
		const size_type remainder=extents[dim]-(pos+num);
		const size_type chunk_size=std::accumulate(extents.begin()+dim+1,extents.end(),1u,mult);
		
		const size_type new_total_size=iterations*(extents[dim]-num)*chunk_size;
		value_type* new_data=allocator_traits::allocate(allocator(),new_total_size); //may throw
		
		value_type* write_pos=new_data;
		value_type* read_pos=data;
		//choose to use move operations only if we can be certain they won't throw
		constexpr bool move=std::is_nothrow_move_constructible<value_type>::value
		  && std::is_nothrow_move_assignable<value_type>::value;
		try{
			for(unsigned int i=0; i<iterations; i++){
				for(size_type j=0; j!=pos*chunk_size; j++,write_pos++,read_pos++)
					allocator_traits::construct(allocator(),write_pos,(move?std::move(*read_pos):*read_pos));
				read_pos+=num*chunk_size; //skip the items being erased
				for(size_type j=0; j!=remainder*chunk_size; j++,write_pos++,read_pos++)
					allocator_traits::construct(allocator(),write_pos,(move?std::move(*read_pos):*read_pos));
			}
		}catch(...){
			//Can only get here if we were copying, in which case we just need to
			//destroy the objects which were successfully copied
			for(value_type* ptr=new_data; ptr!=write_pos; ptr++)
				allocator_traits::destroy(allocator(),ptr);
			allocator_traits::deallocate(allocator(),new_data,new_total_size);
			throw;
		}
		
		destroy_buffer(data, total_size(), allocator());
		data=new_data;
		total_size()=new_total_size;
		extents[dim]-=num;
	}
	
	///Reshape the array to potentially different sizes in each dimension
	///as long as the total size (the product of extents in all dimensions)
	///remains the same
	///\param new_extents an iterable object whose contents will be used as the new sizes
	template<typename Container>
	void reshape(const Container& new_extents){
		if(std::distance(std::begin(new_extents),std::end(new_extents))!=Rank)
			throw std::logic_error("Incorrect number of dimensions when attempting to reshape multidimensional array");
		size_type new_total_size=std::accumulate(new_extents.begin(),new_extents.end(),1u,std::multiplies<size_type>());
		if(new_total_size!=total_size())
			throw std::logic_error("Inconsistent total size when attempting to reshape multidimensional array");
		std::copy_n(std::begin(new_extents),Rank,extents.begin());
	}
	
		///\brief Compute the sum of all entries in the array.
	///
	///This function has no advantage over a the naive approach; it exists for
	///uniformity with the more advanced overloads.
	///\return The scalar sum of the array
	value_type sum() const{
		return(std::accumulate(begin(),end(),value_type(0),std::plus<value_type>()));
	}
	
	//TODO: consider allocator of result array
	///\brief Produce a new array by summing entries along one dimension
	///\param dim The index of the dimension over which to sum.
	///\return An marray with rank reduced by one.
	marray<value_type,Rank-1> sum(size_t dim) const{
		static_assert(Rank>1, "Single dimension sum cannot produce a rank 0 result; use the nullary sum()");
		//copy all extents except the one projected out
		std::array<size_type,Rank-1> sExtents;
		for(size_t i=0,j=0; i<Rank; i++){
			if(i!=dim)
				sExtents[j++]=extents[i];
		}
		marray<value_type,Rank-1> s(sExtents);
		std::fill(s.begin(),s.end(),value_type(0));
		std::array<size_type,Rank-1> indices;
		std::fill(indices.begin(),indices.end(),0);
		size_t dummy_index=0;
		for(auto it=begin(), endIt=end(); it!=endIt; it++){
			s[indices]+=*it;
			//this is messy
			for(size_t i=Rank-1,j=Rank-2; ; i--){
				if(i!=dim){
					if(indices[j]<extents[i]-1){
						indices[j--]++;
						break;
					}
					indices[j--]=0;
				}
				else{
					if(dummy_index<extents[i]-1){
						dummy_index++;
						break;
					}
					dummy_index=0;
				}
				if(i==0) break;
			}
		}
		
		return(s);
	}
	
	//TODO: consider allocator of result array
	///\brief Produce a new array by summing entries along one dimension, but
	///       with the same rank as the original
	///\param dim The index of the dimension over which to sum.
	///\return An marray with the extent in the summed dimension reduced to one.
	marray<value_type,Rank> sum_same_rank(size_t dim) const{
		//copy all extents except the one projected out
		std::array<size_type,Rank> sExtents;
		for(size_t i=0,j=0; i<Rank; i++,j++){
			if(i!=dim)
				sExtents[j]=extents[i];
			else //reduce the projected extent to one
				sExtents[j]=1;
		}
		marray<value_type,Rank> s(sExtents);
		std::fill(s.begin(),s.end(),value_type(0));
		std::array<size_type,Rank> indices;
		std::fill(indices.begin(),indices.end(),0);
		size_t dummy_index=0;
		for(auto it=begin(), endIt=end(); it!=endIt; it++){
			s[indices]+=*it;
			//this is messy
			for(size_t i=Rank-1; ; i--){
				if(i!=dim){
					if(indices[i]<extents[i]-1){
						indices[i]++;
						break;
					}
					indices[i]=0;
				}
				else{
					if(dummy_index<extents[i]-1){
						dummy_index++;
						break;
					}
					dummy_index=0;
				}
				if(i==0) break;
			}
		}
		
		return(s);
	}
	
	iterator begin(){ return(iterator(this,0)); }
	const_iterator begin() const{ return(const_iterator(this,0)); }
	const_iterator cbegin() const{ return(const_iterator(this,0)); }
	iterator end(){ return(iterator(this,total_size())); }
	const_iterator end() const{ return(const_iterator(this,total_size())); }
	const_iterator cend() const{ return(const_iterator(this,total_size())); }
	reverse_iterator rbegin(){ return(reverse_iterator(this,total_size()-1)); }
	const_reverse_iterator rbegin() const{ return(const_reverse_iterator(this,total_size()-1)); }
	const_reverse_iterator crbegin() const{ return(const_reverse_iterator(this,total_size()-1)); }
	reverse_iterator rend(){ return(reverse_iterator(this,-1)); }
	const_reverse_iterator rend() const{ return(const_reverse_iterator(this,-1)); }
	const_reverse_iterator crend() const{ return(const_reverse_iterator(this,-1)); }
	
	reference front(){ return(data[0]); }
	const_reference front() const { return(data[0]); }
	reference back(){ return(data[size()-1]); }
	const_reference back() const { return(size()-1); }
	
	///Get the allocator currently used by the array
	allocator_type get_allocator() const{ return(allocator()); }
	
	///Get the rank (number of dimensions) of the array
	constexpr size_type rank() const noexcept{ return(Rank); }
	///Get the size of the array in the given dimension
	///\param dim the dimension to query
	size_type extent(size_type dim) const{ return(extents[dim]); }
	const std::array<size_type,Rank> get_extents() const{ return(extents); }
	///Get the the total size of the array (the product of sizes in all dimensions)
	size_type size() const noexcept{ return(total_size()); }
	///Check whether the array is empty (has zero size)
	bool empty() const noexcept{ return(total_size()==0); }
	const pointer get_data() const noexcept{ return(data); }
	
	template<typename DerefType, typename DerivedType>
	size_type getCoordinate(iterator_base<DerefType,DerivedType> it, size_type dim){
		size_t offset=it.ptr-data;
		size_type divisor=total_size();
		for(unsigned int d=0; d<dim; d++)
			offset%=(divisor/=extents[d]);
		offset/=(divisor/extents[dim]);
		return(offset);
	}
private:
	std::array<size_type,Rank> extents;
	using size_and_alloc_type=detail::compressed_pair<size_type,allocator_type>;
	size_and_alloc_type total_size_and_alloc;
	value_type* data;
	
	///Verify that an initializer list contains the correct number of sizes
	///and multiply them all together to get the total represented size
	///\param exts the putative list of sizes in all dimensions
	template<typename IndexContainer>
	static size_type compute_total_size(const IndexContainer& exts){
		if(exts.size()!=Rank)
			throw std::logic_error("Incorrect number of dimensions used to initialize multidimensional array");
		return(std::accumulate(exts.begin(),exts.end(),1UL,std::multiplies<size_type>()));
	}
	
	///Get the number of objects which exist in the array
	size_type& total_size() noexcept{ return(total_size_and_alloc.first()); }
	///Get the number of objects which exist in the array
	const size_type& total_size() const noexcept{ return(total_size_and_alloc.first()); }
	///Get the allocator currently used by the array
	allocator_type& allocator() noexcept{ return(total_size_and_alloc.second()); }
	///Get the allocator currently used by the array
	const allocator_type& allocator() const noexcept{ return(total_size_and_alloc.second()); }
	
	///Invoke the destructor of each object in a buffer and deallocate the buffer.
	///Must not throw.
	///\param buffer the buffer to be deallocated
	///\param size the size of the buffer (in objects)
	///\param alloc the allocator to use
	static void destroy_buffer(value_type* buffer, size_type size, allocator_type& alloc) noexcept{
		for(size_type i=0; i!=size; i++)
			allocator_traits::destroy(alloc,buffer+i);
		allocator_traits::deallocate(alloc,buffer,size);
	}
	
	///Initialize a buffer by copy constructing objects from corresponding objects in another buffer
	///
	///May throw if value_type(const value_type&) throws
	///\param buffer start of the buffer to be initialized
	///\param size the size of the buffer (in objects)
	///\param alloc the allocator to use
	///\param source_buffer the existing buffer from whose objects thes objects in the buffer are copy constructed
	///\pre source_buffer is a buffer of at least size valid value_type objects
	static void copy_init_buffer(value_type* buffer, size_type size, allocator_type& alloc, value_type* source_buffer)
	noexcept(noexcept(allocator_traits::construct(std::declval<allocator_type&>(),std::declval<const value_type*>(),std::declval<const value_type&>()))){
		size_type i=0;
		try{
			for(; i!=size; i++)
				allocator_traits::construct(alloc,buffer+i,*(source_buffer+i));
		}catch(...){
			//destroy any objects which had been successfully constructed
			for(size_type j=0; j!=i; j++)
				allocator_traits::destroy(alloc,buffer+j);
			throw;
		}
	}
};

template<typename Iterator, int dummy=Iterator::marray_iterator_tag>
Iterator operator+(typename std::iterator_traits<Iterator>::difference_type n, const Iterator& it){
	return(it+n);
}

//----------------------------------
// Arithmetic operations on marrays
//----------------------------------

//unary negation
template<typename T, unsigned Rank>
marray<T,Rank> operator-(marray<T,Rank>&& a){
	for(auto it=a.begin(), end=a.end(); it!=end; it++)
		*it=-*it;
	return(a);
}
template<typename T, unsigned Rank>
marray<T,Rank> operator-(const marray<T,Rank>& a){
	return(-marray<T,Rank>(a));
}

//Assignment forms are only relevant when types and ranks match

template<typename T, unsigned Rank>
marray<T,Rank>& operator+=(marray<T,Rank>& a1, const marray<T,Rank>& a2){
	assert(a1.size()==a2.size());
	auto i1=a1.begin(), end=a1.end();
	for(auto i2=a2.begin(); i1!=end; i1++,i2++)
		*i1+=*i2;
	return(a1);
}
template<typename T, unsigned Rank>
marray<T,Rank>&& operator+=(marray<T,Rank>&& a1, const marray<T,Rank>& a2){
	assert(a1.size()==a2.size());
	auto i1=a1.begin(), end=a1.end();
	for(auto i2=a2.begin(); i1!=end; i1++,i2++)
		*i1+=*i2;
	return(std::move(a1));
}

template<typename T, unsigned Rank>
marray<T,Rank>& operator-=(marray<T,Rank>& a1, const marray<T,Rank>& a2){
	assert(a1.size()==a2.size());
	auto i1=a1.begin(), end=a1.end();
	for(auto i2=a2.begin(); i1!=end; i1++,i2++)
		*i1-=*i2;
	return(a1);
}
template<typename T, unsigned Rank>
marray<T,Rank>&& operator-=(marray<T,Rank>&& a1, const marray<T,Rank>& a2){
	assert(a1.size()==a2.size());
	auto i1=a1.begin(), end=a1.end();
	for(auto i2=a2.begin(); i1!=end; i1++,i2++)
		*i1-=*i2;
	return(std::move(a1));
}

template<typename T, unsigned Rank>
marray<T,Rank>& operator*=(marray<T,Rank>& a1, const marray<T,Rank>& a2){
	assert(a1.size()==a2.size());
	auto i1=a1.begin(), end=a1.end();
	for(auto i2=a2.begin(); i1!=end; i1++,i2++)
		*i1*=*i2;
	return(a1);
}
template<typename T, unsigned Rank>
marray<T,Rank>&& operator*=(marray<T,Rank>&& a1, const marray<T,Rank>& a2){
	assert(a1.size()==a2.size());
	auto i1=a1.begin(), end=a1.end();
	for(auto i2=a2.begin(); i1!=end; i1++,i2++)
		*i1*=*i2;
	return(std::move(a1));
}

template<typename T, unsigned Rank>
marray<T,Rank>& operator/=(marray<T,Rank>& a1, const marray<T,Rank>& a2){
	assert(a1.size()==a2.size());
	auto i1=a1.begin(), end=a1.end();
	for(auto i2=a2.begin(); i1!=end; i1++,i2++)
		*i1/=*i2;
	return(a1);
}
template<typename T, unsigned Rank>
marray<T,Rank>&& operator/=(marray<T,Rank>&& a1, const marray<T,Rank>& a2){
	assert(a1.size()==a2.size());
	auto i1=a1.begin(), end=a1.end();
	for(auto i2=a2.begin(); i1!=end; i1++,i2++)
		*i1/=*i2;
	return(std::move(a1));
}

//forms with scalars
//addition
template<typename T, unsigned Rank, typename U>
marray<T,Rank>& operator+=(marray<T,Rank>& a, const U& f){
	for(auto& i : a) i+=f;
	return(a);
}
template<typename T, unsigned Rank, typename U>
marray<T,Rank>&& operator+=(marray<T,Rank>&& a, const U& f){
	for(auto& i : a) i+=f;
	return(std::move(a));
}

template<typename T, unsigned Rank, typename U,
typename=typename std::enable_if<!detail::is_marray<U>::value>::type>
marray<T,Rank> operator+(marray<T,Rank>& a, const U& f){
	return(marray<T,Rank>(a)+=f);
}
template<typename T, unsigned Rank, typename U,
typename=typename std::enable_if<!detail::is_marray<U>::value>::type>
marray<T,Rank> operator+(const U& f, marray<T,Rank>& a){
	return(marray<T,Rank>(a)+=f);
}

//subtraction
template<typename T, unsigned Rank, typename U>
marray<T,Rank>& operator-=(marray<T,Rank>& a, const U& f){
	for(auto& i : a) i-=f;
	return(a);
}
template<typename T, unsigned Rank, typename U>
marray<T,Rank>&& operator-=(marray<T,Rank>&& a, const U& f){
	for(auto& i : a) i-=f;
	return(std::move(a));
}

template<typename T, unsigned Rank, typename U,
typename=typename std::enable_if<!detail::is_marray<U>::value>::type>
marray<T,Rank> operator-(marray<T,Rank>& a, const U& f){
	return(marray<T,Rank>(a)-=f);
}
template<typename T, unsigned Rank, typename U,
typename=typename std::enable_if<!detail::is_marray<U>::value>::type>
marray<T,Rank> operator-(const U& f, marray<T,Rank>& a){
	return(marray<T,Rank>(-a)+=f);
}

//multiplication
template<typename T, unsigned Rank, typename U>
marray<T,Rank>& operator*=(marray<T,Rank>& a, const U& f){
	for(auto& i : a) i*=f;
	return(a);
}
template<typename T, unsigned Rank, typename U>
marray<T,Rank>&& operator*=(marray<T,Rank>&& a, const U& f){
	for(auto& i : a) i*=f;
	return(std::move(a));
}

template<typename T, unsigned Rank, typename U,
typename=typename std::enable_if<!detail::is_marray<U>::value>::type>
marray<T,Rank> operator*(marray<T,Rank>& a, const U& f){
	return(marray<T,Rank>(a)*=f);
}
template<typename T, unsigned Rank, typename U,
typename=typename std::enable_if<!detail::is_marray<U>::value>::type>
marray<T,Rank> operator*(const U& f, marray<T,Rank>& a){
	return(marray<T,Rank>(a)*=f);
}

namespace detail{
	//If libdivide is available it can be used to accelerate division of wide integer
	//types, so use a partial specialization to switch that on when appropriate
	
	//general case
	template<typename T, unsigned Rank, typename U,
	bool FastIntDivide=(std::is_same<T,U>::value && (std::is_same<U,int32_t>::value ||
	  std::is_same<U,uint32_t>::value || std::is_same<U,int64_t>::value || std::is_same<U,uint64_t>::value))>
	struct marray_divider{
		static void divide(marray<T,Rank>& a, const U& f){
			for(auto& i : a) i/=f;
		}
	};
#ifdef MARRAY_USE_LIBDIVIDE
	//special case for accelaration
	template<typename T, unsigned Rank, typename U>
	struct marray_divider<T,Rank,U,true>{
		static void divide(marray<T,Rank>& a, const U& f){
			libdivide::divider<U> fast_d(f);
			switch (fast_d.get_algorithm()) {
				case 0: for(auto& i : a) i=i/libdivide::unswitch<0>(fast_d); break;
				case 1: for(auto& i : a) i=i/libdivide::unswitch<1>(fast_d); break;
				case 2: for(auto& i : a) i=i/libdivide::unswitch<2>(fast_d); break;
				case 3: for(auto& i : a) i=i/libdivide::unswitch<3>(fast_d); break;
				case 4: for(auto& i : a) i=i/libdivide::unswitch<4>(fast_d); break;
			}
		}
	};
#endif
}

template<typename T, unsigned Rank, typename U>
marray<T,Rank>& operator/=(marray<T,Rank>& a, const U& f){
	detail::marray_divider<T,Rank,U>::divide(a,f);
	return(a);
}
template<typename T, unsigned Rank, typename U>
marray<T,Rank>&& operator/=(marray<T,Rank>&& a, const U& f){
	detail::marray_divider<T,Rank,U>::divide(a,f);
	return(std::move(a));
}

template<typename T, unsigned Rank, typename U,
typename=typename std::enable_if<!detail::is_marray<U>::value>::type>
marray<T,Rank> operator/(marray<T,Rank>& a, const U& f){
	return(marray<T,Rank>(a)/=f);
}
template<typename T, unsigned Rank, typename U>
marray<T,Rank> operator/(const U& f, marray<T,Rank>& a){
	return(marray<T,Rank>(a)/=f);
}

namespace detail{
	template<typename T, unsigned Rank>
	marray<T,Rank> contruct_min_init(const std::array<size_t,Rank>& extents, std::true_type init_needed){
		return(marray<T,Rank>(extents));
	}
	template<typename T, unsigned Rank>
	marray<T,Rank> contruct_min_init(const std::array<size_t,Rank>& extents, std::false_type init_not_needed){
		return(marray<T,Rank>(extents,no_init));
	}
}

//Generic binary operation between marrays whose types, ranks, or sizes may not match
template<typename T1, unsigned R1, typename T2, unsigned R2, typename F,
		 typename OT=typename std::result_of<F(T1,T2)>::type, //output type
         unsigned OR=detail::Max<R1,R2>::value //output rank
>
marray<OT,OR> applyOp(const marray<T1,R1>& m1, const marray<T2,R2>& m2, const F& op){
	//compute the result extents, or yell at the user if the inputs cannot be reconciled
	std::array<size_t,OR> resultExtents;
	for(unsigned int i1=0,i2=0,j=0; j<OR; j++){
		size_t extent=0;
		if(i1==R1){
			extent=m2.extent(R2-i2-1);
			i2++;
		}
		else if(i2==R2){
			extent=m1.extent(R1-i1-1);
			i1++;
		}
		else{
			size_t e1=m1.extent(R1-i1-1);
			size_t e2=m2.extent(R2-i2-1);
			if(e1==e2)
				extent=e1;
			else if(e1==1)
				extent=e2;
			else if(e2==1)
				extent=e1;
			else{
				std::ostringstream ss;
				ss << "Incompatible array shapes: \n";
				ss << " array 1 dimension " << R1-i1-1 << " has extent " << e1 << '\n';
				ss << " array 2 dimension " << R2-i2-1 << " has extent " << e2 << '\n';
				throw std::logic_error(ss.str());
			}
			i1++;
			i2++;
		}
		resultExtents[OR-j-1]=extent;
	}
	//allocate output
	constexpr bool init_needed=!std::is_trivial<OT>::value;
	marray<OT,OR> result=
	  detail::contruct_min_init<OT,OR>(resultExtents,
                                       std::integral_constant<bool,init_needed>());
	
	std::array<size_t,R1+1> in1; //indices into m1
	std::array<size_t,R2+1> in2; //indices into m2
	//result indices are implicit
	std::fill(in1.begin(),in1.end(),0);
	std::fill(in2.begin(),in2.end(),0);
	auto it1=m1.begin();
	auto it2=m2.begin();
	
	std::array<size_t,R1+1> ex1; //shifted extents of m1
	std::array<size_t,R2+1> ex2; //shifted extents of m1
	for(unsigned int i=0; i<R1; i++)
		ex1[i+1]=m1.extent(i)-1;
	for(unsigned int i=0; i<R2; i++)
		ex2[i+1]=m2.extent(i)-1;
	
	std::array<size_t,R1> st1; //strides for m1
	std::array<size_t,R2> st2; //strides for m2
	{
		size_t stride=1;
		for(size_t i=0; i<R1; i++){
			st1[i]=stride;
			stride*=m1.extent(R1-i-1);
		}
		stride=1;
		for(size_t i=0; i<R2; i++){
			st2[i]=stride;
			stride*=m2.extent(R2-i-1);
		}
	}
	
	//compute the contents of the result
	for(auto it=result.begin(), endIt=result.end(); it!=endIt; it++){
		*it=op(*it1,*it2);
		
		//update the indices
		for(size_t i=0; i<OR; i++){
			size_t d1=R1-i;
			size_t d2=R2-i;
			//if either index is beyond its corresponding rank,
			//we must just need to work on the other
			if(i>=R1){
				if(in2[d2]<ex2[d2]){
					in2[d2]++;
					break;
				}
				in2[d2]=0;
			}
			else if(i>=R2){
				if(in1[d1]<ex1[d1]){
					in1[d1]++;
					break;
				}
				in1[d1]=0;
			}
			else{ //both are in range
				//successfully incrementing either index indicates completion
				bool done=false;
				//if an extent is one we always set the corresponding index to 0
				if(ex1[d1]>0 && in1[d1]<ex1[d1]){
					in1[d1]++;
					done=true;
				}
				else
					in1[d1]=0;
				if(ex2[d2]>0 && in2[d2]<ex2[d2]){
					in2[d2]++;
					done=true;
				}
				else
					in2[d2]=0;
				if(done)
					break;
			}
		}
		//obtain new iterators
		it1=m1.begin();
		for(size_t i=0; i<R1; i++)
			it1+=st1[i]*in1[R1-i];
		it2=m2.begin();
		for(size_t i=0; i<R2; i++)
			it2+=st2[i]*in2[R2-i];
	}
	
	//done
	return(result);
}

namespace detail{
	//Core functions for binary operations on arrays with mathcing types, ranks, and shapes.
	//If the contained type is trivial initialization can be skipped, for a few percent of
	//time saved. Otherwise, copy-construct from the left operand and use an assignment form.
	#define same_shape_binary_op(name,op) \
	template<typename T, unsigned Rank> \
	marray<T,Rank> name ## _matching(const marray<T,Rank>& a1, const marray<T,Rank>& a2, std::true_type init_needed){ \
		return(marray<T,Rank>(a1) op ## = a2); \
	} \
	template<typename T, unsigned Rank> \
	marray<T,Rank> name ## _matching(const marray<T,Rank>& a1, const marray<T,Rank>& a2, std::false_type init_not_needed){ \
		marray<T,Rank>r(a1.get_extents(),detail::no_init); \
		auto i1=a1.cbegin(), i2=a2.cbegin(); \
		for(auto& e : r) e=*i1++ op *i2++; \
		return(r); \
	}
	
	same_shape_binary_op(sum,+)
	same_shape_binary_op(difference,-)
	same_shape_binary_op(product,*)
	same_shape_binary_op(quotient,/)
	
	#undef same_shape_binary_op
}

//Addition
//Fully general case: different types or ranks
template<typename T1, unsigned R1, typename T2, unsigned R2,
         typename OT=decltype(std::declval<T1>()+std::declval<T2>()),
         unsigned OR=detail::Max<R1,R2>::value //output rank
        >
marray<OT,OR> operator+(const marray<T1,R1>& a1, const marray<T2,R2>& a2){
	return(applyOp(a1,a2,[](const T1& t1, const T2& t2){return(t1+t2);}));
}

//special case: same type, same rank
template<typename T, unsigned Rank>
marray<T,Rank> operator+(const marray<T,Rank>& a1, const marray<T,Rank>& a2){
	//if shapes match, use fast path
	if(a1.get_extents()==a2.get_extents()){
		constexpr bool init_needed=!std::is_trivial<T>::value;
		return(detail::sum_matching(a1,a2,std::integral_constant<bool,init_needed>()));
	}
	//otherwise, use slow, general method
	return(applyOp(a1,a2,[](const T& t1, const T& t2){return(t1+t2);}));
}

//Subtraction
//Fully general case: different types or ranks
template<typename T1, unsigned R1, typename T2, unsigned R2,
         typename OT=decltype(std::declval<T1>()-std::declval<T2>()),
         unsigned OR=detail::Max<R1,R2>::value //output rank
        >
marray<OT,OR> operator-(const marray<T1,R1>& a1, const marray<T2,R2>& a2){
	return(applyOp(a1,a2,[](const T1& t1, const T2& t2){return(t1-t2);}));
}

//special case: same type, same rank
template<typename T, unsigned Rank>
marray<T,Rank> operator-(const marray<T,Rank>& a1, const marray<T,Rank>& a2){
	//if shapes match, use fast path
	if(a1.get_extents()==a2.get_extents()){
		constexpr bool init_needed=!std::is_trivial<T>::value;
		return(detail::difference_matching(a1,a2,std::integral_constant<bool,init_needed>()));
	}
	//otherwise, use slow, general method
	return(applyOp(a1,a2,[](const T& t1, const T& t2){return(t1-t2);}));
}

//Multiplication
//Fully general case: different types or ranks
template<typename T1, unsigned R1, typename T2, unsigned R2,
         typename OT=decltype(std::declval<T1>()*std::declval<T2>()),
         unsigned OR=detail::Max<R1,R2>::value //output rank
        >
marray<OT,OR> operator*(const marray<T1,R1>& a1, const marray<T2,R2>& a2){
	return(applyOp(a1,a2,[](const T1& t1, const T2& t2){return(t1*t2);}));
}

//special case: same type, same rank
template<typename T, unsigned Rank>
marray<T,Rank> operator*(const marray<T,Rank>& a1, const marray<T,Rank>& a2){
	//if shapes match, use fast path
	if(a1.get_extents()==a2.get_extents()){
		constexpr bool init_needed=!std::is_trivial<T>::value;
		return(detail::product_matching(a1,a2,std::integral_constant<bool,init_needed>()));
	}
	//otherwise, use slow, general method
	return(applyOp(a1,a2,[](const T& t1, const T& t2){return(t1*t2);}));
}

//Division
//Fully general case: different types or ranks
template<typename T1, unsigned R1, typename T2, unsigned R2,
         typename OT=decltype(std::declval<T1>()/std::declval<T2>()),
         unsigned OR=detail::Max<R1,R2>::value //output rank
        >
marray<OT,OR> operator/(const marray<T1,R1>& a1, const marray<T2,R2>& a2){
	return(applyOp(a1,a2,[](const T1& t1, const T2& t2){return(t1/t2);}));
}

//special case: same type, same rank
template<typename T, unsigned Rank>
marray<T,Rank> operator/(const marray<T,Rank>& a1, const marray<T,Rank>& a2){
	//if shapes match, use fast path
	if(a1.get_extents()==a2.get_extents()){
		constexpr bool init_needed=!std::is_trivial<T>::value;
		return(detail::quotient_matching(a1,a2,std::integral_constant<bool,init_needed>()));
	}
	//otherwise, use slow, general method
	return(applyOp(a1,a2,[](const T& t1, const T& t2){return(t1/t2);}));
}

//fancy operations

using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::cosh;
using std::sinh;
using std::tanh;
using std::acosh;
using std::asinh;
using std::atanh;
using std::exp;
using std::exp2;
using std::log;
using std::log2;
using std::log10;
using std::sqrt;
using std::cbrt;
using std::erf;
using std::erfc;
using std::lgamma;
using std::ceil;
using std::floor;
using std::trunc;
using std::round;
using std::abs;
	
#define UNARY_MARRAY_FUNC(fname) \
template<typename T, unsigned R> \
marray<T,R> fname(marray<T,R>&& a){ \
	for(auto i=a.begin(), e=a.end(); i!=e; i++) \
		*i=fname(*i); \
	return(a); \
} \
template<typename T, unsigned R> \
marray<T,R> fname(const marray<T,R>& a){ \
	return(fname(marray<T,R>(a))); \
}

UNARY_MARRAY_FUNC(cos);
UNARY_MARRAY_FUNC(sin);
UNARY_MARRAY_FUNC(tan);
UNARY_MARRAY_FUNC(acos);
UNARY_MARRAY_FUNC(asin);
UNARY_MARRAY_FUNC(atan);
UNARY_MARRAY_FUNC(cosh);
UNARY_MARRAY_FUNC(sinh);
UNARY_MARRAY_FUNC(tanh);
UNARY_MARRAY_FUNC(acosh);
UNARY_MARRAY_FUNC(asinh);
UNARY_MARRAY_FUNC(atanh);
UNARY_MARRAY_FUNC(exp);
UNARY_MARRAY_FUNC(exp2);
UNARY_MARRAY_FUNC(log);
UNARY_MARRAY_FUNC(log2);
UNARY_MARRAY_FUNC(log10);
UNARY_MARRAY_FUNC(sqrt);
UNARY_MARRAY_FUNC(cbrt);
UNARY_MARRAY_FUNC(erf);
UNARY_MARRAY_FUNC(erfc);
UNARY_MARRAY_FUNC(tgmama);
UNARY_MARRAY_FUNC(lgamma);
UNARY_MARRAY_FUNC(ceil);
UNARY_MARRAY_FUNC(floor);
UNARY_MARRAY_FUNC(trunc);
UNARY_MARRAY_FUNC(round);
UNARY_MARRAY_FUNC(abs);

#undef UNARY_MARRAY_FUNC

using std::atan2;
using std::pow;
using std::hypot;
	
#define BINARY_MARRAY_FUNC(fname) \
template<typename T1, unsigned R1, typename T2, unsigned R2, \
         typename OT=decltype(fname(std::declval<T1>(),std::declval<T2>())), \
         unsigned OR=detail::Max<R1,R2>::value \
        > \
marray<OT,OR> fname(const marray<T1,R1>& a1, const marray<T2,R2>& a2){ \
	return(applyOp(a1,a2,[](const T1& t1, const T2& t2){return(fname(t1,t2));})); \
} \
template<typename T, unsigned Rank> \
marray<T,Rank> fname(const marray<T,Rank>& a1, const marray<T,Rank>& a2){ \
	if(a1.get_extents()==a2.get_extents()){ \
		marray<T,Rank> result(a1); \
		for(auto ri=result.begin(),i1=a1.begin(),i2=a2.begin(),re=result.end; \
			ri!=re; ri++,i1++,i2++) \
			*ri=fname(*i1,*i2); \
		return(result); \
	} \
	return(applyOp(a1,a2,[](const T& t1, const T& t2){return(fname(t1,t2));})); \
}

BINARY_MARRAY_FUNC(atan2);
BINARY_MARRAY_FUNC(pow);
BINARY_MARRAY_FUNC(hypot);

#undef BINARY_MARRAY_FUNC

} // end namespace nusquids

#endif //MARRAY_H
