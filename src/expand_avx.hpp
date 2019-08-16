/*
 *  Authors:
 *      Hoang-Son Pham <hoang-son.pham@irisa.fr>
 *      Axlexander Termier <alexander.termier@irisa.fr>
 *      Dominique Lavenier <dominique.lavenier@irisa.fr>
 *
 *  Copyright:
 *      Hoang-Son Pham, 2017
 *
 *  Revision information:
 *      $Id: expand_avs.hpp 2017-04-07 13:37:27Z pham $
 *
 *  This file is part of SSDPS,
 *  Statistically Significant Discriminative Pattern Search
 *
 *  Permission is hereby granted, free of charge, to any person obtaining
 *  a copy of this software and associated documentation files (the
 *  "Software"), to deal in the Software without restriction, including
 *  without limitation the rights to use, copy, modify, merge, publish,
 *  distribute, sublicense, and/or sell copies of the Software, and to
 *  permit persons to whom the Software is furnished to do so, subject to
 *  the following conditions:
 *
 *  The above copyright notice and this permission notice shall be
 *  included in all copies or substantial portions of the Software.
 *
 *  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 *  LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 *  OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 *  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 */

#ifndef EXPAND_AVX_HPP_
#define EXPAND_AVX_HPP_

#include <immintrin.h>
#include <iostream>
#include <vector>
#include <bitset>
#include <stdlib.h>

#include "utils.hpp"

using namespace std;


/**
 * Allocator for aligned data.
 * Modified from the Mallocator from Stephan T. Lavavej.
 * <http://blogs.msdn.com/b/vcblog/archive/2008/08/28/the-mallocator.aspx>
 */
template <typename T, std::size_t Alignment>
class aligned_allocator
{
	public:
		// The following will be the same for virtually all allocators.
		typedef T * pointer;
		typedef const T * const_pointer;
		typedef T& reference;
		typedef const T& const_reference;
		typedef T value_type;
		typedef std::size_t size_type;
		typedef ptrdiff_t difference_type;
 
		T * address(T& r) const
		{
			return &r;
		}
 
		const T * address(const T& s) const
		{
			return &s;
		}
 
		std::size_t max_size() const
		{
			// The following has been carefully written to be independent of
			// the definition of size_t and to avoid signed/unsigned warnings.
			return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
		}
 
 
		// The following must be the same for all allocators.
		template <typename U>
		struct rebind
		{
			typedef aligned_allocator<U, Alignment> other;
		} ;
 
		bool operator!=(const aligned_allocator& other) const
		{
			return !(*this == other);
		}
 
		void construct(T * const p, const T& t) const
		{
			void * const pv = static_cast<void *>(p);
 
			new (pv) T(t);
		}
 
		void destroy(T * const p) const
		{
			p->~T();
		}
 
		// Returns true if and only if storage allocated from *this
		// can be deallocated from other, and vice versa.
		// Always returns true for stateless allocators.
		bool operator==(const aligned_allocator& other) const
		{
			return true;
		}
 
 
		// Default constructor, copy constructor, rebinding constructor, and destructor.
		// Empty for stateless allocators.
		aligned_allocator() { }
 
		aligned_allocator(const aligned_allocator&) { }
 
		template <typename U> aligned_allocator(const aligned_allocator<U, Alignment>&) { }
 
		~aligned_allocator() { }
 
 
		// The following will be different for each allocator.
		T * allocate(const std::size_t n) const
		{
			// The return value of allocate(0) is unspecified.
			// Mallocator returns NULL in order to avoid depending
			// on malloc(0)'s implementation-defined behavior
			// (the implementation can define malloc(0) to return NULL,
			// in which case the bad_alloc check below would fire).
			// All allocators can return NULL in this case.
			if (n == 0) {
				return NULL;
			}
 
			// All allocators should contain an integer overflow check.
			// The Standardization Committee recommends that std::length_error
			// be thrown in the case of integer overflow.
			if (n > max_size())
			{
			  //throw std::cout("aligned_allocator<T>::allocate() - Integer overflow.");
			}
 
			// Mallocator wraps malloc().
			void * const pv = _mm_malloc(n * sizeof(T), Alignment);
 
			// Allocators should throw std::bad_alloc in the case of memory allocation failure.
			if (pv == NULL)
			{
				throw std::bad_alloc();
			}
 
			return static_cast<T *>(pv);
		}
 
		void deallocate(T * const p, const std::size_t n) const
		{
			_mm_free(p);
		}
 
 
		// The following will be the same for all allocators that ignore hints.
		template <typename U>
		T * allocate(const std::size_t n, const U * /* const hint */) const
		{
			return allocate(n);
		}
 
 
		// Allocators are not required to be assignable, so
		// all allocators should have a private unimplemented
		// assignment operator. Note that this will trigger the
		// off-by-default (enabled under /Wall) warning C4626
		// "assignment operator could not be generated because a
		// base class assignment operator is inaccessible" within
		// the STL headers, but that warning is useless.
	private:
		aligned_allocator& operator=(const aligned_allocator&);
};


///////////////////////////////////////////////////////////////////////////////////////////////////
/*
Note: the algorithm uses a transposition matrix as input. Therefore, the notation:
Tidlist = column of the matrix
Itemset = row of the matrix
 */

typedef std::vector<__m256i, aligned_allocator<__m256i, sizeof(__m256i)> > aligned_vector;

const int nb_bits = 256;//nb of bits per register (AVX2=256bits)
const int nb_chunks = 8; //nb of integer values per register = 256/32

typedef std::vector<int> Tidlist; //store ids of items in transposition matrix

struct ITEM{
  int id;
  int label;
};
typedef std::vector<ITEM> Itemset;
//set of items:
//each item is a pair of id-label
//it stores the original ids of an items when using reduced dataset

typedef __m256i Tidset; //tidset of data which have transactions (number of individuals) that are less than 256 bits

typedef aligned_vector Tidset_vector; //tidset of data which have transactions (number of individuals) that are larger than 256 bits
//transposition transaction dataset is a set of vector
//each vector is an set of tids
//the size of transtraction table = the number of items
struct Transaction : std::vector<Tidset_vector> //define data struture
{
  int nb_sample;
  int nb_case;
  int nb_control;
  Tidset_vector case_itemset;
  Tidset_vector control_itemset;
  Itemset tidset;
};

/////////////////////////////////////////////////////////////////////////////////////////

int popcount(int v);

int get_size(Tidset_vector& a);

void SetBit(Tidset & vector, size_t position, bool value);

int check_empty_avx (Tidset_vector& a);

int min_tidset_avx(Tidset_vector& a, Transaction& att);

int max_tidset_avx(Tidset_vector& a, Transaction& att);

Tidset_vector add_tidset_avx(Tidset_vector& a, Tidset_vector& b);

Tidset_vector remove_tidset_avx(Tidset_vector& a, Tidset_vector& b);

int check_itemset_score(Tidset_vector& p, Transaction& att, float or_threshold, float rr_threshold, float arr_threshold, int minCase);

std::vector<int> get_bitset_pos(Tidset_vector& a,  Transaction& att);

void print_itemset(Tidset_vector& a, Transaction& att);

Tidlist compute_tidlist_avx(Tidset_vector& p, Transaction& att);

Tidset_vector compute_closure_avx(Tidlist tid, Transaction& att, int nb_registers, int option);

Transaction reduced_dataset_avx(Tidlist tid, Transaction& att);

int predict_expand_avx(Tidlist tid, float threshold, Transaction& att, int nb_registers );


#endif /* EXPAND_AVX_HPP_ */
