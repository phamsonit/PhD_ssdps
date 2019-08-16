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
 *      $Id: expand_avx.cpp 2017-04-07 13:37:27Z pham $
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
#include <immintrin.h>
#include <iostream>
#include <vector>
#include <bitset>

#include "expand_avx.hpp"
#include "utils.hpp"


using namespace std;

/////////////////////////////////////////////////
/////////////////////AVX2///////////////////////
////////////////////////////////////////////////

//count number of bits in an integer
int popcount(int v)
 {
    v = v - ((v >> 1) & 0x55555555);                // put count of each 2 bits into those 2 bits
    v = (v & 0x33333333) + ((v >> 2) & 0x33333333); // put count of each 4 bits into those 4 bits  
    return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24;
 }

//return the size of an Tidset = number of bits 1 in an tidset
int get_size(Tidset_vector& a)
{
  int result = 0;
  for(int i=0; i<a.size();i++){
      int* f = (int*)&a[i];
      for(int j=0; j<nb_chunks; j++)
    	  result += popcount(f[j]);
    }
  return result;
}

//set a value (0 or 1) at a position of tidset
void SetBit(Tidset & vector, size_t position, bool value)
{
  //    assert(position <= 255);

    uint8_t lut[32] = { 0 };
    lut[position >> 3] = 1 << (position & 7);
    __m256i mask = _mm256_loadu_si256((__m256i*)lut);
    if (value)
        vector = _mm256_or_si256(mask, vector);
    else
        vector = _mm256_andnot_si256(mask, vector);
}

//check tidset empty? return 1: empty; 0: not empty
int check_empty_avx (Tidset_vector& a)
{
  int check = 1;
  for(int i=0; i<a.size(); i++){
      int* f = (int*)&a[i];
      for(int l=0; l<nb_chunks; l++)
    	  if (f[l]!=0){
    		  check = 0;
    		  break;
    	  }
    	}
  return check;
}


//return minimal tid of an tidset  = find the first position of tidset that is set by 1
int min_tidset_avx(Tidset_vector& a, Transaction& att)
{
  int pos = 0;
  for(int i=0; i<a.size(); i++){
      int* f = (int*)&a[i];
      for(int l=0; l<nb_chunks; l++){
    	  for(int x=0; x<32; x++) {
    		  int bit = (f[l] >> x) & 1;
    		  if(bit) return att.nb_sample - (a.size()*256-pos);
    		  else pos++;
    	  }
      }
    }
	return pos;
}

//return maximal tid of an tidset = find the last position of tidset that is set by 1
int max_tidset_avx(Tidset_vector& a,Transaction& att)
{
  int pos = att.nb_sample-1;
  for(int i=a.size()-1; i>=0; i++){
      int* f = (int*)&a[i];
      for(int l=nb_chunks-1; l>=0; l--)	{
    	  for(int x=31; x>=0; x--) {
    		  int bit = (f[l] >> x) & 1;
    		  if(bit) return pos;
    		  else pos--;
    	  }
      }
    }
	return pos;
}

//compute intersection of two tidset
Tidset_vector add_tidset_avx(Tidset_vector& a, Tidset_vector& b)
{
  Tidset_vector tmp;
  for(int i=0; i<a.size();i++)
    tmp.push_back(_mm256_or_si256(a[i],b[i]));
  return tmp;
}

//compute subsection of two tidsets
Tidset_vector remove_tidset_avx(Tidset_vector& a, Tidset_vector& b)
{
  Tidset_vector tmp;
  for(int i=0; i<a.size();i++)
    tmp.push_back(_mm256_xor_si256(a[i],b[i]));
  return tmp;
}


//find the positions which are set by 1 in an tidset = list of transaction ids
std::vector<int> get_bitset_pos(Tidset_vector& a, Transaction& att)
{
  std::vector<int> result;
  int pos = a.size()*256;
  for(int i=0; i<a.size(); i++)
    {
      int* f = (int*)&a[i];
      for(int l=0; l<nb_chunks; l++)
	{
	  for(int x=0; x<32; x++) {
	    int bit = (f[l] >> x) & 1;
	    if(bit) {result.push_back(att.nb_sample-pos); pos--;}
	    else pos--;
	  }
	}
    }
  return result;
}

////////////////////////////////
void print_itemset(Tidset_vector& a, Transaction& att)
{
  int n = a.size()*256-att.nb_sample;
  int count = 0;
  // cout<<"(";
  for(int i=0; i<a.size(); i++)
    {
      int* f = (int*)&a[i];
      for(int l=0; l<nb_chunks; l++)
	{
	  for(int x=0; x<32; x++) 
	    {
	      int bit = (f[l] >> x) & 1;
	      if(bit) cout<<count-n<<" ";
	      count++;
	    }
	}
    }
  //  cout<<")";
}

//compute intersection of two tidsets
Tidlist compute_tidlist_avx(Tidset_vector& p, Transaction& att)
{
  Tidlist tid;
  // #pragma omp parallel for num_threads(nb_threads)
  for(int i=0; i<att.tidset.size(); i++)
    {
      //if( _mm256_testc_si256(_mm256_loadu_si256(&att[i][0]) , p ) ) tid.push_back(i);
      int founded[p.size()];
      for(int j=0; j<p.size(); j++)
	//founded[j] = _mm256_testc_si256(_mm256_loadu_si256(&att[i][j]) , p ); //px[j] is subset of att[][] ???
	founded[j] = _mm256_testc_si256(att[i][j], p[j] ); //px[j] is subset of att[][] ???
      bool found = true;
      for(int k=0; k<p.size(); k++) if(founded[k]==0) { found=false; break;}
      if(found) tid.push_back(i);
     }   
  return tid;
}

//compute closure of an tidset
Tidset_vector compute_closure_avx(Tidlist tid, Transaction& att, int nb_registers, int option)
//option: 0 - case, 1 - control, 2 - all data
{
  //compute intersection
  Tidset_vector result;
  for(int i=0; i<nb_registers; i++) result.push_back(att[tid[0]][i]);

  for(int i=1; i<tid.size(); i++)
    for(int j=0; j<nb_registers; j++)
      {
	result[j] = _mm256_and_si256(result[j] , att[tid[i]][j]) ;
      }
 
  int n = nb_registers*256;
  switch (option)
    {
    case 0:
      for(int i=(n-1); i>= (n-att.nb_control); i--)
	{
	  int v = i/256;
	  int pos = i%256;
	  SetBit(result[v], pos, 0);
	}
      return result;
      break;

    case 1:
      for(int i=(n-att.nb_control-1); i>= (n-att.nb_sample); i--)
	{
	  int v = i/256;
	  int pos = i%256;
	  SetBit(result[v], pos, 0);
	}
      return result;
      break;

    case 2:
      return result;
      break;
    }

  return result;
}

//reduce dataset att with regard to a given tidset
Transaction reduced_dataset_avx(Tidlist tid, Transaction& att)
{
  Transaction dtt;
  for(int i=0; i< tid.size();i++)
    {
      dtt.push_back(att[tid[i]]);
      ITEM tid_tmp;
      tid_tmp.id = i;
      tid_tmp.label = att.tidset[tid[i]].label;
      dtt.tidset.push_back(tid_tmp);
    }
  dtt.nb_sample = att.nb_sample;
  dtt.nb_case = att.nb_case;
  dtt.nb_control = att.nb_control;
  dtt.case_itemset = att.case_itemset;
  dtt.control_itemset = att.control_itemset;

  return dtt;
}

//heuristic: predict expanding based on odds ratio
int predict_expand_avx(Tidlist tid, float threshold, Transaction& att, int nb_registers )
{
  int nb_case_ext = 0;
  int nb_control_ext = 0;
  Tidset_vector p_ext_all = compute_closure_avx(tid,att,nb_registers,2);

  int n =  nb_registers*256;
  int count = nb_registers*256-1;
  for(int i=nb_registers-1; i>=0; i--)
    {
      int* f = (int*)&p_ext_all[i];
      for(int l=nb_chunks-1; l>=0; l--)
	{
	  for(int k=31; k>=0; k--) 
	    {
	      if( f[l]&(1<<k) ){
			  if(count >= n-att.nb_control ){ nb_control_ext++;}
				else if(count >= n-att.nb_sample ) {nb_case_ext++;}
		}
	      count--;	      
	    }
	}
    }
  //cout<<nb_case_ext<<" "<<nb_control_ext<<endl;
  //  cout<<"# case "<<nb_case_ext<<"; # control "<<nb_control_ext<<endl;
  int p_ext_control_size = (threshold*nb_control_ext*att.nb_case) / (nb_control_ext*(threshold-1) + att.nb_control);
  if(nb_case_ext >= p_ext_control_size) return true;  else return false;
}

//check discriminative scores of an itemset ~ pruning
int check_itemset_score(Tidset_vector& p, Transaction& att, float or_threshold, float rr_threshold, float arr_threshold, int min_case_out)
{
  int a = 0;
  int c = 0;
  int count = p.size()*256;
  for(int i=0; i<p.size(); i++){
      int* f = (int*)&p[i];
      for(int l=0; l<nb_chunks; l++){
    	  for(int k=0;k<32;k++){
    		  if( f[l]&(1<<k) ){ //check if the value at position k of f[l] is set 1
    			  if(count <= att.nb_control ){ c++;}
    			  else if(count <= att.nb_sample ) {a++;}
    		  }
	      count--;
    	  }
      }
  }

  //cout<<a<<","<<c<<endl;
  //c > 0:  exist tids in the control group
  if(c > 0){
      int b = att.nb_case - a;
      int d = att.nb_control - c;
      float odd = odd_ratio(a,b,c,d);
      float rr = risk_ratio(a,b,c,d);
      float arr = difference_risk(a,b,c,d);
      float lci_or = LCI(odd,a,b,c,d);
      float lci_rr = R_LCI(rr,a,b,c,d);
      //float uci = UCI(odd,a,b,c,d);
      ////////////////////////////////////////////////////////////////////////////
      if(
    		  (odd >= or_threshold)  &&
			  (rr >= rr_threshold)   &&
			  (arr >= arr_threshold) &&
			  (lci_or >= or_threshold) &&
			  (lci_rr >= rr_threshold)
	  	  )return true;
      else
    	  return false;///
    }else
    	if((a-c) > min_case_out)
    		return true;
    	else
    		return false;
}
