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
 *      $Id: heuristicSearch.cpp 2017-04-07 13:37:27Z pham $
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

#include "heuristicSearch.hpp"
#include "expand_avx.hpp"
#include "utils.hpp"


using namespace std;

////////////////////////////////////////////////////
////Search the largest patterns/////////////////////
////////////////////////////////////////////////////

//print scores of itemset (heuristic search)
void print_itemset_score(Tidset_vector& p, Transaction& att, float& or_threshold, float& rr_threshold, float & arr_threshold)
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

  if(c > 0){
      int b = att.nb_case - a;
      int d = att.nb_control - c;
      //float chi = chi2(a,b,c,d);// cout<<"(chi2: "<<chi<<")";
      float odd = odd_ratio(a,b,c,d);
      float rr = risk_ratio(a,b,c,d);
      float arr = difference_risk(a,b,c,d);
      float lci = LCI(odd,a,b,c,d);
      float uci = UCI(odd,a,b,c,d);
	  //float p_val = p_value(a,b,c,d);
	  cout<<(100*float(a)/(a+b)) <<" : "<<(100*float(c)/(c+d)) <<" : "<<odd<<" : "<<rr<<" : "<<arr<<" : "<<lci<<"-"<<uci;
	  //update thresholds after successfully print a pattern
      ///////////////////////////////////////////////////////////////////////////////////
      if( (odd - or_threshold) > 0.1) or_threshold += 0.1; else or_threshold = odd;/////
      ////////////////////////////////////////////////////////////////////////////////
    }else
		cout<<(100*float(a)/att.nb_case) <<" : "<<(100*float(c)/att.nb_control);
}


void expand_control_heu(Tidset_vector p, Tidlist tid_p, float& or_threshold, float& rr_threshold, float& arr_threshold, int min_case_out, int& nb_patterns, Transaction& att, int nb_registers, int& nb_it)
{
  Tidlist tid = compute_tidlist_avx(p, att);
  Tidset_vector p_ext_control = compute_closure_avx(tid, att, nb_registers, 1);
  int n = nb_registers*256;
  Tidset_vector p_tmp = p;
  for(int i=(n-att.nb_control-1); i>= (n-att.nb_sample); i--) {
      int v = i/256;
      int pos = i%256;
      SetBit(p_tmp[v], pos, 0);
    }
  p_ext_control = remove_tidset_avx(p_ext_control,p_tmp);

  if(!check_empty_avx(p_ext_control)) {
      Tidset_vector q = add_tidset_avx(p, p_ext_control); //q = p U {e} U p_ext
      Tidlist tid_q = compute_tidlist_avx(q, att);
      Tidset_vector p_ext_all = compute_closure_avx(tid_q, att, nb_registers, 2);
      //cout<<"p_ext_all:"; print_itemset(p_ext_all,att); cout<<endl;
      p_ext_all = remove_tidset_avx(q, p_ext_all);

      if( (check_empty_avx(p_ext_all)) && check_itemset_score(q, att, or_threshold, rr_threshold, arr_threshold, min_case_out) ){
	  //cout<<"in"<<endl;
	  nb_patterns++;
	  for(int i=0;i<tid_q.size()-1;i++) cout<<att.tidset[tid_q[i]].label<<",";  cout<<att.tidset[tid_q[tid_q.size()-1]].label;
	  cout<<"(";
	  //print_itemset(q,att);
	  //cout<<": ";
	  print_itemset_score(q, att, or_threshold, rr_threshold, arr_threshold);
	  cout<<")"<<endl;
	  nb_it=0;
	}
   }  else {
      //Itemset_vector p_ext_all = compute_closure_avx(tid, att, nb_registers, 2);
      //p_ext_all = remove_tidset_avx(p, p_ext_all);
      //if( (check_empty_avx(p_ext_all)) && check_itemset_score(p, att, or_threshold, rr_threshold, arr_threshold, min_case_out) )
      if( check_itemset_score(p, att, or_threshold, rr_threshold, arr_threshold, min_case_out) ){
	  nb_patterns++;
	  for(int i=0;i<tid.size()-1;i++) cout<<att.tidset[tid[i]].label<<","; cout<<att.tidset[tid[tid.size()-1]].label;
	  cout<<"(";
	  //print_itemset(p,att);
	  //cout<<": ";
	  print_itemset_score(p, att, or_threshold, rr_threshold, arr_threshold);
	  //cout<<": "<<nb_it;
	  cout<<")"<<endl;
	  nb_it=0;
	}
    }
}

/////////////////////
void expand_case_heu(Tidset_vector p, int e, float& or_threshold, float& rr_threshold, float& arr_threshold, int min_case_out, int it_threshold,  int& nb_patterns, int& nb_prunes, Transaction& att, int nb_registers, int& nb_it)
{
  //count number iteration. If it is equal to the iteration threshold then stop searching
  nb_it++;
  if(nb_it==it_threshold) {
	   cout<<"#nb_patterns: "<<nb_patterns<<endl;
	   exit(1);
   }

  //cout<<endl<<"expand case:"<<e<<endl;
  int n = nb_registers*256;
  int v = (n-att.nb_sample+e) / 256;
  int pos = (n-att.nb_sample+e) % 256;
  SetBit(p[v],pos,true); //p=p U {e}
  //cout<<"new items: ";print_itemset(p,att);cout<<endl;
  Tidlist tid = compute_tidlist_avx(p,att);
  //cout<<"tidlist avx:"; for(int i=0;i<tid.size();i++) cout<<att.tidset[tid[i]].label<<" ";  cout<<endl;
  if(tid.size()>1){
      if(predict_expand_avx(tid, or_threshold, att, nb_registers)){
	  Tidset_vector p_ext_case = compute_closure_avx(tid, att, nb_registers, 0);
	  p_ext_case = remove_tidset_avx(p, p_ext_case);
	  if(!check_empty_avx(p_ext_case)){
	      std::vector<int> max_item = get_bitset_pos(p_ext_case,att);
	      if(max_item[max_item.size()-1] < e){
		 Tidset_vector q = add_tidset_avx(p, p_ext_case); //Q = p U {e} U p_ext
		 Tidlist tid_q = compute_tidlist_avx(q, att);
		 Transaction ratt = reduced_dataset_avx(tid_q, att); //reduced data set

		 Tidset_vector k = remove_tidset_avx(q, att.case_itemset); //K = I+ \Q
		 std::vector<int> k_ext = get_bitset_pos(k,att);

		 //if(k_ext.size()>0)
		 //for(int i=0; i<k_ext.size(); i++) //expanding from small to larger id

		 // #pragma omp parallel for// num_threads(4)

		 for(int i=k_ext.size()-1; i>=0; i--) //reverse version
		     if(k_ext[i]<e)
		       expand_case_heu(q, k_ext[i], or_threshold, rr_threshold, arr_threshold, min_case_out,it_threshold, nb_patterns, nb_prunes, ratt, nb_registers, nb_it);

		 //find discriminative pattern
		 Tidset_vector p_ext_case2 = compute_closure_avx(tid_q, att, nb_registers, 0);
		 p_ext_case2 = remove_tidset_avx(q, p_ext_case2);
		 if(check_empty_avx(p_ext_case2))
		  if(get_size(q) >= min_case_out)
		    expand_control_heu(q, tid_q, or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, ratt, nb_registers, nb_it);
		  ///////////////////////////////////////////////////////////////
		}
	    }
	  else
	    {
	      //expand p with all row ids smaller than min row_id in case
	      Transaction ratt = reduced_dataset_avx(tid, att); //reduced data set
	      
		  //expanding from smaller id to larger id
	      //int min = min_tidset_avx(p);
	      //if(min > 0)
	      //for(int i=0; i<min; i++)
		  //expand_case_heu(p, i, or_threshold, rr_threshold, nb_patterns, nb_prunes, ratt, nb_registers, nb_it);

	      //reverse version: expanding from large id to small id
	      Tidset_vector k = remove_tidset_avx(p, att.case_itemset); //K = I+ \p
	      std::vector<int> k_ext = get_bitset_pos(k,att);

	      for(int i=k_ext.size()-1; i>=0; i--)
	         if(k_ext[i]<e)
	           expand_case_heu(p, k_ext[i], or_threshold, rr_threshold, arr_threshold, min_case_out, it_threshold, nb_patterns, nb_prunes, ratt, nb_registers, nb_it);

	      //find discriminative pattern
	      if(get_size(p) >= min_case_out)
	    	  expand_control_heu(p, tid, or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, ratt, nb_registers, nb_it);
	    }
	}
	else { nb_prunes++; }
    }
}



