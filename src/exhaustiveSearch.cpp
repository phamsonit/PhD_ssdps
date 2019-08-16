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
 *      $Id: exhaustiveSearch.cpp 2017-04-07 13:37:27Z pham $
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

#include "exhaustiveSearch.hpp"
#include "expand_avx.hpp"
#include "utils.hpp"


using namespace std;


//////////////////////////////////////////////////////
//////////exhaustive search//////////////////////////
/////////////////////////////////////////////////////

//print scores of itemset in exhaustive search
void print_itemset_score_exh(Tidset_vector& p, Transaction& att, float& or_threshold, float& rr_threshold, float & arr_threshold)
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

  if(c>0) {
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
    } else 
		cout<<(100*float(a)/att.nb_case) <<" : "<<(100*float(c)/att.nb_control);
}


//expand pattern to tids in control group
void expand_control_exh(Tidset_vector p, int e, float& or_threshold, float& rr_threshold, float& arr_threshold, int min_case_out, int& nb_patterns, Transaction& att, int nb_registers, int& nb_pruning_control)
{
  //cout<<endl<<"expand control: "<<e<<endl;
  int n = nb_registers*256;
  int v = (n-att.nb_sample+e) / 256;
  int pos = (n-att.nb_sample+e) % 256;
  SetBit(p[v],pos,true); //p=p U {e}

  //cout<<"itemset :";print_itemset(p,att);cout<<endl;
  Tidlist tid = compute_tidlist_avx(p, att);
  // cout<<"tidlist avx:"; for(int i=0;i<tid.size();i++) cout<<att.tidset[tid[i]].label<<" ";  cout<<endl;
  if(tid.size()>1) {
    if(check_itemset_score(p, att, or_threshold, rr_threshold, arr_threshold, min_case_out)){
	  nb_pruning_control++; //count the number of positive procedure calling
	  Tidset_vector p_ext_control = compute_closure_avx(tid, att, nb_registers, 1);
	  Tidset_vector p_tmp = p;
	  for(int i=(n-att.nb_control-1); i>= (n-att.nb_sample); i--) {
	      int v = i/256;
	      int pos = i%256;
	      SetBit(p_tmp[v], pos, 0);
	    }
	  p_ext_control = remove_tidset_avx(p_ext_control,p_tmp);
	  //cout<<"p_ext_control: ";print_itemset(p_ext_control,att);cout<<endl;

	  if(!check_empty_avx(p_ext_control)) {
	      std::vector<int> max_item = get_bitset_pos(p_ext_control,att);
	      //if(max_itemset_avx(p_ext_control) < e)
	      if(max_item[max_item.size()-1] < e){
		  Tidset_vector q = add_tidset_avx(p, p_ext_control); //q = p U {e} U p_ext
		  Tidlist tid_q = compute_tidlist_avx(q, att);
		  Tidset_vector p_ext_all = compute_closure_avx(tid_q, att, nb_registers, 2);
		  p_ext_all = remove_tidset_avx(q, p_ext_all);
		  //cout<<"p_ext_all inside";print_itemset(p_ext_all,att);
		  if( (check_empty_avx(p_ext_all)) && check_itemset_score(q, att, or_threshold, rr_threshold, arr_threshold, min_case_out) ) {
		    if(get_size(q) >= min_case_out)  {
			//print discriminative pattern
			//cout<<"in  : "<<e<<" : ";
		       nb_patterns++;
		       //cout<<tid.size()<<" : ";
		       for(int i=0; i<tid_q.size(); i++) cout<<att.tidset[tid_q[i]].label<<" "; 
		       //cout<<"itemset:";
		       //cout<<": ";
		       //print_itemset(q,att);
		       //cout<<"odds ratio:";
		       cout<<"(";
		       print_itemset_score_exh(q, att, or_threshold, rr_threshold, arr_threshold);
		       //cout<<": "<<nb_it;
		       cout<<")"<<endl;
		      }
		       
		       /////////////////////////find all discriminative patterns///////////////////////
		       Tidset_vector t = add_tidset_avx(p_tmp,p_ext_control);
		       Tidset_vector k = remove_tidset_avx(t, att.control_itemset); //k = I- \
		       //cout<<"k: "; print_itemset(k);cout<<endl;
		       Transaction ratt = reduced_dataset_avx(tid_q, att); //reduced data set
		       std::vector<int> k_ext = get_bitset_pos(k,att);
		       //cout<<"k_ext: "; for(int i=0;i<k_ext.size();i++) cout<<k_ext[i]<<" "; cout<<endl;
		       for(int i=0;i<k_ext.size();i++)
			 if( (k_ext[i] >= att.nb_control) && (k_ext[i]<e) )
			   expand_control_exh(q, k_ext[i], or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, ratt, nb_registers, nb_pruning_control);
		       /////////////////////////////////////////////////////////////////////////////////////
		     }
		}
	    }else {
	      Tidset_vector p_ext_all = compute_closure_avx(tid, att, nb_registers, 2);
	      p_ext_all = remove_tidset_avx(p, p_ext_all);
	      //cout<<"p_ext_all: ";print_itemset(p_ext_all,att);cout<<endl;
	      //set all control = 0
	      if(check_empty_avx(p_ext_all))
		if(check_itemset_score(p, att, or_threshold, rr_threshold, arr_threshold, min_case_out)){
		  if(get_size(p) >= min_case_out){
		      //cout<<"out : "<<e<<" : ";
		      nb_patterns++;
		      //cout<<"tidlist:"; 
		      for(int i=0; i<tid.size(); i++) cout<<att.tidset[tid[i]].label<<" ";
		      //cout<<"itemset:";
		      //cout<<": ";
		      //print_itemset(p,att);
		      //cout<<"odds ratio:";
		      cout<<"(";
		      print_itemset_score_exh(p, att, or_threshold, rr_threshold, arr_threshold);
		      //cout<<": "<<nb_it;
		      cout<<")"<<endl;
		    }
		
		  ////find all discriminative patterns////////////////////////////////////////
		  Tidset_vector k = remove_tidset_avx(p_tmp, att.control_itemset); //k = I- \p
		  //cout<<"k: "; print_itemset(k);cout<<endl;
		  std::vector<int> k_ext = get_bitset_pos(k,att);
		  //cout<<"k_ext: ";  for(int i=0;i<k_ext.size();i++) cout<<k_ext[i]<<" "; cout<<endl;
		  Transaction ratt = reduced_dataset_avx(tid, att); //reduced data set
		  for(int i=0; i<k_ext.size(); i++)
		    if( (k_ext[i] >= att.nb_control) && (k_ext[i]<e) )
		      expand_control_exh(p, k_ext[i], or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, ratt, nb_registers,nb_pruning_control);
		  ////////////////////////////////////////////////////////////////////
		}
	    }
	}
    }//  if(tidlist.size()>1)
}

/////////////////////////////////
void expand_case_exh(Tidset_vector p, int e, float& or_threshold, float& rr_threshold, float& arr_threshold, int min_case_out, int& nb_patterns, int& nb_pruning_case, Transaction& att, int nb_registers, int& nb_pruning_control)
{
  //cout<<endl<<"expand case:"<<e<<endl;
  int n = nb_registers*256;
  int v = (n-att.nb_sample+e) / 256;
  int pos = (n-att.nb_sample+e) % 256;
  SetBit(p[v],pos,true); //p=p U {e}
  //cout<<"new items: ";print_itemset(p,att);cout<<endl;
  Tidlist tid = compute_tidlist_avx(p,att);
  //cout<<"tidlist avx:"; for(int i=0;i<tid.size();i++) cout<<att.tidset[tid[i]].label<<" ";  cout<<endl;
  if(tid.size()>1) {
     //if(predict_expand_avx(tid, or_threshold, att, nb_registers)){
	  Tidset_vector p_ext_case = compute_closure_avx(tid, att, nb_registers, 0); 
	  p_ext_case = remove_tidset_avx(p, p_ext_case);
	  //cout<<"p_ext_case: ";print_itemset(p_ext_case,att);cout<<endl;
	  //cout<<"max p_ext = "<<max_itemset_avx(p_ext_case)<<endl;
	  if(!check_empty_avx(p_ext_case)){
	      std::vector<int> max_item = get_bitset_pos(p_ext_case,att);
	      if(max_item[max_item.size()-1] < e){
		 Tidset_vector q = add_tidset_avx(p, p_ext_case); //Q = p U {e} U p_ext
		 Tidset_vector k = remove_tidset_avx(q, att.case_itemset); //K = I+ \Q
		 Tidlist tid_q = compute_tidlist_avx(q, att);
		 Transaction ratt = reduced_dataset_avx(tid_q, att); //reduced data set

		 std::vector<int> k_ext = get_bitset_pos(k,att);
		 //if(k_ext.size()>0)
		 for(int i=0; i<k_ext.size(); i++)
		     if(k_ext[i]<e)
		       expand_case_exh(q, k_ext[i], or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, nb_pruning_case, ratt, nb_registers,nb_pruning_control);


		 //check if closure of q is empty in control then output q		 
		 Tidset_vector p_ext_control = compute_closure_avx(tid_q, att, nb_registers, 1);
		 if(check_empty_avx(p_ext_control) && (get_size(q) >= min_case_out)) {
		     for(int i=0;i<tid_q.size();i++) cout<<att.tidset[tid_q[i]].label<<" ";
		     cout<<"(";
		     print_itemset_score_exh(q, att, or_threshold, rr_threshold, arr_threshold);
		     cout<<")"<<endl;
		     nb_patterns++;
		   }
		 
		 //expand q with all row ids in control
		   for(int i=att.nb_case; i<att.nb_sample; i++)
		     expand_control_exh(q, i, or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, ratt, nb_registers,nb_pruning_control);
		  ///////////////////////////////////////////////////////////////
		}
	    } else {
	      //expand p with all row ids smaller than min row_id in case
	      Transaction ratt = reduced_dataset_avx(tid, att); //reduced data set
	      //find from small id to larger id
	      std::vector<int> min_item = get_bitset_pos(p,att);
	      int min = min_item[0];
	      //int min = min_tidset_avx(p);
	      if(min > 0)
		  for(int i=0; i<min; i++)
		  expand_case_exh(p, i, or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, nb_pruning_case, ratt, nb_registers,nb_pruning_control);

	      //check if closure of p is empty in control then print p    
	      Tidset_vector p_ext_control = compute_closure_avx(tid, att, nb_registers, 1);
	      if(check_empty_avx(p_ext_control) && (get_size(p) >= min_case_out)){
		  for(int i=0;i<tid.size();i++) cout<<att.tidset[tid[i]].label<<" ";
		  cout<<"(";
		  print_itemset_score_exh(p, att, or_threshold, rr_threshold, arr_threshold);
		  cout<<")"<<endl;
		  nb_patterns++;
		}  
	      //expand q with all row ids in control
	      // if(get_size(p) >= min_case_out)
		for(int i=att.nb_case; i<att.nb_sample; i++)
		  expand_control_exh(p, i, or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, ratt, nb_registers,nb_pruning_control);
	    }
		//}
	//else { nb_pruning_case++; }
    }
}
