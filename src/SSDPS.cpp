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
 *      $Id: SSDPS.cpp 2017-04-07 13:37:27Z pham $
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

//============================================================================
//============================================================================
#include <immintrin.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <bitset> 
#include <map>
#include <stdlib.h>

#include "utils.hpp"
#include "exhaustiveSearch.hpp"
#include "heuristicSearch.hpp"

using namespace std;

  //////////////////////////
  /////global variants//////
  //////////////////////////
  string input_file;   //binary matrix file name
  int nb_case = 0;           //number of cases (positive individuals)
  int nb_control = 0;        //number of controls (negative individuals)
  float or_threshold = 1;    //odds ratio threshold
  float rr_threshold = 1;    //risk ratio threshold
  float arr_threshold = 0;   //absolute risk reduction threshold
  float lci_threshold = 1;   //lower confidence interval of odd ratio
  float p_val = 0;           //p-value threshold 
  float max_control = 0;     //maximal number of control individuals containing items
  float min_case = 0;        //minimal number of case individuals containing items
  float min_case_out = 0;    //minimal number of case individual containing patterns
  int it_threshold = 1000000; //number of searching steps threshold 
  int nb_sample = 0;        //number of sample
  int method = 0;           //searching method. 0: exhaustive search, 1: heuristic search (searching the largest patterns)

  ///////////////////////////////////////////////////////////
  //convert string to int
  constexpr unsigned int str2int(const char* str, int h = 0){
      return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
  }

  //split string by delimiter
  vector<string> split(string str, char delimiter) {
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;
    while(getline(ss, tok, delimiter)) {
      internal.push_back(tok);
    }
    return internal;
  }


int main(int argc, char* argv[]) {
  
  /////////////////////////////////
  //////parse input parameters///
  /////////////////////////////////
  if(argc < 3)    {
      cout<<"USAGE:"<<endl;
	  cout<<"./SSDPS [OPTION] INPUT"<<endl;
	  cout<<"OPTION:"<<endl;
      cout<<"-or: odds ratio threshold (default 1)" <<endl;
      cout<<"-rr: risk ratio threshold (default 1)"<<endl;
      cout<<"-ar: absolute risk reduction threshold (default 0)"<<endl;
//      cout<<"-pval: p-val threshold (default 0)"<<endl;
      cout<<"-min: minimal support in the 1st class (default 0%)"<<endl;
      cout<<"-max: maximal support in the 2nd class (default 100%)"<<endl;
//      cout<<"-min: minimal case support of output pattern (default 0%)"<<endl;
	  cout<<"-heuristics: mining the largest patterns (default exhaustive mining)"<<endl;
      cout<<"-iteration: number of iteration x 1,000,000 (default 1,000,000 interations)"<<endl;
      return 0;
    }else{
      input_file = argv[argc-1];
      ifstream para (input_file.c_str() , ifstream::in);
	  //read the first line of input file to find the number of case and control
	  string first_line;
      getline(para,first_line);
      vector<string> sample_size = split(first_line,' ');
      nb_case    = std::stoi(sample_size[1]);
      nb_control = std::stoi(sample_size[2]); 
      nb_sample  = nb_case+nb_control;
  	  max_control =  nb_control;
	  para.close();
	  /////////////////////////
	  float temp;	  
  	  for(int i=1; i<argc-2; ++i){
		  switch(str2int(argv[i])){
			case str2int("-or"):
			  temp = strtof(argv[i+1],&argv[i+1]);
			  if(temp!=0)
				  or_threshold = temp;
			  else
				  cout <<"parameter of "<< argv[i]<< " is invalid"<< endl;
			  break;

			case str2int("-rr"):
			  temp = strtof(argv[i+1],&argv[i+1]);
			  if(temp!=0)
				  rr_threshold = temp;
			  else
				  cout <<"parameter of "<< argv[i]<< " is invalid"<< endl;
			  break;

			case str2int("-ar"):
			  temp = strtof(argv[i+1],&argv[i+1]);
			  if(temp!=0)
				  arr_threshold = temp;
			  else
				  cout <<"parameter of "<< argv[i]<< " is invalid"<< endl;
			  break;

			case str2int("-min"):
			  temp = strtof(argv[i+1],&argv[i+1]);
			  if(temp!=0)
				  min_case = (temp/100)*nb_case;
			  else
				  cout <<"parameter of "<< argv[i]<< " is invalid"<< endl;
			  break;

			case str2int("-max"):
			  temp = strtof(argv[i+1],&argv[i+1]);
			  if(temp!=0)
				  max_control = (temp/100)*nb_control;
			  else
				  cout <<"parameter of "<< argv[i]<< " is invalid"<< endl;
			  break;
				
			  case str2int("-heuristics"):
  				  method = 1;
				  break;
			
  			  case str2int("-iteration"):
			  temp = strtof(argv[i+1],&argv[i+1]);
			  if(temp!=0)
				  it_threshold = temp*it_threshold;
			  else
				  cout <<"parameter of "<< argv[i]<< " is invalid"<< endl;	  
  			  break;  		  
		  }  
  	  }	
  }
      
	  
	/*

  //////////////////////////////////////
  /////load p-val file//////////////////
  // in case, calculating the p-value of large sample size we need create a p-val file
  //////////////////////////////////////
  std::vector< vector<string> > p_value;
  ifstream p_val_data (p_val_file.c_str(), ifstream::in);
  for(int i=0; i<nb_case+1; i++)
    {
      std::vector<string> t;
      for(int j=0; j<nb_control+1; j++)
	{
	  string line;
	  getline(p_val_data,line);
	  if(!line.empty())
	    {
	      std::vector<string> tmp = split(line,'\t');
	      //double value = std::stod(tmp[2]);
	      //p_value[i][j]=value;
	      t.push_back(tmp[2]);
	    }
	}
      p_value.push_back(t);
    }
*/


  //////////load input data into transactionTable/////
    Transaction transaction;//tt_v;
    int nb_trans = 0;
    ifstream database (input_file.c_str() , ifstream::in);

 //initial registers and masks
  int nb_registers = 0; //nb of register using
  if ( (nb_sample % nb_bits)==0 ) 
	  nb_registers = (nb_sample / nb_bits);
  else 
	  nb_registers = (nb_sample / nb_bits) + 1;

  //cout<<"nb_registers "<<nb_registers<<endl;
   
  while(database.good())    {
      //read line
      string line;
      getline(database,line);
      if(!line.empty() && (line[0]!='#'))	{
	  int bb_case = 0;
	  int bb_control = 0;
	  for(int i=0; i<nb_sample; i++)
	    if(line[i]=='1'){
	      if (i<nb_case) bb_case++;
	      else bb_control++;}
	  //double p = std::stod(p_value[bb_case][bb_control]); //get p-value of the current item
	  float p = p_value(bb_case,nb_case-bb_case,bb_control,nb_control-bb_control);
      //select items based on p-value, case and control support
	  if(p_val!=0)  {
	      if( (p <= p_val) && (p > 0) && (bb_control <= max_control) && (bb_case >= min_case) ){
		  ITEM tid;
		  tid.id = nb_trans;
		  tid.label = nb_trans;
		  transaction.tidset.push_back(tid);	  
		  //add line to tt_v
		  int kk=0;
		  __m256i mask = _mm256_set1_epi32(-1);
		  int i=0; //nb_registers
		  int n=nb_bits*nb_registers-1; //total bits reading
		  aligned_vector tmp;

		  while(i<nb_registers){
		      int j=0; //nb_chunks
		      int arr[nb_chunks];
		      while(j<nb_chunks){
			  std::bitset<32> bs;
			  int l=0;
			  while(l<32){ //transform 32 bits into int and put it into arr[j]
			      if(n<nb_sample) {int a = (int)line[kk]-48; bs[l]= a; kk++;}
			      else bs[l]=0;
			      l++;
			      n--;
			    }
			  //cout<<bs<<endl;
			  arr[j]=bs.to_ulong();
			  j++;
			}
		      __m256i x_tmp = _mm256_maskload_epi32(arr,mask);
		      tmp.push_back(x_tmp);
		      i++;
		    }
		  transaction.push_back(tmp);
		}
		 //select items based on case and control supports
	    }else {
	    if( (bb_control <= max_control) && (bb_case >= min_case) ){
		  //add tid to to tt_v
		  ITEM tid;
		  tid.id = nb_trans;
		  tid.label = nb_trans;
		  transaction.tidset.push_back(tid);	  

		  //add line to tt_v
		  int kk=0;
		  __m256i mask = _mm256_set1_epi32(-1);
		  int i=0; //nb_registers
		  int n=nb_bits*nb_registers-1; //total bits reading
		  aligned_vector tmp;

		  while(i<nb_registers){
		      int j=0; //nb_chunks
		      int arr[nb_chunks];
		      while(j<nb_chunks){
			  std::bitset<32> bs;
			  int l=0;
			  while(l<32){ //transform 32 bits into int and put it into arr[j]
			      if(n<nb_sample) {int a = (int)line[kk]-48; bs[l]= a; kk++;}
			      else bs[l]=0;
			      l++;
			      n--;
			    }
			  //cout<<bs<<endl;
			  arr[j]=bs.to_ulong();
			  j++;
			}
		      __m256i x_tmp = _mm256_maskload_epi32(arr,mask);
		      tmp.push_back(x_tmp);
		      i++;
		    }
		  transaction.push_back(tmp);
		}
	  }
	  nb_trans++;
	}
    }

  //init other values for tt_v
  transaction.nb_case = nb_case;//number of case samples
  transaction.nb_control = nb_control; //number of control samples
  transaction.nb_sample = nb_case+nb_control; //total number of samples

  //init itemset of case group
  int n = nb_bits*nb_registers;
  int pos = 0;
  for(int i=0; i<nb_registers; i++) {
    __m256i t = _mm256_setzero_si256();
    for(int j=0; j<256; j++) {
      if(pos >= (n-nb_sample) && pos < (n-nb_control) )
	 SetBit(t, pos-i*256, 1);
      pos++;
      }
    transaction.case_itemset.push_back(t);
    }

  //init itemset of control group
  pos = 0;
  for(int i=0; i<nb_registers; i++) {
    __m256i t = _mm256_setzero_si256();
    for(int j=0; j<256; j++) {
      if(pos >= (n-nb_control) )
	 SetBit(t, pos-i*256, 1);
      pos++;
      }
    transaction.control_itemset.push_back(t);
    }


 /////////////////////////////////////////
 //////discriminative pattern mining//////
 /////////////////////////////////////////
  clock_t begin = clock();
  int nb_patterns = 0;//number of output patterns
  int nb_pruning_case = 0; //number of pruning nodes in case group
  int nb_pruning_control=0; //number of pruning nodes in control group
  int nb_it = 0; //number of iterations (running steps)
  if(method==1){
      cout<<"#Heuristic mining statistically significant discriminative patterns"<<endl;
      cout<<"#size of data: "<<nb_trans<<" x "<<transaction.nb_sample<<endl;
      cout<<"#size of reduced data: "<<transaction.size()<<" x "<<transaction.nb_sample<<endl;
      cout<<"#risk thresholds (OR, RR, AR): "<<or_threshold<<", "<<rr_threshold<<", "<<arr_threshold<<endl;
      if(p_val!=0) cout<<"#p_value_threshold: "<<p_val<<endl;
      cout<<"#min case support: "<<(min_case/nb_case)*100<<"%"<<endl;
      cout<<"#max control support: "<<(max_control/nb_control)*100<<"%"<<endl;
      cout<<"#min case output: "<<(min_case_out/nb_case)*100<<"%"<<endl;
      cout<<"#stopping steps: "<<it_threshold<<endl;
      cout<<endl<<"Output:"<<endl;
      cout<<"#patterns ( % classs1 : % class2 : OR : RR : AR : CI(lci-uci) )"<<endl;
	  //start from largest tid
      for(int e=nb_case-1; e>=min_case; e--){
	  Tidset_vector p; //creat an empty transaction set (tidset)
	  //init p
	  for(int i=0; i<nb_registers; i++){
	      __m256i t = _mm256_setzero_si256();
	      p.push_back(t);
	    }
	  //expand p wich each of tid in case (e)
		expand_case_heu(p, e, or_threshold, rr_threshold, arr_threshold, min_case_out, it_threshold, nb_patterns, nb_pruning_case, transaction, nb_registers, nb_it);
	}
	 //Find all statistically significant discriminative patterns (exhaustive search)
    } else {
      cout<<"#Exhaustive mining statistically significant discriminative patterns"<<endl;
      cout<<"#size of data: "<<nb_trans<<" x "<<transaction.nb_sample<<endl;
      cout<<"#size of reduced data: "<<transaction.size()<<" x "<<transaction.nb_sample<<endl;
      cout<<"#risk thresholds (OR, RR, AR): "<<or_threshold<<", "<<rr_threshold<<", "<<arr_threshold<<endl;
      if(p_val!=0) cout<<"#p_value_threshold: "<<p_val<<endl;
      cout<<"#min case support: "<<(min_case/nb_case)*100<<"%"<<endl;
      cout<<"#max control support: "<<(max_control/nb_control)*100<<"%"<<endl;
      cout<<"#min case output: "<<(min_case_out/nb_case)*100<<"%"<<endl;
      cout<<endl<<"Output:"<<endl;
      cout<<"#patterns ( % class1 : % class2 : OR : RR : AR : CI(lci-uci) )"<<endl;
      //start from the smallest tid
      for(int e=min_case; e<transaction.nb_case; e++){
	  Tidset_vector p; //creat an empty transaction set (tidset)
	  for(int i=0; i<nb_registers; i++){
	      __m256i t = _mm256_setzero_si256();
	      p.push_back(t);
	    }
	  //expand p with each tid in case group
	  expand_case_exh(p, e, or_threshold, rr_threshold, arr_threshold, min_case_out, nb_patterns, nb_pruning_case, transaction, nb_registers, nb_pruning_control);
	}
    }
  //cout<<endl<<"#nb_prunes "<<nb_prunes<<endl;
  cout<<endl<<"#nb_patterns "<<nb_patterns<<endl;
  //cout<<endl<<"#nb_pruning_control "<<nb_pruning_control<<endl;
  clock_t end = clock();
  cout<<"#running time "<<(float)(end-begin)/CLOCKS_PER_SEC<<" s"<<endl;

  /////////////////////////////////////

  return 0;
}
