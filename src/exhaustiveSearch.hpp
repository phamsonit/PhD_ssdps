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
 *      $Id: exhaustiveSearch.hpp 2017-04-07 13:37:27Z pham $
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

#ifndef EXHAUSTIVESEARCH_HPP_
#define EXHAUSTIVESEARCH_HPP_

#include <immintrin.h>
#include <iostream>
#include <vector>
#include <bitset>
#include <stdlib.h>

#include "expand_avx.hpp"
#include "utils.hpp"

using namespace std;

/////////////////////////////////////////////////////////////////////////////
//////////////EXHAUTIVE: FIND ALL STATISTICAL SIGNIFICANT PATTERNS///////////
////////////////////////////////////////////////////////////////////////////
void print_itemset_score_exh(Tidset_vector& p, Transaction& att, float& or_threshold, float& rr_threshold, float & arr_threshold);
void expand_control_exh(Tidset_vector p, int e, float& or_threshold, float& rr_threshold, float& arr_threshold, int min_case_out, int& nb_patterns, Transaction& att, int nb_registers, int& nb_pruning_control);
void expand_case_exh(Tidset_vector p, int e, float& or_threshold, float& rr_threshold,float& arr_threshold, int min_case_out, int& nb_patterns, int& nb_pruning_case, Transaction& att, int nb_registers, int& nb_pruning_control);

#endif /* EXHAUSTIVESEARCH_HPP_ */
