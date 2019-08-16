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
 *      $Id: utils.cpp 2017-04-07 13:37:27Z pham $
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
#include <iostream>
#include "utils.hpp"

using namespace std;


//TODO: p-value: Pearson extract test cannot calculate large number !
float p_value(int a, int b, int c, int d)
{
	int n = a+b+c+d;
	float p = ( std::tgamma(a+b+1)*std::tgamma(c+d+1)*std::tgamma(a+c+1)*std::tgamma(b+d+1) ) /
		      ( std::tgamma(a+1)*std::tgamma(b+1)*std::tgamma(c+1)*std::tgamma(d+1)*std::tgamma(n+1) );
	return p;
}

//return odd ratio
float odd_ratio(int a, int b, int c, int d)
{
  float odd1;
  float odd2;
  float odd;

  if(b!=0) odd1 = a/b; else odd1=0;
  if(d!=0) odd2 = c/d; else odd2=0;
  //if(odd2!=0) odd = odd1/odd2; else odd=999;
  odd = odd1/odd2;

  return float(a*d) / float(b*c);
}

//return upper confidence interval of odds ratio
float UCI(float odd, int a, int b, int c, int d)
{
  //return exp( log(odd) + 1.96 * sqrt(1/a+1/b+1/c+1/d));
  return exp( log(odd) + 1.96 * sqrt(1/float(a) + 1/float(b) + 1/float(c) + 1/float(d)) );
}

//return lower confidence interval of odds ratio
float LCI(float odd, int a, int b, int c, int d)
{
  return exp( log(odd) - 1.96 * sqrt(1/float(a) + 1/float(b) + 1/float(c) + 1/float(d)) );
}

//////////////////////////////////////////////////
//return risk ratio
float risk_ratio(int a, int b, int c, int d)
{
  float case_support = (float)a/(a+b);
  float control_support = (float)c/(c+d);
  return (case_support/control_support);
}

//return upper confidence interval of risk ratio
float R_UCI(float odd, int a, int b, int c, int d)
{
  return exp( log(odd) + 1.96 * sqrt(1/float(a) - 1/float(a+b) + 1/float(c) - 1/float(c+d)) );
}

//return lower confidence interval of odds ratio
float R_LCI(float odd, int a, int b, int c, int d)
{
  return exp( log(odd) - 1.96 * sqrt(1/float(a) - 1/float(a+b) + 1/float(c) - 1/float(c+d)) );
}

/////////////////////////////////////////////////////////////////
//return risk difference
float difference_risk(int a, int b, int c, int d)
{
  float case_support = (float)a/(a+b);
  float control_support = (float)c/(c+d);

  return (case_support - control_support);

}

////////////////////////////////////////////////////////////////////////////////////
//return chi^2
/// chi^2 function for itemsets
//inline int convex_function(int posTot, int pos, int negTot, int neg, int precision) {
//inline int X2(int X, int a, int Y, int b)
float chi2(int a, int b, int c, int d)
 {
   int   X = a+b;
   int   Y = c+d;
    // Too big: float calc = ((Y*a - X*b)*(Y*a - X*b)*(X+Y)) / (float)((a+b)*(X+Y-a-b)*X*Y);
    float yaminxb = Y*a - X*c;
    float one = yaminxb / ((a+c)*(X+Y-a-c));
    if (isnan(one)) return 0; // possible division by 0, return 0
    float two = yaminxb / (X*Y);

    return one*two*(X+Y);
}

////////////////////////////////////////////////////////
/// information gain function for itemsets
float nlogn ( float n ) {
    // Calculate -n*log(n) with 0 if n=0
    if ( !n )
        return 0;
    return -log(n) * n;
}
//inline int convex_function(int posTot, int pos, int negTot, int neg, int precision) {
//return info_gain
float info_gain(int a, int b, int c, int d)
{
  int posTot = a+b;
  int negTot = c+d;
  int pos = a;
  int neg = c;

  float tot = posTot + negTot;
  float base = nlogn(posTot/tot) + nlogn(negTot/tot);
    
  float covered = pos + neg; 
  float notCovered = tot - covered; 
  float calc = covered    * ( nlogn(pos/covered) + nlogn(neg/covered) ) +
    notCovered * ( nlogn((posTot-pos)/notCovered) + nlogn((negTot-neg)/notCovered) );
  if (isnan(calc)) calc = -1;
  return (base - calc/tot);
}
