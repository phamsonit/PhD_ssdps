# SSDPS - Statistically Significant Discriminative Patterns Search
SSDPS is an algorithm to mine statistically significant discrimininative patterns
in two-class datasets. It allows to exhaustively discover patterns which satisfy
both discriminative scores and confidence intervals or using heuristic strategies to
mine the largest statistically significant discriminative patterns.
Please refer to the full paper at https://arxiv.org/abs/1906.01581 for more details.


# Usage: #

## Compile ssdps ##
Using the follow command-line to compile the program:

`g++ -std=c++11 -mavx2 src/*.cpp -o SSDPS`

### Running SSDPS

SSDPS is run using the following command-line format:

`./SSDPS [options] INPUT `

The options, input file format and output are described below.

#### Options:

`-or <n>`

Odds ratio threshold. Default is 0.

`-gr <n>`

Grown rate support threshold. Default is 0.

`-ds <n>`

Difference support threshold. Default is 0.

`-p-value <n>`

p-value threshold. Default is 0.

`-min <n>`

Minimal support of the pattern in the first class. Default is 0.0%	

`-max <n>`

Maximal support of the pattern in the 2nd class. Default is 1.0%

`-heuristic`

Using heuristic search.

`-iteration <n>`

Number of searching iterations x 1,000,000. Default is 1x1,000,000. This option is only used with -heuristic option. 



#### Input data
The input data of SSDPS can be stored in a plain text file. The folowing example shows an input data including 16 transactions (8 transactions of 1st class, 8 transactions of 2nd class) and 10 items.

\# 8   8

1110010011100100

1110101011101010

1110010111101010

0101101001101011

0001011101101011

0110101101111001

1110011101111001

0111100101011100

1001101101011100

0101110001011100

The format of this file as follows:
- The first row contains the number of transactions of the 1st class and the 2nd class, respectively. Note that the character \# is a sign indicated the first line. These numbers are separated by a space.
- The following lines present the set of items. Each line corresponds to an item id. In a line, the value of 1 at the column i^th means that the corresponding item id occurs in the i^th transaction id.
- Columns correspond to transaction ids. Transaction ids of the first class are presented first. For examle, the first 8 columns present the transaction ids of the 1st class, and the last 8 columns correspond to the transaction ids of the 2nd class.

#### Output format

Each line of the output file presents a discriminative pattern.
For example:
1,3,6,10,22 (15% : 4% : 4.92857 : 3.75 : 0.22 : 1.50337-16.1576)

This pattern can be interpreted as follow:
- 1,3,6,10,22 : pattern (the set of items, or itemset),
- 15%: the support of the pattern in the 1st class,
- 4%: the support of the pattern in the 2nd class,
- 4.92857: the value of odd ratio of supports,
- 3.75: the value of grow rate supports,
- 0.22: the difference of supports,
- 1.50337-16.1576: lower confidence interval - uper confidence interval.

For questions, bug or reports, please mail to: hoang-s.pham@uclouvain.be

## Authors

Pham Hoang Son: hoang.s.pham@uclouvain.be

Alexandre Termier: Alexandre.Termier@irisa.fr

Dominique Lavenier: Dominique.Lavenier@irisa.fr

## Copy right

This program is free software; you can redistribute it and/or modify it under the terms of
the GNU Lesser General Public License as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

