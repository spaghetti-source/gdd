# Group Decision Diagram

Group Decision Diagram (GDD) is a data structure to maintain a subset of a permutation group.
This is a C++ implementation of GDD.

## Install

This is a header-only library. Please copy the files in ./include

## Usage

````
#include <bits/stdc++.h>
#include "gdd.hh"
#include "permutation.hh"
#include "memoryusage.hh"
#include "util.hh"

using namespace std;

int main() {
  Permutation R = Permutation::parse(24, "(12,13,15,14)(9,1,18,21)(11,3,16,23)");
  Permutation U = Permutation::parse(24, "(0,1,3,2)(8,4,16,12)(9,5,17,13)");
  Permutation F = Permutation::parse(24, "(8,9,11,10)(2,12,21,7)(3,14,20,5)");

  int n = U.n;
  vector<Permutation> gen = {R, U, F};
  
  GroupDecisionDiagram gdd(gen);
  auto order = gdd.groupOrder();
  cout << "Group Order = " << order << endl;

  int z = gdd.top;
  for (Permutation g: gen) {
    z = gdd.cup(z, gdd.singleton(g));
    z = gdd.cup(z, gdd.singleton(g.inv()));
  }

  cout << "i: time | #GDDSize | #TotalGDDSize | Cardinality" << endl;
  auto prev = order; prev = 1;
  int w = gdd.top;

  time_t begin = clock();
  for (int iter = 1; iter < 15; ++iter) {
    w = gdd.cartesianProduct(z, w);
    time_t end = clock();
    double diff = 1.0 * (begin - end) / CLOCKS_PER_SEC;
    cout << iter << ": " << diff << " | ";
    auto count = gdd.cardinality(w);
    cout << gdd.numberOfNodes(w) << " | " << gdd.node.size() << " | " << count << " || " << count - prev << endl;
    if (count == order) break;
    prev = count;
  }
  return 0;
}
````

## Performance

### Rubik Cube's puzzle:

|         |        |         |           |            |               |
| :-----: | -----: | ------: | --------: | ---------: | ------------: |
|  \(k\)  |      3 |       4 |         5 |          6 |             7 |
| \#Conf. |  1,068 |  10,011 |    93,840 |    878,880 |     8,221,632 |
| GDD (size) |  6,537 |  58,850 |   476,497 |  3,641,049 |    26,721,270 |
| GDD (time) |  0.02s |   0.30s |     3.54s |     35.04s |        338.2s |
| piDD (size) | 26,166 | 247,770 | 2,319,674 | 21,062,444 |             — |
| piDD (time) |  0.11s |   1.80s |     23.4s |     247.3s | \(\le\) 7200s |
| rhoDD (size)  | 36,799 | 331,549 | 2,978,733 | 26,715,003 |             — |
| rhoDD (time) |  0.51s |   6.36s |     58.5s |     572.9s | \(\le\) 7200s |

piDD, rhoDD are existing data structure for the same purpose.

## Reference 

- Takanori Maehara and Yuma Inoue (2019): "Group Decision Diagram (GDD): A Compact Representation for Permutations", in Proceedings of the 33rd AAAI Conference on Artificial Intelligence (AAAI'19), Honolulu, Hawaii, United States, January 27--February 01, 2019.
