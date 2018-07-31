//
// Rubik's Cube
//
#include <bits/stdc++.h>
#include "gdd.hh"
#include "permutation.hh"
#include "memoryusage.hh"
#include "util.hh"

using namespace std;

int main() {
  Permutation U = Permutation::permutation({2,4,7,1,6,0,3,5,32,33,34,11,12,13,14,15,8,9,10,19,20,21,22,23,16,17,18,27,28,29,30,31,24,25,26,35,36,37,38,39,40,41,42,43,44,45,46,47});
  Permutation D = Permutation::permutation({0,1,2,3,4,5,6,7,8,9,10,11,12,21,22,23,16,17,18,19,20,29,30,31,24,25,26,27,28,37,38,39,32,33,34,35,36,13,14,15,42,44,47,41,46,40,43,45});
  Permutation L = Permutation::permutation({16,1,2,19,4,21,6,7,10,12,15,9,14,8,11,13,40,17,18,43,20,45,22,23,24,25,26,27,28,29,30,31,32,33,5,35,3,37,38,0,39,41,42,36,44,34,46,47});
  Permutation R = Permutation::permutation({0,1,37,3,35,5,6,32,8,9,10,11,12,13,14,15,16,17,2,19,4,21,22,7,26,28,31,25,30,24,27,29,47,33,34,44,36,42,38,39,40,41,18,43,20,45,46,23});
  Permutation F = Permutation::permutation({0,1,2,3,4,24,27,29,8,9,7,11,6,13,14,5,18,20,23,17,22,16,19,21,42,25,26,41,28,40,30,31,32,33,34,35,36,37,38,39,10,12,15,43,44,45,46,47});
  Permutation B = Permutation::permutation({13,11,8,3,4,5,6,7,45,9,10,46,12,47,14,15,16,17,18,19,20,21,22,23,24,25,0,27,1,29,30,2,34,36,39,33,38,32,35,37,40,41,42,43,44,31,28,26});
  int n = U.n;
  vector<Permutation> gen = {U, D, L, R, F, B};
  
  GroupDecisionDiagram gdd(gen);
  auto order = gdd.groupOrder();
  cout << "Group Order = " << order << endl;

  int z = gdd.bot;
  for (Permutation g: gen) 
    z = gdd.cup(z, gdd.singleton(g));

  cout << "i: time | #GDDSize | #TotalGDDSize | Cardinality" << endl;
  auto prev = order; prev = 1;
  int w = gdd.top;
  for (int iter = 1; iter < 7; ++iter) {
    tick();
    w = gdd.cartesianProduct(z, w);
    cout << iter << ": " << tick() << " | ";
    auto count = gdd.cardinality(w);
    cout << gdd.numberOfNodes(w) << " | " << gdd.node.size() << " | " << count << " || " << count - prev << endl;
    if (count == order) break;
    prev = count;
  }
  return 0;
}
