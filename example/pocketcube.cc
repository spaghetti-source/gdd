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
  for (int iter = 1; iter < 15; ++iter) {
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
