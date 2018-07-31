//
// Pancake Problem
//
// https://en.wikipedia.org/wiki/Pancake_sorting
//
#include <bits/stdc++.h>
#include "gdd.hh"
#include "permutation.hh"
#include "memoryusage.hh"
#include "util.hh"

using namespace std;

int main(int argc, char *argv[]) {
  int n = atoi(argv[1]);
  vector<Permutation> gen;
  for (int k = 1; k < n; ++k) {
    Permutation g(n);
    for (int i = 0; i <= k; ++i) 
      g.p[i] = k-i;
    gen.push_back(g);
  }
  
  GroupDecisionDiagram gdd(gen);
  gdd.node.reserve(4000000);
  auto order = gdd.groupOrder();
  cout << "Group Order = " << order << endl;

  int z = gdd.singleton(Permutation(n));
  for (Permutation g: gen) 
    z = gdd.cup(z, gdd.singleton(g));

  cout << "i: time | #GDDSize | #TotalGDDSize | Cardinality" << endl;
  auto prev = order; prev = 1;
  int w = gdd.top;
  for (int iter = 1; ; ++iter) {
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
