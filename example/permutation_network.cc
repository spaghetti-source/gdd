//
// Permutation Network
//
#include <bits/stdc++.h>
#include "gdd.hh"
#include "permutation.hh"
#include "memoryusage.hh"
#include "util.hh"

using namespace std;

int main() {
  int n; cin >> n;
  vector<Permutation> gen = {Permutation(n)};
  for (int i = 0; i < n-1; ++i) {
    Permutation g(n);
    swap(g.p[i], g.p[i+1]);
    gen.push_back(g);
  }
  GroupDecisionDiagram gdd(gen);
  cout << "Group Order = " << gdd.groupOrder() << endl;

  int x = gdd.bot, z = gdd.top;
  for (Permutation g: gen) 
    x = gdd.cup(x, gdd.singleton(g));

  int m = 0;
  tick();
  cout << "m: #GDDSize | #TotalGDDSize | Cardinality" << endl;
  for (; ; ++m) {
    int w = gdd.cartesianProduct(x, z);
    cout << m << ": " ;
    cout << gdd.numberOfNodes(w) << " | " << gdd.node.size() << " | " << gdd.cardinality(w) << endl;
    if (w == z) break;
    z = w;
  }
  cout << "time: " << tick() << "[sec]" << endl;
  cout << "used memory: " << getPeakRSS()/(1024*1024) << "[MB]" << endl;
  return 0;
}
