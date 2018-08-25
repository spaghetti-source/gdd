//
// Conjugacy
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
  gen.push_back(Permutation::parse(n, "(0,1)"));
  Permutation pi(n);
  for (int i = 0; i < n; ++i) pi.p[i] = (i+1) % n;
  gen.push_back(pi);

  GroupDecisionDiagram gdd(gen);
  vector<int> gs, is;
  for (Permutation g: gen) {
    gs.push_back(gdd.singleton(g));
    is.push_back(gdd.singleton(g.inv()));
  }
  int z = gdd.fullGroup();
  for (int iter = 0; z != gdd.bot; ++iter) {
    int x = gdd.sample(z);
    while (1) {
      int prev = x;
      for (int i = 0; i < gs.size(); ++i) {
        int l = gs[i], r = is[i];
        x = gdd.cup(x, gdd.cartesianProduct(gdd.cartesianProduct(l, x), r));
      }
      if (x == prev) break;
    }
    cout << gdd.cardinality(x) << endl;
    //cout << gdd.enumerate(x) << endl;
    z = gdd.diff(z, x);
  }

  return 0;
}
