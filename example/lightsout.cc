//
// 4x4 Lights Out
//
#include <bits/stdc++.h>
#include "gdd.hh"
#include "permutation.hh"
#include "memoryusage.hh"
#include "util.hh"

using namespace std;

int main() {
  int width = 4, turn = width*width;
  int n = 2*width*width;
  vector<Permutation> gen, elem;

  /*
  for (int i = 0; i < n; ++i) {
    for (int j = i+1; j < n; ++j) {
      gen.push_back(Permutation::transposition(n, i, j));
    }
  }
  */

  int dx[] = {0,1,0,-1,0}, dy[] = {0,0,-1,0,1};
  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < width; ++j) {
      Permutation g(n);
      for (int k = 0; k < 5; ++k) {
        int x = i+dx[k], y = j+dy[k];
        if (x < 0 || y < 0 || x >= width || y >= width) continue;
        g *= Permutation::transposition(n, x+width*y, x+width*y+turn);
      }
      elem.push_back(g);
    }
  }
  gen = elem;
  GroupDecisionDiagram gdd(gen);
  cout << "Group Order = " << gdd.groupOrder() << endl;

  int z = gdd.singleton(Permutation(n));
  for (int i = 0; i < elem.size(); ++i) {
    z = gdd.cup(z, gdd.singleton(elem[i]));
  }
  int w = gdd.top;
  for (int k = 0; k < 20; ++k) {
    w = gdd.cartesianProduct(z, w);
    cout << gdd.cardinality(w) << " | " << gdd.numberOfNodes(w) << " | " << gdd.node.size() << endl;
  }

}
