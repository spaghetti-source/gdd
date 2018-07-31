//
// 4x4 Lights Out
//
#include <bits/stdc++.h>
#include "gdd.hh"
#include "permutation.hh"
#include "memoryusage.hh"
#include "util.hh"

using namespace std;

int main(int argc, char *argv[]) {
  int width, height; 
  if (argc == 2) {
    width = height = atoi(argv[1]);
  } else {
    width = atoi(argv[1]);
    height = atoi(argv[2]);
  }
  int turn = width*height;
  int n = 2*width*height;
  vector<Permutation> gen, elem;

  int dx[] = {0,1,0,-1,0}, dy[] = {0,0,-1,0,1};
  for (int i = 0; i < width; ++i) {
    for (int j = 0; j < height; ++j) {
      Permutation g(n);
      for (int k = 0; k < 5; ++k) {
        int x = i+dx[k], y = j+dy[k];
        if (x < 0 || y < 0 || x >= width || y >= height) continue;
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
  GroupDecisionDiagram::BigInt prev = 0;
  cout << "cardinality | number of nodes | total nodes" << endl;
  for (int k = 0; k < 20; ++k) {
    auto ww = w;
    w = gdd.cartesianProduct(z, w);
    if (w == ww) break;
    GroupDecisionDiagram::BigInt curr = gdd.cardinality(w);
    cout << curr - prev << " | " << gdd.numberOfNodes(w) << " | " << gdd.node.size() << endl;
    prev = curr;
  }
}
