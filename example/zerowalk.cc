#include <bits/stdc++.h>
#include "gdd.hh"
#include "permutation.hh"
#include "memoryusage.hh"
#include "util.hh"

using namespace std;

int main(int argc, char *argv[]) {
  int n = atoi(argv[1]), k = atoi(argv[2]);
  vector<Permutation> gen;
  gen.push_back(Permutation::parse(n, "(0,1)"));
  Permutation pi(n);
  for (int i = 0; i < n; ++i) pi.p[i] = (i+1) % n;
  gen.push_back(pi);
  GroupDecisionDiagram gdd(gen);

  struct Edge {
    int src, dst;
    Permutation g;
  };
  int N = k*k;
  vector<vector<Edge>> adj(N);
  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < k; ++j) {
      int dx[] = {1,0,-1,0}, dy[] = {0,1,0,-1};
      for (int p = 0; p < 4; ++p) {
        int I = i + dx[p], J = j + dy[p];
        if (I < 0 || I >= k || J < 0 || J >= k) continue;
        int src = i*k + j;
        int dst = I*k + J;
        Permutation pi = gen[rand() % 2];
        adj[src].push_back({src, dst, pi});
      }
    }
  }
  int start = 0, end = N-1;

  vector<int> curr(N);
  curr[start] = 1; // unit
  int dist = 0;
  for (int i = 0; ; ++i) {
    cout << i << ": ";
    for (int u = 0; u < N; ++u) {
      cout << gdd.cardinality(curr[u]) << " ";
    }
    cout << endl;
    if (gdd.member(Permutation(n), curr[end])) break;

    vector<int> next = curr;
    for (int u = 0; u < N; ++u) {
      for (Edge &e: adj[u]) {
        next[e.dst] = gdd.cup(next[e.dst], gdd.leftMultiplication(e.g, curr[e.src]));
      }
    }
    if (curr == next) break;
    curr = next;
    ++dist;
  }
  cout << "distance = " << dist << endl;
  return 0;
}
