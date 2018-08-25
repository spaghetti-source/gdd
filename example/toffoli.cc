//
// Toffoli gate synthesis
//
// Toffoli gate (P, N, t) is a boolean function that determined by
//   Gate(x)_i = x_i for all i \neq t
//   Gate(x)_t = lnot x_t iff x_p = 1 for all p in P and x_n = 0 for all n in N.
// Each gate is a permutation on 2^n, and they generate G_{2^n}.
//
// The problem is to find the diameter of the Cayley graph of the gates.
// 

#include <bits/stdc++.h>
#include "gdd.hh"
#include "permutation.hh"
#include "memoryusage.hh"
#include "util.hh"

using namespace std;

vector<vector<int>> binary;
map<vector<int>,int> idx;

void makeAllBits(int n) {
  vector<int> x(n);
  function<void(int)> rec = [&](int i) {
    if (i == n) {
      idx[x] = binary.size();
      binary.push_back(x);
      return;
    }
    rec(i+1);
    x[i] = 1;
    rec(i+1);
    x[i] = 0;
  };
  rec(0);
}

Permutation Toffoli(vector<int> P, vector<int> N, int t) {
  int m = binary.size();
  Permutation pi(m);
  for (int i = 0; i < m; ++i) {
    vector<int> x = binary[i];
    bool check = true;
    for (int j = 0; j < x.size(); ++j) {
      if (P[j] == 1 && x[j] == 0) check = false;
      if (N[j] == 1 && x[j] == 1) check = false;
    }
    if (check) x[t] ^= 1;
    pi.p[i] = idx[x];
  }
  return pi;
}
vector<Permutation> makeAllToffoli(int n) {
  vector<Permutation> gen;
  for (int t = 0; t < n; ++t) {
    vector<int> arr;
    for (int i = 0; i < n; ++i)
      if (i != t) arr.push_back(i);
    vector<int> P(n), N(n);
    function<void(int)> tripartite = [&](int i) {
      if (i == arr.size()) {
        Permutation pi = Toffoli(P, N, t);
        if (!pi.identity()) {
          if (find(gen.begin(), gen.end(), pi) == gen.end()) 
            gen.push_back(pi);
        }
        return;
      }
      int j = arr[i];
      tripartite(i+1);
      P[j] = 1;
      tripartite(i+1);
      P[j] = 0;
      N[j] = 1;
      tripartite(i+1);
      N[j] = 0;
    };
    tripartite(0);
  }
  return gen;
}

int main(int argc, char *argv[]) {
  int n = atoi(argv[1]);

  makeAllBits(n);
  vector<Permutation> gen = makeAllToffoli(n);
  cout << gen.size() << endl;
  int m = gen[0].n;

  /*
  vector<Permutation> pidd;
  for (int i = 0; i < m; ++i) {
    for (int j = i+1; j < m; ++j) {
      Permutation pi(m);
      swap(pi.p[i], pi.p[j]);
      pidd.push_back(pi);
    }
  }
  GroupDecisionDiagram gdd(pidd);
  */
  GroupDecisionDiagram gdd(gen);

  auto order = gdd.groupOrder();
  cout << order << endl;

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
    double diff = 1.0 * (end - begin) / CLOCKS_PER_SEC;
    cout << iter << ": " << diff << " | ";
    auto count = gdd.cardinality(w);
    cout << gdd.numberOfNodes(w) << " | " << gdd.node.size() << " | " << count << " || " << count - prev << endl;
    if (count == prev) break;
    prev = count;
  }
  return 0;
}
