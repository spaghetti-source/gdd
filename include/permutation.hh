// 
#pragma once

#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>

struct Permutation {
  static const int nmax = 100;

  int n, p[nmax];
  Permutation(int n = 0) : n(n) { // {, p(n) {
    assert(n <= nmax);
    for (int i = 0; i < n; ++i) p[i] = i;
  }
  int operator()(int i) const { 
    return p[i]; 
  }
  bool empty() const { 
    return n == 0; 
  }
  bool identity() const {
    for (int i = 0; i < n; ++i) 
      if (p[i] != i) return false;
    return true;
  }
  int size() const { return n; }

  bool operator==(const Permutation &pi) const {
    for (int i = 0; i < n; ++i) 
      if (p[i] != pi.p[i]) return false;
    return true;
  }
  bool operator<(const Permutation &pi) const {
    for (int i = 0; i < n; ++i) 
      if (p[i] != pi.p[i]) return p[i] < pi.p[i];
    return false;
  }
  Permutation &operator*=(const Permutation &pi) {
    int q[nmax];
    for (int i = 0; i < n; ++i) {
      q[i] = p[pi.p[i]];
    }
    for (int i = 0; i < n; ++i) {
      p[i] = q[i];
    }
    return *this;
  }
  Permutation inv() const {
    Permutation pi(n);
    for (int i = 0; i < n; ++i) {
      pi.p[p[i]] = i;
    }
    return pi;
  }


  static Permutation transposition(int n, int i, int j) {
    Permutation g(n);
    std::swap(g.p[i], g.p[j]);
    return g;
  }
  static Permutation permutation(std::vector<int> q) {
    Permutation g(q.size());
    for (int i = 0; i < q.size(); ++i) 
      g.p[i] = q[i];
    return g;
  }
  static Permutation parse(int n, std::string s) {
    Permutation pi(n);
    for (int i = 0; i < s.size(); ++i) {
      int j = i+1; 
      while (1) {
        if (s[j] == ')') break;
        if (s[j] == ',') s[j] = ' ';
        ++j;
      }
      std::stringstream ss(s.substr(i+1, j-i-1));
      std::vector<int> vs;
      for (int k; ss >> k; ) vs.push_back(k);
      Permutation rho(n);
      for (int i = 0; i < vs.size(); ++i) {
        rho.p[vs[i]] = vs[(i+1)%vs.size()];
      }
      pi *= rho;
      i = j;
    }
    return pi;
  }
};
Permutation operator*(Permutation x, Permutation y) { 
  return x *= y;
}
std::ostream &operator<<(std::ostream &os, Permutation g) {
  os << "[";
  for (int i = 0; i < g.n; ++i) {
    if (i > 0) os << " ";
    os << g(i);
  }
  os << "]";
  return os;
}


