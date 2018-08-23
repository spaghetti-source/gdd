#pragma once

#include <vector>
#include <queue>
#include <list>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <functional>
#include <cmath>

// #include <boost/multiprecision/cpp_int.hpp>

#include "permutation.hh"
#include "util.hh"

struct GroupDecisionDiagram {
  using BigInt = long long; 
  //using BigInt = boost::multiprecision::cpp_int;

  static const int nmax = Permutation::nmax;
  static const int rmax = ceil(nmax * log(nmax) / log(2)); 

  // g = tr[i][alpha] is a representative of G_{i-1}/G_i such that g(beta[i]) = alpha
  int n, r;
  Permutation tr[rmax][nmax], trinv[rmax][nmax];
  int beta[rmax];

  struct Node { int i, alpha, lo, hi; };
  const int bot = 0, top = 1;
  std::vector<Node> node;

  using HashTable = std::unordered_map<std::pair<int,int>,int>;
  HashTable getNodeCache[rmax][nmax]; 
  int getNode(char i, char alpha, int lo, int hi) {
    assert(alpha != beta[i]);
    if (hi == bot) return lo;
    auto key = std::make_pair(lo, hi);
    if (!getNodeCache[i][alpha].count(key)) {
      getNodeCache[i][alpha][key] = node.size();
      node.push_back({i, alpha, lo, hi});
    }
    return getNodeCache[i][alpha][key];
  }

  int member(Permutation h, int u) {
    if (h.identity() && u == top) return true;
    if (u == bot) return false;
    int alpha = h(beta[node[u].i]);
    if (alpha == beta[node[u].i]) {
      return member(h, node[u].lo); 
    } else {
      while (1) {
        if (node[u].alpha == alpha) 
          return member(tr[node[u].i][alpha].inv() * h, node[u].hi);
        if (node[u].i != node[node[u].lo].i) return false;
        u = node[u].lo;
      }
    }
  }

  int singleton(Permutation g, int i = 0) {
    int alpha;
    for (; i < r; ++i) {
      alpha = g(beta[i]);
      if (alpha != beta[i]) break;
    }
    if (i == r && g.identity()) return top;
    if (tr[i][alpha].empty())   return bot; // g is not in G
    return getNode(i, alpha, bot, singleton(trinv[i][alpha]*g, i+1));
  }

  HashTable cupCache;
  int cup(int u, int v) {
    if (u == bot) return v;
    if (v == bot) return u;
    if (u == v) return u;
    auto uv = std::minmax(u, v);
    if (!cupCache.count(uv)) {
      if (node[u].alpha > node[v].alpha) std::swap(u, v);
      if (node[u].i > node[v].i) std::swap(u, v);
      if (node[u].i == node[v].i && node[u].alpha == node[v].alpha) {
        cupCache[uv] = getNode(node[u].i, node[u].alpha, cup(node[u].lo, node[v].lo), cup(node[u].hi, node[v].hi));
      } else {
        cupCache[uv] = getNode(node[u].i, node[u].alpha, cup(node[u].lo, v), node[u].hi);
      }
    }
    return cupCache[uv];
  }

  HashTable lMCache;
  int leftMultiplication(Permutation g, int v, int i = 0) {
    if (v == bot) return bot;
    if (g.identity()) return v;
    int u = singleton(g, i);
    if (v == top) return u;
    auto uv = std::make_pair(u, v);
    if (lMCache.count(uv)) return lMCache[uv];

    while (i < node[v].i && g(beta[i]) == beta[i]) ++i;
    if (i < node[v].i) {
      int alpha = g(beta[i]);
      return getNode(i, alpha, bot, leftMultiplication(trinv[i][alpha] * g, v, i+1));
    } else {
      int child[n]; std::fill(child, child+n, 0);
      do {
        int gamma = g(node[v].alpha);
        Permutation h = trinv[i][gamma] * g * tr[i][node[v].alpha];
        child[gamma] = leftMultiplication(h, node[v].hi, i+1);
        v = node[v].lo;
      } while (node[v].i == i);
      int gamma = g(beta[i]);
      Permutation h = trinv[i][gamma] * g;
      child[gamma] = leftMultiplication(h, v, i+1);

      int v = child[beta[i]];
      for (int alpha = n-1; alpha > beta[i]; --alpha) 
        if (child[alpha]) v = getNode(i, alpha, v, child[alpha]);
      return lMCache[uv] = v;
    }
  }
  // 
  // X = g Y
  //
  // Y の元は [1, ..., i] を fix する．g は i を動かす．
  HashTable cartesianProductCache;
  int cartesianProduct(int u, int v) {
    if (u == bot || v == bot) return bot;
    if (v == top) return u;
    if (u == top) return v;
    auto uv = std::make_pair(u,v);
    if (!cartesianProductCache.count(uv)) {
      cartesianProductCache[uv] = cup(cartesianProduct(node[u].lo, v), leftMultiplication(tr[node[u].i][node[u].alpha], cartesianProduct(node[u].hi, v)));
    }
    return cartesianProductCache[uv];
  }

  std::unordered_map<int, BigInt> cardinalityCache;
  BigInt cardinality(int u) {
    if (u <= 1) return u;
    if (!cardinalityCache.count(u)) {
      BigInt a = 0;
      a += cardinality(node[u].lo) + cardinality(node[u].hi);
      cardinalityCache[u] = a;
    }
    return cardinalityCache[u];
  }
  int numberOfNodes(int u) {
    std::unordered_set<int> visited;
    std::vector<int> stack = {u};
    visited.insert(u);
    while (!stack.empty()) {
      int v = stack.back(); stack.pop_back();
      if (v <= 1) continue;
      for (int w: {node[v].lo, node[v].hi}) {
        if (!visited.count(w)) {
          visited.insert(w);
          stack.push_back(w);
        }
      }
    }
    return visited.size();
  }

  void enumerateRec(int u, Permutation pi, std::vector<Permutation> &out) {
    if (u == 0) return;
    if (u == 1) {
      out.push_back(pi);
      return;
    }
    enumerateRec(node[u].lo, pi, out);
    enumerateRec(node[u].hi, pi*tr[node[u].i][node[u].alpha], out);
  }
  std::vector<Permutation> enumerate(int u) { 
    std::vector<Permutation> out;
    enumerateRec(u, Permutation(n), out);
    sort(out.begin(), out.end());
    return out;
  }

  GroupDecisionDiagram(std::vector<Permutation> gen) {
    n = gen[0].size(); r = 0;
    node.assign(2, Node({rmax,nmax}));
    SchreierSims(gen);
  }

  BigInt groupOrder() {
    BigInt size = 1;
    for (int i = 0; i < r; ++i) {
      int c = 0;
      for (int alpha = 0; alpha < n; ++alpha) 
        if (!tr[i][alpha].empty()) ++c;
      size *= c;
    }
    return size;
  }

  // Schreier-Sims with Jerrum filter
  void SchreierSims(std::vector<Permutation> gen) {
    while (1) {
      int alpha, length; 
      std::tie(alpha, length) = longestOrbit(gen);
      if (length == 1) break;

      std::unordered_map<int, Permutation> trans;
      std::tie(trans, gen) = reduce(gen, alpha);
      for (auto x: trans) {
        tr[r][x.first] = x.second;
        trinv[r][x.first] = x.second.inv();
      }
      beta[r++] = alpha;
    }
  }
  // Compute the logest orbit of given generators
  std::pair<int, int> longestOrbit(std::vector<Permutation> gen) {
    std::vector<int> stack, visited(n);
    int max_alpha = -1, max_length = 0;
    for (int alpha = 0; alpha < n; ++alpha) {
      if (visited[alpha]) continue;
      stack.push_back(alpha);
      visited[alpha] = true;
      int length = 1;
      while (!stack.empty()) {
        int gamma = stack.back(); stack.pop_back();
        for (Permutation g: gen) {
          int delta = g(gamma);
          if (!visited[delta]) {
            stack.push_back(delta);
            visited[delta] = true;
            ++length;
          }
        }
      }
      if (length > max_length) {
        max_alpha = alpha;
        max_length = length;
      }
    }
    return std::make_pair(max_alpha, max_length);
  }
  // Compute orbit(alpha) and generators of stabilizer(alpha)
  std::pair<std::unordered_map<int, Permutation>, std::vector<Permutation>> reduce(std::vector<Permutation> gen, int alpha) {

    // construct a Schreier tree
    std::unordered_map<int, Permutation> tr;
    std::queue<int> que;
    que.push(alpha);
    tr[alpha] = Permutation(n);
    std::vector<int> order = {alpha};
    while (!que.empty()) {
      int a = que.front(); que.pop();
      for (Permutation g: gen) {
        int b = g(a);
        if (!tr.count(b)) {
          tr[b] = g * tr[a];
          order.push_back(b);
          que.push(b);
        }
      }
    }

    // construct a Jerrum tree
    struct Edge { 
      int src, dst;
      Permutation pi;
      std::list<Edge>::iterator rev;
    };
    std::vector<std::list<Edge>> adj(n);
    auto addEdge = [&](int u, Permutation g) {
      int v = g(u);
      Permutation h = g.inv();
      adj[u].push_front({u, v, g});
      adj[v].push_front({v, u, h});
      adj[u].front().rev = adj[v].begin();
      adj[v].front().rev = adj[u].begin();
    };
    auto eraseEdge = [&](std::list<Edge>::iterator e) {
      auto r = e->rev;
      adj[e->src].erase(e);
      adj[r->src].erase(r);
    };
    std::function<void(Permutation)> insert = [&](Permutation g) {
      while (!g.identity()) {
        int s, t;
        for (s = 0; (t = g(s)) == s; ++s);

        int min_u = n, min_v;
        std::vector<std::list<Edge>::iterator> next(n); 
        // can assume the graph is tree
        std::function<bool(int,int)> findPath = [&](int u, int p) {
          if (u == t) return true;
          for (auto e = adj[u].begin(); e != adj[u].end(); ++e) {
            int v = e->pi(u);
            if (v == p) continue;
            if (findPath(v, u)) {
              if (min_u > u) {
                min_u = u;
                min_v = v;
              }
              next[u] = e;
              return true;
            }
          }
          return false;
        };
        bool found = findPath(s, -1);
        addEdge(s, g); 
        if (!found) break;

        auto prod = [&](int u, int v) {
          Permutation pi(n);
          for (; u != v; u = next[u]->dst) 
            pi = next[u]->pi * pi;
          return pi;
        };
        g = prod(s, min_u) * g.inv() * prod(min_u, t);
        eraseEdge(next[min_u]);
      }
    };
    // from far to near: keep generators close to the original ones
    for (int i = order.size()-1; i >= 0; --i) {
      int q = order[i];
      Permutation &h = tr[q];
      for (Permutation g: gen) 
        insert(tr[g(q)].inv() * g * h);
    }

    std::vector<Permutation> subgen;
    for (int u = 0; u < n; ++u) {
      for (Edge e: adj[u]) {
        if (e.src < e.dst) subgen.push_back(e.pi);
      }
    }
    return std::make_pair(tr, subgen);
  }
  bool isGroupMember(Permutation g) {
    for (int i = 0; i < r; ++i) {
      int alpha = g(beta[i]);
      if (tr[i][alpha].empty()) return false;
      g = trinv[i][alpha] * g;
    }
    return g.identity();
  }

  std::string outputDot(int u) {
    std::stringstream ss;
    ss << "digraph graphname {\n";
    
    std::unordered_set<int> visited;
    std::vector<int> stack = {u};
    visited.insert(u);
    while (!stack.empty()) {
      u = stack.back(); stack.pop_back();
      for (int v: {node[u].lo, node[u].hi}) {
        if (v) ss << "  v" << u << " -> " << "v" << v << "\n";
        if (!visited.count(v)) {
          visited.insert(v);
          stack.push_back(v);
        }
      }
    }
    ss << "}\n";
    return ss.str();
  }
};
