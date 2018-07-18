#pragma once

#include <cassert>
#include <ctime>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

double tick() {
  static clock_t oldtick;
  clock_t newtick = std::clock();
  double diff = 1.0*(newtick - oldtick) / CLOCKS_PER_SEC;
  oldtick = newtick;
  return diff;
}


template <class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  os << "[";
  for (int i = 0; i < v.size(); os << v[i++]) 
    if (i > 0) os << " ";
  os << "]";
  return os;
}

template <class S, class T>
std::ostream &operator<<(std::ostream &os, const std::pair<S, T> &x) {
  os << "(" << x.first << "," << x.second << ")";
  return os;
}

std::ostream &operator<<(std::ostream &os, const __int128_t &value) {
  if (std::ostream::sentry(os)) {
    __uint128_t tmp = value < 0 ? -value : value;
    char buffer[64];
    char *d = std::end(buffer);
    do {
      --d;
      *d = "0123456789"[tmp % 10];
      tmp /= 10;
    } while (tmp != 0);
    if (value < 0) {
      --d;
      *d = '-';
    }
    int len = std::end(buffer) - d;
    if (os.rdbuf()->sputn(d, len) != len) {
      os.setstate(std::ios_base::badbit);
    }
  }
  return os;
}

template <class T>
inline void hash_combine(std::size_t & seed, const T & v) {
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}
namespace std {
  template<typename S, typename T> struct hash<pair<S, T>> {
    inline size_t operator()(const pair<S, T> & v) const {
      size_t seed = 0;
      ::hash_combine(seed, v.first);
      ::hash_combine(seed, v.second);
      return seed;
    }
  };
}
