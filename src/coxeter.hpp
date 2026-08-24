#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

// Finite Coxeter groups via the geometric representation.
// H3/H4 use Z[phi] with phi^2 = phi+1; Weyl types use ordinary integers.

static constexpr int COX_RMAX = 6;

struct Phi {
  int64_t a = 0;  // a + b*phi
  int64_t b = 0;
};

inline Phi ph_add(Phi x, Phi y) { return {x.a + y.a, x.b + y.b}; }
inline Phi ph_sub(Phi x, Phi y) { return {x.a - y.a, x.b - y.b}; }
inline Phi ph_mul(Phi x, Phi y) {
  return {x.a * y.a + x.b * y.b, x.a * y.b + x.b * y.a + x.b * y.b};
}
inline Phi ph_scale(int k, Phi x) { return {k * x.a, k * x.b}; }
inline bool ph_eq(Phi x, Phi y) { return x.a == y.a && x.b == y.b; }
inline bool ph_zero(Phi x) { return x.a == 0 && x.b == 0; }

// a + b*phi > 0 with phi = (1+sqrt(5))/2, i.e. 2a+b + b*sqrt(5) > 0.
inline bool ph_pos(Phi x) {
  const __int128 s = 2 * static_cast<__int128>(x.a) + x.b;
  const __int128 b = x.b;
  if (b == 0) return s > 0;
  const __int128 s2 = s * s;
  const __int128 fiveb2 = 5 * b * b;
  if (b > 0) return s >= 0 || s2 < fiveb2;
  return s > 0 && s2 > fiveb2;
}

struct CVec {
  std::array<Phi, COX_RMAX> c{};
  int rank = 0;
};

inline bool cvec_eq(const CVec& a, const CVec& b) {
  if (a.rank != b.rank) return false;
  for (int i = 0; i < a.rank; ++i)
    if (!ph_eq(a.c[static_cast<size_t>(i)], b.c[static_cast<size_t>(i)])) return false;
  return true;
}

inline bool cvec_pos(const CVec& v) {
  for (int i = 0; i < v.rank; ++i) {
    if (ph_zero(v.c[static_cast<size_t>(i)])) continue;
    return ph_pos(v.c[static_cast<size_t>(i)]);
  }
  return false;
}

inline uint64_t cvec_hash(const CVec& v) {
  uint64_t h = 1469598103934665603ULL;
  for (int i = 0; i < v.rank; ++i) {
    h ^= static_cast<uint64_t>(v.c[static_cast<size_t>(i)].a);
    h *= 1099511628211ULL;
    h ^= static_cast<uint64_t>(v.c[static_cast<size_t>(i)].b);
    h *= 1099511628211ULL;
  }
  return h;
}

struct CMat {
  std::array<CVec, COX_RMAX> col{};
  int rank = 0;
};

inline uint64_t cmat_hash(const CMat& m) {
  uint64_t h = 0xcbf29ce484222325ULL;
  for (int j = 0; j < m.rank; ++j) {
    h ^= cvec_hash(m.col[static_cast<size_t>(j)]);
    h *= 0x100000001b3ULL;
  }
  return h;
}

inline bool cmat_eq(const CMat& a, const CMat& b) {
  if (a.rank != b.rank) return false;
  for (int j = 0; j < a.rank; ++j)
    if (!cvec_eq(a.col[static_cast<size_t>(j)], b.col[static_cast<size_t>(j)])) return false;
  return true;
}

inline CVec cmat_apply(const CMat& m, const CVec& v) {
  CVec out;
  out.rank = m.rank;
  for (int j = 0; j < m.rank; ++j) {
    if (ph_zero(v.c[static_cast<size_t>(j)])) continue;
    for (int i = 0; i < m.rank; ++i) {
      out.c[static_cast<size_t>(i)] = ph_add(
          out.c[static_cast<size_t>(i)],
          ph_mul(v.c[static_cast<size_t>(j)], m.col[static_cast<size_t>(j)].c[static_cast<size_t>(i)]));
    }
  }
  return out;
}

inline CVec e_simple(int rank, int i) {
  CVec v;
  v.rank = rank;
  v.c[static_cast<size_t>(i)] = {1, 0};
  return v;
}

inline CMat id_mat(int rank) {
  CMat m;
  m.rank = rank;
  for (int i = 0; i < rank; ++i) m.col[static_cast<size_t>(i)] = e_simple(rank, i);
  return m;
}

struct CoxeterGroup {
  std::string name;
  int rank = 0;
  int n_elt = 0;
  int n_roots = 0;
  int w0_len = 0;
  std::array<std::array<Phi, COX_RMAX>, COX_RMAX> A{};
  std::vector<CMat> elt;
  std::vector<uint8_t> length;
  std::vector<std::array<uint32_t, COX_RMAX>> rs;
  std::vector<uint64_t> inv_lo, inv_hi;
  std::vector<uint32_t> inv;  // inv[w] = w^{-1}
  std::vector<CVec> proot;
  int row_words = 0;
  std::vector<uint64_t> below;

  uint32_t mul(uint32_t w, int s) const { return rs[w][static_cast<size_t>(s)]; }
  bool rdesc(uint32_t w, int s) const { return length[mul(w, s)] < length[w]; }
  int first_rdesc(uint32_t w) const {
    for (int s = 0; s < rank; ++s)
      if (rdesc(w, s)) return s;
    return -1;
  }
  int last_rdesc(uint32_t w) const {
    for (int s = rank - 1; s >= 0; --s)
      if (rdesc(w, s)) return s;
    return -1;
  }
  bool le(uint32_t u, uint32_t v) const {
    const uint64_t word = below[static_cast<size_t>(v) * static_cast<size_t>(row_words) + (u >> 6)];
    return (word >> (u & 63u)) & 1ull;
  }
  uint32_t word(const int* gens, int n) const {
    uint32_t w = 0;
    for (int i = 0; i < n; ++i) w = mul(w, gens[i]);
    return w;
  }
  // First-right-descent peeling: w * s0 * s1 * ... = 1, so the printed word is reduced
  // for w^{-1}. Certificate lines use this encoding; multiplying a cert word recovers w^{-1}.
  std::string show(uint32_t w) const {
    if (w == 0) return "e";
    std::string o;
    while (length[w]) {
      const int s = first_rdesc(w);
      if (!o.empty()) o += ',';
      o += std::to_string(s + 1);
      w = mul(w, s);
    }
    return o;
  }
};

inline CVec reflect(const CoxeterGroup& g, int s, const CVec& x) {
  Phi pair{0, 0};
  for (int j = 0; j < g.rank; ++j)
    pair = ph_add(pair, ph_mul(g.A[static_cast<size_t>(s)][static_cast<size_t>(j)],
                              x.c[static_cast<size_t>(j)]));
  CVec y = x;
  y.c[static_cast<size_t>(s)] = ph_sub(y.c[static_cast<size_t>(s)], pair);
  return y;
}

inline CMat s_mat(const CoxeterGroup& g, int s) {
  CMat m;
  m.rank = g.rank;
  for (int j = 0; j < g.rank; ++j) m.col[static_cast<size_t>(j)] = reflect(g, s, e_simple(g.rank, j));
  return m;
}

inline CMat cmat_mul(const CMat& a, const CMat& b) {
  CMat c;
  c.rank = a.rank;
  for (int j = 0; j < a.rank; ++j) c.col[static_cast<size_t>(j)] = cmat_apply(a, b.col[static_cast<size_t>(j)]);
  return c;
}

struct MatKey {
  CMat m;
  bool operator==(const MatKey& o) const { return cmat_eq(m, o.m); }
};
struct MatKeyHash {
  size_t operator()(const MatKey& k) const { return static_cast<size_t>(cmat_hash(k.m)); }
};

inline void fill_cartan(CoxeterGroup& g, const int m[COX_RMAX][COX_RMAX], bool htype) {
  for (int i = 0; i < g.rank; ++i) {
    for (int j = 0; j < g.rank; ++j) {
      if (i == j) {
        g.A[static_cast<size_t>(i)][static_cast<size_t>(j)] = {2, 0};
        continue;
      }
      const int mij = m[i][j];
      if (mij <= 2) g.A[static_cast<size_t>(i)][static_cast<size_t>(j)] = {0, 0};
      else if (mij == 3) g.A[static_cast<size_t>(i)][static_cast<size_t>(j)] = {-1, 0};
      else if (mij == 4) {
        g.A[static_cast<size_t>(i)][static_cast<size_t>(j)] = {0, 0};
      } else if (mij == 5) {
        if (!htype) throw std::runtime_error("m=5 needs H-type");
        g.A[static_cast<size_t>(i)][static_cast<size_t>(j)] = {0, -1};  // -phi
      } else {
        throw std::runtime_error("bad m_ij");
      }
    }
  }
}

// A_{long,short}=-1, A_{short,long}=-2 so that the product is 2 (m=4).
inline void set_double_bond(CoxeterGroup& g, int long_idx, int short_idx) {
  g.A[static_cast<size_t>(long_idx)][static_cast<size_t>(short_idx)] = {-1, 0};
  g.A[static_cast<size_t>(short_idx)][static_cast<size_t>(long_idx)] = {-2, 0};
}

inline bool le_recursive(const CoxeterGroup& g, uint32_t u, uint32_t v) {
  if (u == v) return true;
  if (g.length[u] > g.length[v]) return false;
  const int s = g.first_rdesc(v);
  if (s < 0) return false;
  const uint32_t vs = g.mul(v, s);
  if (g.rdesc(u, s)) return le_recursive(g, g.mul(u, s), vs);
  return le_recursive(g, u, vs);
}

inline uint32_t inverse_of(const CoxeterGroup& g, uint32_t w) {
  int word[64];
  int n = 0;
  uint32_t x = w;
  while (g.length[x]) {
    const int s = g.first_rdesc(x);
    if (s < 0 || n >= 64) throw std::runtime_error("inverse_of: bad word");
    word[n++] = s;
    x = g.mul(x, s);
  }
  uint32_t inv = 0;
  for (int i = 0; i < n; ++i) inv = g.mul(inv, word[i]);
  return inv;
}

inline void build_inversions(CoxeterGroup& g) {
  g.proot.clear();
  auto consider = [&](const CVec& v) {
    if (!cvec_pos(v)) return;
    for (const CVec& p : g.proot)
      if (cvec_eq(p, v)) return;
    g.proot.push_back(v);
  };
  for (int i = 0; i < g.rank; ++i) consider(e_simple(g.rank, i));
  for (size_t k = 0; k < g.proot.size(); ++k) {
    if (g.proot.size() > 128) throw std::runtime_error(g.name + ": root generation diverged");
    for (int s = 0; s < g.rank; ++s) consider(reflect(g, s, g.proot[k]));
  }
  g.n_roots = static_cast<int>(g.proot.size());
  if (g.n_roots != g.w0_len)
    throw std::runtime_error(g.name + ": n_roots=" + std::to_string(g.n_roots) +
                             " != ell(w0)=" + std::to_string(g.w0_len));
  g.inv.resize(static_cast<size_t>(g.n_elt));
  for (int w = 0; w < g.n_elt; ++w)
    g.inv[static_cast<size_t>(w)] = inverse_of(g, static_cast<uint32_t>(w));
  for (int w = 0; w < g.n_elt; ++w)
    if (g.inv[g.inv[static_cast<size_t>(w)]] != static_cast<uint32_t>(w))
      throw std::runtime_error(g.name + ": inverse is not an involution at " + std::to_string(w));
  g.inv_lo.assign(static_cast<size_t>(g.n_elt), 0);
  g.inv_hi.assign(static_cast<size_t>(g.n_elt), 0);
  for (int w = 0; w < g.n_elt; ++w) {
    uint64_t lo = 0, hi = 0;
    const CMat& winv = g.elt[g.inv[static_cast<size_t>(w)]];
    for (int r = 0; r < g.n_roots; ++r) {
      const CVec im = cmat_apply(winv, g.proot[static_cast<size_t>(r)]);
      if (cvec_pos(im)) continue;
      if (r < 64) lo |= 1ULL << r;
      else hi |= 1ULL << (r - 64);
    }
    g.inv_lo[static_cast<size_t>(w)] = lo;
    g.inv_hi[static_cast<size_t>(w)] = hi;
    const int pop = __builtin_popcountll(lo) + __builtin_popcountll(hi);
    if (pop != static_cast<int>(g.length[static_cast<size_t>(w)]))
      throw std::runtime_error(g.name + ": left-inv != length at " + std::to_string(w));
  }
}

inline void build_bruhat(CoxeterGroup& g) {
  const int n = g.n_elt;
  g.row_words = (n + 63) / 64;
  g.below.assign(static_cast<size_t>(n) * static_cast<size_t>(g.row_words), 0);
  auto setbit = [&](uint32_t u, uint32_t v) {
    g.below[static_cast<size_t>(v) * static_cast<size_t>(g.row_words) + (u >> 6)] |= 1ull
                                                                                    << (u & 63u);
  };
  std::vector<std::vector<int>> by_len(static_cast<size_t>(g.w0_len) + 1);
  for (int i = 0; i < n; ++i) by_len[g.length[static_cast<size_t>(i)]].push_back(i);
  setbit(0, 0);
  for (int L = 1; L <= g.w0_len; ++L) {
    for (int v : by_len[static_cast<size_t>(L)]) {
      const int s = g.first_rdesc(static_cast<uint32_t>(v));
      const uint32_t vs = g.mul(static_cast<uint32_t>(v), s);
      for (int u = 0; u < n; ++u) {
        if (u == v) {
          setbit(static_cast<uint32_t>(u), static_cast<uint32_t>(v));
          continue;
        }
        if (g.length[static_cast<size_t>(u)] > L) continue;
        const bool ans = g.rdesc(static_cast<uint32_t>(u), s)
                             ? g.le(g.mul(static_cast<uint32_t>(u), s), vs)
                             : g.le(static_cast<uint32_t>(u), vs);
        if (ans) setbit(static_cast<uint32_t>(u), static_cast<uint32_t>(v));
      }
    }
  }
}

inline CoxeterGroup build_coxeter(const std::string& name, int rank, const int m[COX_RMAX][COX_RMAX],
                                  bool htype, int expect_w, int expect_w0, int bond_long = -1,
                                  int bond_short = -1) {
  CoxeterGroup g;
  g.name = name;
  g.rank = rank;
  fill_cartan(g, m, htype);
  if (bond_long >= 0) set_double_bond(g, bond_long, bond_short);
  std::array<CMat, COX_RMAX> sm{};
  for (int s = 0; s < rank; ++s) sm[static_cast<size_t>(s)] = s_mat(g, s);

  std::unordered_map<MatKey, uint32_t, MatKeyHash> id_of;
  g.elt.push_back(id_mat(rank));
  g.length.push_back(0);
  g.rs.push_back({});
  id_of[MatKey{g.elt[0]}] = 0;

  const int cap = expect_w ? expect_w + 1 : 50000;
  for (uint32_t i = 0; i < g.elt.size(); ++i) {
    if (static_cast<int>(g.elt.size()) > cap)
      throw std::runtime_error(name + ": group enumeration exceeded cap " + std::to_string(cap));
    for (int s = 0; s < rank; ++s) {
      CMat nxt = cmat_mul(g.elt[i], sm[static_cast<size_t>(s)]);
      MatKey key{nxt};
      uint32_t id;
      const auto it = id_of.find(key);
      if (it == id_of.end()) {
        id = static_cast<uint32_t>(g.elt.size());
        id_of[key] = id;
        g.elt.push_back(nxt);
        g.length.push_back(static_cast<uint8_t>(g.length[i] + 1));
        g.rs.push_back({});
      } else {
        id = it->second;
      }
      g.rs[i][static_cast<size_t>(s)] = id;
    }
  }
  g.n_elt = static_cast<int>(g.elt.size());
  g.w0_len = 0;
  for (uint8_t L : g.length) g.w0_len = std::max(g.w0_len, static_cast<int>(L));
  if (expect_w && g.n_elt != expect_w)
    throw std::runtime_error(name + ": |W|=" + std::to_string(g.n_elt) + " want " +
                             std::to_string(expect_w));
  if (expect_w0 && g.w0_len != expect_w0)
    throw std::runtime_error(name + ": ell(w0)=" + std::to_string(g.w0_len) + " want " +
                             std::to_string(expect_w0));
  build_inversions(g);
  build_bruhat(g);
  return g;
}

inline CoxeterGroup make_H3() {
  int m[COX_RMAX][COX_RMAX] = {};
  m[0][1] = m[1][0] = 5;
  m[1][2] = m[2][1] = 3;
  m[0][2] = m[2][0] = 2;
  return build_coxeter("H3", 3, m, true, 120, 15);
}

inline CoxeterGroup make_H4() {
  int m[COX_RMAX][COX_RMAX] = {};
  m[0][1] = m[1][0] = 5;
  m[1][2] = m[2][1] = 3;
  m[2][3] = m[3][2] = 3;
  m[0][2] = m[2][0] = 2;
  m[0][3] = m[3][0] = 2;
  m[1][3] = m[3][1] = 2;
  return build_coxeter("H4", 4, m, true, 14400, 60);
}

inline CoxeterGroup make_F4() {
  int m[COX_RMAX][COX_RMAX] = {};
  m[0][1] = m[1][0] = 3;
  m[1][2] = m[2][1] = 4;
  m[2][3] = m[3][2] = 3;
  m[0][2] = m[2][0] = 2;
  m[0][3] = m[3][0] = 2;
  m[1][3] = m[3][1] = 2;
  return build_coxeter("F4", 4, m, false, 1152, 24, 1, 2);
}

inline CoxeterGroup make_B(int n) {
  if (n < 2 || n > 5) throw std::runtime_error("Bn for n=2..5");
  int m[COX_RMAX][COX_RMAX] = {};
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      m[i][j] = (i == j) ? 1 : (std::abs(i - j) == 1 ? 3 : 2);
  m[0][1] = m[1][0] = 4;
  const int order[6] = {0, 0, 8, 48, 384, 3840};
  const int w0[6] = {0, 0, 4, 9, 16, 25};
  return build_coxeter("B" + std::to_string(n), n, m, false, order[n], w0[n], 1, 0);
}
