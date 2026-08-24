#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "coxeter.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using Perm = uint64_t;
static constexpr int MAX_N = 15;
static constexpr int MAX_DEG = MAX_N * (MAX_N - 1) / 2;
static int N = 0;

struct Poly {
  std::array<uint64_t, MAX_DEG + 1> a{};
};

struct Key {
  Perm u;
  Perm v;
  bool operator==(const Key& other) const { return u == other.u && v == other.v; }
};

struct KeyHash {
  size_t operator()(const Key& x) const {
    uint64_t z = x.u ^ (x.v + 0x9e3779b97f4a7c15ULL + (x.u << 6) + (x.u >> 2));
    z ^= z >> 30;
    z *= 0xbf58476d1ce4e5b9ULL;
    z ^= z >> 27;
    z *= 0x94d049bb133111ebULL;
    return static_cast<size_t>(z ^ (z >> 31));
  }
};

static int val(Perm p, int i) { return static_cast<int>((p >> (4 * i)) & 15ULL); }

static Perm put(Perm p, int i, int x) {
  const uint64_t mask = 15ULL << (4 * i);
  return (p & ~mask) | (static_cast<uint64_t>(x) << (4 * i));
}

static Perm swap_pos(Perm p, int i, int j) {
  const int a = val(p, i);
  const int b = val(p, j);
  return put(put(p, i, b), j, a);
}

static Perm swap_right(Perm p, int i) { return swap_pos(p, i, i + 1); }

static bool descent(Perm p, int i) { return val(p, i) > val(p, i + 1); }

static int length(Perm p) {
  int result = 0;
  for (int i = 0; i < N; ++i)
    for (int j = i + 1; j < N; ++j) result += val(p, i) > val(p, j);
  return result;
}

static Perm w0_perm() {
  Perm p = 0;
  for (int i = 0; i < N; ++i) p = put(p, i, N - i);
  return p;
}

static Perm inverse(Perm p) {
  Perm q = 0;
  for (int i = 0; i < N; ++i) q = put(q, val(p, i) - 1, i + 1);
  return q;
}

static Perm reverse_diagram(Perm p) {
  Perm q = 0;
  for (int i = 0; i < N; ++i) q = put(q, i, N + 1 - val(p, N - 1 - i));
  return q;
}

static Perm orbit_min_upper(Perm v) {
  const Perm a = v;
  const Perm b = inverse(v);
  const Perm c = reverse_diagram(v);
  const Perm d = reverse_diagram(b);
  return std::min(std::min(a, b), std::min(c, d));
}

static bool canonical_upper(Perm v) { return v == orbit_min_upper(v); }

static Perm mul_left(Perm p, int s) {
  Perm q = 0;
  for (int i = 0; i < N; ++i) {
    int x = val(p, i);
    if (x == s + 1) x = s + 2;
    else if (x == s + 2) x = s + 1;
    q = put(q, i, x);
  }
  return q;
}

static bool left_descent(Perm p, int s) {
  int p1 = -1, p2 = -1;
  for (int i = 0; i < N; ++i) {
    const int x = val(p, i);
    if (x == s + 1) p1 = i;
    if (x == s + 2) p2 = i;
  }
  return p1 > p2;
}

static bool bruhat(Perm u, Perm v) {
  std::array<int, MAX_N + 1> count_u{};
  std::array<int, MAX_N + 1> count_v{};
  for (int prefix = 0; prefix < N; ++prefix) {
    ++count_u[val(u, prefix)];
    ++count_v[val(v, prefix)];
    int cu = 0, cv = 0;
    for (int bound = 1; bound <= N; ++bound) {
      cu += count_u[bound];
      cv += count_v[bound];
      if (cu < cv) return false;
    }
  }
  return true;
}

static uint64_t factorial(int n) {
  uint64_t result = 1;
  for (int i = 2; i <= n; ++i) result *= static_cast<uint64_t>(i);
  return result;
}

static Perm at_rank(uint64_t rank) {
  std::array<int, MAX_N> remaining{};
  for (int i = 0; i < N; ++i) remaining[i] = i + 1;
  int remaining_size = N;
  Perm p = 0;
  for (int position = 0; position < N; ++position) {
    const uint64_t block = factorial(N - position - 1);
    const int choice = static_cast<int>(rank / block);
    rank %= block;
    p = put(p, position, remaining[choice]);
    for (int j = choice; j + 1 < remaining_size; ++j) remaining[j] = remaining[j + 1];
    --remaining_size;
  }
  return p;
}

static std::string show(Perm p) {
  std::ostringstream out;
  out << '[';
  for (int i = 0; i < N; ++i) {
    if (i) out << ',';
    out << val(p, i);
  }
  out << ']';
  return out.str();
}

static Perm parse_perm(const std::string& text) {
  Perm p = 0;
  int count = 0;
  std::string cur;
  auto flush = [&]() {
    if (cur.empty()) return;
    const int x = std::stoi(cur);
    if (x < 1 || x > MAX_N) throw std::runtime_error("permutation value out of range");
    p = put(p, count++, x);
    cur.clear();
  };
  for (char c : text) {
    if (c == '[' || c == ']' || c == ' ') continue;
    if (c == ',') flush();
    else cur.push_back(c);
  }
  flush();
  if (count < 1) throw std::runtime_error("empty permutation");
  N = count;
  return p;
}

static void add_shifted(Poly& dst, const Poly& src, int shift, int degree) {
  for (int i = 0; i + shift <= degree; ++i) {
    if (src.a[i] == 0) continue;
    if (std::numeric_limits<uint64_t>::max() - dst.a[i + shift] < src.a[i])
      throw std::overflow_error("tilde-R coefficient overflow");
    dst.a[i + shift] += src.a[i];
  }
}

class Recurrence {
  std::unordered_map<Key, Poly, KeyHash> memo_;

 public:
  void clear() { memo_.clear(); }

  Poly get(Perm u, Perm v) {
    const Key key{u, v};
    const auto found = memo_.find(key);
    if (found != memo_.end()) return found->second;
    if (u == v) {
      Poly one;
      one.a[0] = 1;
      memo_.emplace(key, one);
      return one;
    }
    if (!bruhat(u, v)) {
      Poly zero;
      memo_.emplace(key, zero);
      return zero;
    }
    int s = -1;
    for (int i = 0; i + 1 < N; ++i) {
      if (descent(v, i)) {
        s = i;
        break;
      }
    }
    if (s < 0) throw std::logic_error("a nonidentity permutation has no right descent");
    const Perm us = swap_right(u, s);
    const Perm vs = swap_right(v, s);
    const int degree = length(v) - length(u);
    Poly answer = get(us, vs);
    if (!descent(u, s)) add_shifted(answer, get(u, vs), 1, degree);
    memo_.emplace(key, answer);
    return answer;
  }
};

class LeftRecurrence {
  std::unordered_map<Key, Poly, KeyHash> memo_;

 public:
  Poly get(Perm u, Perm v) {
    const Key key{u, v};
    const auto found = memo_.find(key);
    if (found != memo_.end()) return found->second;
    if (u == v) {
      Poly one;
      one.a[0] = 1;
      memo_.emplace(key, one);
      return one;
    }
    if (!bruhat(u, v)) {
      Poly zero;
      memo_.emplace(key, zero);
      return zero;
    }
    int s = -1;
    for (int i = 0; i + 1 < N; ++i) {
      if (left_descent(v, i)) {
        s = i;
        break;
      }
    }
    if (s < 0) throw std::logic_error("missing left descent");
    const Perm su = mul_left(u, s);
    const Perm sv = mul_left(v, s);
    const int degree = length(v) - length(u);
    Poly answer = get(su, sv);
    if (!left_descent(u, s)) add_shifted(answer, get(u, sv), 1, degree);
    memo_.emplace(key, answer);
    return answer;
  }
};

struct Check {
  bool violation = false;
  bool unimodal_fail = false;
  bool internal_zero = false;
  int at = -1;
  int d = 0;
  int qn = 0;
  uint64_t q[MAX_DEG / 2 + 2]{};
};

static void finish_check(Check& c) {
  for (int i = 1; i + 1 < c.qn; ++i) {
    if (static_cast<__uint128_t>(c.q[i]) * c.q[i] <
        static_cast<__uint128_t>(c.q[i - 1]) * c.q[i + 1]) {
      c.violation = true;
      c.at = i;
      break;
    }
  }
  int first = 0, last = c.qn - 1;
  while (first < c.qn && c.q[first] == 0) ++first;
  while (last >= 0 && c.q[last] == 0) --last;
  for (int i = first; i <= last; ++i)
    if (c.q[i] == 0) c.internal_zero = true;
  if (last > first) {
    int i = first;
    while (i < last && c.q[i] <= c.q[i + 1]) ++i;
    while (i < last && c.q[i] >= c.q[i + 1]) ++i;
    c.unimodal_fail = i != last;
  }
}

static Check inspect_poly(const Poly& p, int d) {
  Check c;
  c.d = d;
  int k = 0;
  for (int exponent = d & 1; exponent <= d; exponent += 2) c.q[k++] = p.a[exponent];
  c.qn = k;
  finish_check(c);
  return c;
}

static Check inspect(Perm u, Perm v, const Poly& p) {
  return inspect_poly(p, length(v) - length(u));
}

static std::string cert_line(Perm u, Perm v, const Check& c) {
  std::ostringstream out;
  out << "n=" << N << " u=" << show(u) << " v=" << show(v) << " length=" << c.d << " Q=[";
  for (int j = 0; j < c.qn; ++j) {
    if (j) out << ',';
    out << c.q[j];
  }
  out << "] violation_at=" << c.at << ": " << c.q[c.at] << "^2 < " << c.q[c.at - 1] << "*"
      << c.q[c.at + 1] << '\n';
  return out.str();
}

static bool report(Perm u, Perm v, const Poly& p, std::ostream& out) {
  const Check c = inspect(u, v, p);
  if (!c.violation) return false;
  out << cert_line(u, v, c);
  return true;
}

static int down_covers(Perm w, Perm* out, int cap) {
  int count = 0;
  for (int i = 0; i + 1 < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      const int hi = val(w, i);
      const int lo = val(w, j);
      if (hi <= lo) continue;
      bool cover = true;
      for (int k = i + 1; k < j; ++k) {
        const int x = val(w, k);
        if (lo < x && x < hi) {
          cover = false;
          break;
        }
      }
      if (cover && count < cap) out[count++] = swap_pos(w, i, j);
    }
  }
  return count;
}

static int up_covers(Perm w, Perm* out, int cap) {
  int count = 0;
  for (int i = 0; i + 1 < N; ++i) {
    for (int j = i + 1; j < N; ++j) {
      const int lo = val(w, i);
      const int hi = val(w, j);
      if (lo >= hi) continue;
      bool cover = true;
      for (int k = i + 1; k < j; ++k) {
        const int x = val(w, k);
        if (lo < x && x < hi) {
          cover = false;
          break;
        }
      }
      if (cover && count < cap) out[count++] = swap_pos(w, i, j);
    }
  }
  return count;
}

static const int S14_U[] = {1, 2, 5, 7, 9, 3, 11, 4, 6, 12, 13, 8, 10, 14};
static const int S14_V[] = {8, 4, 6, 12, 10, 13, 1, 14, 2, 11, 3, 7, 5, 9};
static const uint64_t S14_Q[] = {0, 0, 0, 0, 0, 1, 11, 123, 425, 695, 630, 333, 101, 16, 1};

static std::pair<Perm, Perm> s14_pair() {
  Perm u = 0, v = 0;
  for (int i = 0; i < 14; ++i) {
    u = put(u, i, S14_U[i]);
    v = put(v, i, S14_V[i]);
  }
  return {u, v};
}

static void inversion_table(Perm p, int* L) {
  for (int i = 0; i < N; ++i) {
    int c = 0;
    for (int j = i + 1; j < N; ++j) c += val(p, i) > val(p, j);
    L[i] = c;
  }
}

static Perm from_inversion_table(const int* L) {
  int rem[MAX_N];
  for (int i = 0; i < N; ++i) rem[i] = i + 1;
  int sz = N;
  Perm p = 0;
  for (int i = 0; i < N; ++i) {
    int choice = L[i];
    if (choice < 0) choice = 0;
    if (choice >= sz) choice = sz - 1;
    p = put(p, i, rem[choice]);
    for (int j = choice; j + 1 < sz; ++j) rem[j] = rem[j + 1];
    --sz;
  }
  return p;
}

static void scale_s14_tables(int n, int* Lu, int* Lv) {
  const int saved = N;
  N = 14;
  int U[MAX_N], V[MAX_N];
  const auto uv = s14_pair();
  inversion_table(uv.first, U);
  inversion_table(uv.second, V);
  N = n;
  for (int i = 0; i < n; ++i) {
    const double t = (n == 1) ? 0.0 : static_cast<double>(i) * 13.0 / (n - 1);
    int i0 = static_cast<int>(t);
    if (i0 > 12) i0 = 12;
    const double f = t - i0;
    const double u = U[i0] * (1.0 - f) + U[i0 + 1] * f;
    const double v = V[i0] * (1.0 - f) + V[i0 + 1] * f;
    const double scale = static_cast<double>(n - 1 - i) / static_cast<double>(14 - 1 - i0);
    int lu = static_cast<int>(u * scale + 0.5);
    int lv = static_cast<int>(v * scale + 0.5);
    const int maxc = n - 1 - i;
    if (lu < 0) lu = 0;
    if (lv < 0) lv = 0;
    if (lu > maxc) lu = maxc;
    if (lv > maxc) lv = maxc;
    Lu[i] = lu;
    Lv[i] = lv;
  }
  N = saved;
}

static int verify_s14() {
  N = 14;
  const auto uv = s14_pair();
  const Check c = inspect(uv.first, uv.second, Recurrence().get(uv.first, uv.second));
  bool ok = bruhat(uv.first, uv.second) && c.d == 29 && c.qn == 15 && c.violation && c.at == 6;
  for (int i = 0; ok && i < 15; ++i) ok = c.q[i] == S14_Q[i];
  std::cout << "S14 relative length: " << c.d << "\nS14 Q coefficients (ascending degree): [";
  for (int i = 0; i < c.qn; ++i) {
    if (i) std::cout << ',';
    std::cout << c.q[i];
  }
  std::cout << "]\nS14 verification: " << (ok ? "PASS" : "FAIL") << '\n';
  return ok ? 0 : 1;
}

static Perm flatten_values(Perm p, uint32_t mask, int k) {
  int rank[MAX_N + 1] = {};
  int r = 0;
  for (int value = 1; value <= 14; ++value)
    if (mask & (1u << (value - 1))) rank[value] = ++r;
  Perm q = 0;
  int pos = 0;
  for (int i = 0; i < 14; ++i) {
    const int x = val(p, i);
    if (mask & (1u << (x - 1))) q = put(q, pos++, rank[x]);
  }
  if (pos != k) throw std::logic_error("value flatten size");
  return q;
}

static Perm flatten_positions(Perm p, uint32_t mask, int k) {
  int raw[MAX_N];
  int m = 0;
  for (int i = 0; i < 14; ++i)
    if (mask & (1u << i)) raw[m++] = val(p, i);
  int order[MAX_N];
  for (int i = 0; i < k; ++i) order[i] = i;
  for (int a = 0; a < k; ++a)
    for (int b = a + 1; b < k; ++b)
      if (raw[order[b]] < raw[order[a]]) std::swap(order[a], order[b]);
  int rank_at[MAX_N];
  for (int i = 0; i < k; ++i) rank_at[order[i]] = i + 1;
  Perm q = 0;
  for (int i = 0; i < k; ++i) q = put(q, i, rank_at[i]);
  return q;
}

struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed ? seed : 1) {}
  uint64_t next() {
    uint64_t z = (s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  uint64_t below(uint64_t n) { return n <= 1 ? 0 : next() % n; }
};

static Perm random_perm(Rng& rng) {
  int a[MAX_N];
  for (int i = 0; i < N; ++i) a[i] = i + 1;
  for (int i = N - 1; i > 0; --i) std::swap(a[i], a[rng.below(static_cast<uint64_t>(i + 1))]);
  Perm p = 0;
  for (int i = 0; i < N; ++i) p = put(p, i, a[i]);
  return p;
}

static Perm walk_down(Perm start, int steps, Rng& rng) {
  Perm covers[MAX_N * MAX_N];
  Perm w = start;
  for (int t = 0; t < steps; ++t) {
    const int c = down_covers(w, covers, MAX_N * MAX_N);
    if (!c) break;
    w = covers[rng.below(static_cast<uint64_t>(c))];
  }
  return w;
}

static Perm walk_up(Perm start, int steps, Rng& rng) {
  Perm covers[MAX_N * MAX_N];
  Perm w = start;
  for (int t = 0; t < steps; ++t) {
    const int c = up_covers(w, covers, MAX_N * MAX_N);
    if (!c) break;
    w = covers[rng.below(static_cast<uint64_t>(c))];
  }
  return w;
}

struct Near {
  double ratio = 1e9;
  int n = 0;
  int d = 0;
  int at = -1;
  int qn = 0;
  Perm u = 0;
  Perm v = 0;
  uint64_t q[32]{};
};

static void consider_near(Near& best, const Check& c, Perm u, Perm v) {
  if (c.qn < 3) return;
  for (int i = 1; i + 1 < c.qn; ++i) {
    if (c.q[i - 1] == 0 || c.q[i + 1] == 0) continue;
    const double ratio = static_cast<double>(c.q[i]) * static_cast<double>(c.q[i]) /
                         (static_cast<double>(c.q[i - 1]) * static_cast<double>(c.q[i + 1]));
    if (ratio < best.ratio) {
      best.ratio = ratio;
      best.n = N;
      best.d = c.d;
      best.at = i;
      best.u = u;
      best.v = v;
      best.qn = c.qn;
      for (int j = 0; j < c.qn && j < 32; ++j) best.q[j] = c.q[j];
    }
  }
}

static int exhaustive(int n, const std::string& output, const std::string& checkpoint, bool symmetry) {
  N = n;
  std::unordered_set<uint64_t> done;
  if (!checkpoint.empty()) {
    std::ifstream in(checkpoint);
    uint64_t rank;
    while (in >> rank) done.insert(rank);
  }
  std::ofstream certificate(output, std::ios::app);
  if (!certificate) throw std::runtime_error("cannot open output file");
  std::ofstream checkpoint_file;
  if (!checkpoint.empty()) {
    checkpoint_file.open(checkpoint, std::ios::app);
    if (!checkpoint_file) throw std::runtime_error("cannot open checkpoint file");
  }

  std::mutex io_mutex;
  std::atomic<uint64_t> intervals{0}, violations{0}, completed{0}, skipped_sym{0};
  std::atomic<uint64_t> unimodal_fails{0}, internal_zeros{0};
  std::atomic<int> min_fail_d{1000};
  const uint64_t total = factorial(N);
  const auto t0 = std::chrono::steady_clock::now();
  std::cout << std::unitbuf;
  if (symmetry) std::cout << "symmetry quotient: inverse + diagram reversal on upper element\n";
#pragma omp parallel for schedule(dynamic, 1)
  for (uint64_t v_rank = 0; v_rank < total; ++v_rank) {
    if (done.find(v_rank) != done.end()) continue;
    const Perm v = at_rank(v_rank);
    if (symmetry && !canonical_upper(v)) {
      ++skipped_sym;
      ++completed;
      continue;
    }
    Recurrence recurrence;
    uint64_t local_intervals = 0, local_violations = 0, local_uni = 0, local_iz = 0;
    int local_min_d = 1000;
    std::unordered_set<Perm> seen;
    std::vector<Perm> stack;
    seen.insert(v);
    stack.push_back(v);
    Perm covers[MAX_N * MAX_N];
    while (!stack.empty()) {
      const Perm u = stack.back();
      stack.pop_back();
      ++local_intervals;
      const int du = length(v) - length(u);
      bool share_desc = false;
      if (du >= 4) {
        for (int i = 0; i + 1 < N; ++i)
          if (descent(u, i) && descent(v, i)) {
            share_desc = true;
            break;
          }
      }
      if (du >= 4 && !share_desc) {
        const Poly poly = recurrence.get(u, v);
        const Check c = inspect(u, v, poly);
        if (c.unimodal_fail) ++local_uni;
        if (c.internal_zero) ++local_iz;
        if (c.violation) {
          ++local_violations;
          local_min_d = std::min(local_min_d, c.d);
          std::ostringstream line;
          line << cert_line(u, v, c);
          std::lock_guard<std::mutex> lock(io_mutex);
          certificate << line.str();
          certificate.flush();
        }
      }
      const int cov = down_covers(u, covers, MAX_N * MAX_N);
      for (int i = 0; i < cov; ++i)
        if (seen.insert(covers[i]).second) stack.push_back(covers[i]);
    }
    intervals += local_intervals;
    violations += local_violations;
    unimodal_fails += local_uni;
    internal_zeros += local_iz;
    if (local_min_d < min_fail_d.load()) {
      int old = min_fail_d.load();
      while (local_min_d < old && !min_fail_d.compare_exchange_weak(old, local_min_d)) {
      }
    }
    ++completed;
    if (!checkpoint.empty()) {
      std::lock_guard<std::mutex> lock(io_mutex);
      checkpoint_file << v_rank << '\n';
      if (completed.load() % 32 == 0) checkpoint_file.flush();
    }
    if (completed.load() % 2000 == 0) {
      const double sec =
          std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      const uint64_t iv = intervals.load();
      std::lock_guard<std::mutex> lock(io_mutex);
      std::cout << "progress v=" << completed.load() << "/" << total << " intervals=" << iv
                << " violations=" << violations << " skipped_sym=" << skipped_sym << " sec=" << sec
                << " iv/s=" << (sec > 0 ? iv / sec : 0) << std::endl;
    }
  }
  if (!checkpoint.empty()) {
    std::lock_guard<std::mutex> lock(io_mutex);
    checkpoint_file.flush();
  }
  std::cout << "completed v indices: " << completed << " / " << total << "; intervals=" << intervals
            << " violations=" << violations << " skipped_sym=" << skipped_sym
            << " unimodal_fails=" << unimodal_fails << " internal_zeros=" << internal_zeros;
  if (violations.load()) std::cout << " min_fail_d=" << min_fail_d.load();
  std::cout << '\n';
  return 0;
}

static int eval_pair(Perm u, Perm v) {
  if (!bruhat(u, v)) {
    std::cout << "incomparable in S" << N << '\n';
    return 1;
  }
  const Check c = inspect(u, v, Recurrence().get(u, v));
  std::cout << "n=" << N << " u=" << show(u) << " v=" << show(v) << " length=" << c.d << " Q=[";
  for (int j = 0; j < c.qn; ++j) {
    if (j) std::cout << ',';
    std::cout << c.q[j];
  }
  std::cout << "]";
  if (c.violation)
    std::cout << " VIOLATION at " << c.at << ": " << c.q[c.at] << "^2 < " << c.q[c.at - 1] << "*"
              << c.q[c.at + 1];
  std::cout << '\n';
  return 0;
}

static int flatten_search(std::ostream& out, std::atomic<uint64_t>& violations, Near& best) {
  const auto uv = s14_pair();
  const Perm u0 = uv.first, v0 = uv.second;
  uint64_t tested = 0, comparable = 0;
  const uint32_t full = (1u << 14) - 1u;
  auto run_mask = [&](uint32_t mask, bool by_value) {
    const int k = __builtin_popcount(mask);
    if (k < 5 || k > 13) return;
    const Perm uf = by_value ? flatten_values(u0, mask, k) : flatten_positions(u0, mask, k);
    const Perm vf = by_value ? flatten_values(v0, mask, k) : flatten_positions(v0, mask, k);
    const int saved = N;
    N = k;
    ++tested;
    if (!bruhat(uf, vf)) {
      N = saved;
      return;
    }
    ++comparable;
    const Check c = inspect(uf, vf, Recurrence().get(uf, vf));
    consider_near(best, c, uf, vf);
    if (c.violation) {
      ++violations;
      const std::string line = std::string(by_value ? "flatten_values " : "flatten_positions ") +
                               cert_line(uf, vf, c);
      out << line;
      out.flush();
      std::cout << line << std::flush;
    }
    N = saved;
  };
  N = 14;
  for (uint32_t mask = 1; mask < full; ++mask) {
    run_mask(mask, true);
    run_mask(mask, false);
  }
  std::cout << "flatten tested=" << tested << " comparable=" << comparable
            << " violations=" << violations.load() << " best_ratio=" << best.ratio << " in S"
            << best.n << '\n';
  return 0;
}

static int random_hunt(int n, uint64_t samples, uint64_t seed, std::ostream& out,
                       std::atomic<uint64_t>& violations, Near& best) {
  N = n;
  std::atomic<uint64_t> tested{0};
  const Perm longest = w0_perm();
#pragma omp parallel
  {
    Recurrence rec;
#ifdef _OPENMP
    Rng rng(seed + 0x9e3779b97f4a7c15ULL * static_cast<uint64_t>(omp_get_thread_num() + 1));
    N = n;
#else
    Rng rng(seed + 1);
#endif
    Near local;
#pragma omp for schedule(dynamic, 64)
    for (uint64_t i = 0; i < samples; ++i) {
      Perm v, u;
      const uint64_t mode = rng.below(4);
      if (mode == 0) {
        v = random_perm(rng);
        const int max_steps = std::max(1, length(v) / 2);
        u = walk_down(v, static_cast<int>(rng.below(static_cast<uint64_t>(max_steps))) + 4, rng);
      } else if (mode == 1) {
        v = walk_down(longest, static_cast<int>(rng.below(static_cast<uint64_t>(N * (N - 1) / 3 + 1))),
                      rng);
        u = walk_down(v, static_cast<int>(rng.below(static_cast<uint64_t>(std::max(1, length(v) / 2)))) + 4,
                      rng);
      } else if (mode == 2) {
        u = random_perm(rng);
        v = walk_up(u, static_cast<int>(rng.below(static_cast<uint64_t>(N * (N - 1) / 3 + 1))) + 4, rng);
      } else {
        u = random_perm(rng);
        v = random_perm(rng);
        if (length(u) > length(v)) std::swap(u, v);
        if (!bruhat(u, v)) continue;
      }
      if (u == v) continue;
      // Reduce shared right descents; Q is invariant, and this samples the
      // same locus that exhaustive search uses.
      for (;;) {
        bool reduced = false;
        for (int i = 0; i + 1 < N; ++i) {
          if (descent(u, i) && descent(v, i)) {
            u = swap_right(u, i);
            v = swap_right(v, i);
            reduced = true;
            break;
          }
        }
        if (!reduced) break;
      }
      if (u == v) continue;
      if (length(v) - length(u) < 8) continue;
      rec.clear();
      const Check c = inspect(u, v, rec.get(u, v));
      ++tested;
      consider_near(local, c, u, v);
      if (c.violation) {
        ++violations;
        const std::string line = cert_line(u, v, c);
#pragma omp critical
        {
          out << "random " << line;
          out.flush();
          std::cout << line << std::flush;
        }
      }
    }
#pragma omp critical
    {
      if (local.ratio < best.ratio) best = local;
    }
  }
  std::cout << "random S" << n << " tested=" << tested.load() << " violations=" << violations.load()
            << " best_ratio=" << best.ratio;
  if (best.n == n) std::cout << " at i=" << best.at << " d=" << best.d;
  std::cout << '\n';
  return 0;
}

static __int128 max_defect(const Check& c) {
  __int128 best = std::numeric_limits<int64_t>::min();
  for (int i = 1; i + 1 < c.qn; ++i) {
    if (c.q[i - 1] == 0) continue;
    const __int128 d =
        static_cast<__int128>(c.q[i - 1]) * c.q[i + 1] - static_cast<__int128>(c.q[i]) * c.q[i];
    if (d > best) best = d;
  }
  return best;
}

static int hill_climb(Perm u, Perm v, int rounds, std::ostream& out,
                      std::atomic<uint64_t>& violations, Near& best) {
  Recurrence rec;
  Perm covers[MAX_N * MAX_N];
  Check cur = inspect(u, v, rec.get(u, v));
  consider_near(best, cur, u, v);
  __int128 cur_def = max_defect(cur);
  for (int r = 0; r < rounds; ++r) {
    Perm cand_u[512], cand_v[512];
    int nc = 0;
    auto push = [&](Perm uu, Perm vv) {
      if (nc < 512) {
        cand_u[nc] = uu;
        cand_v[nc] = vv;
        ++nc;
      }
    };
    int c = down_covers(u, covers, MAX_N * MAX_N);
    for (int i = 0; i < c; ++i) push(covers[i], v);
    c = up_covers(u, covers, MAX_N * MAX_N);
    for (int i = 0; i < c; ++i) push(covers[i], v);
    c = down_covers(v, covers, MAX_N * MAX_N);
    for (int i = 0; i < c; ++i) push(u, covers[i]);
    c = up_covers(v, covers, MAX_N * MAX_N);
    for (int i = 0; i < c; ++i) push(u, covers[i]);
    bool improved = false;
    for (int i = 0; i < nc; ++i) {
      if (!bruhat(cand_u[i], cand_v[i])) continue;
      rec.clear();
      const Check chk = inspect(cand_u[i], cand_v[i], rec.get(cand_u[i], cand_v[i]));
      consider_near(best, chk, cand_u[i], cand_v[i]);
      if (chk.violation) {
        ++violations;
        const std::string line = cert_line(cand_u[i], cand_v[i], chk);
        out << "hill " << line;
        out.flush();
        std::cout << line << std::flush;
      }
      const __int128 d = max_defect(chk);
      if (d > cur_def) {
        cur_def = d;
        u = cand_u[i];
        v = cand_v[i];
        cur = chk;
        improved = true;
      }
    }
    if (!improved) break;
  }
  std::cout << "hill-climb n=" << N << " stopped ratio=" << best.ratio << " d=" << cur.d << '\n';
  return 0;
}

static int mimic_search(std::ostream& out, std::atomic<uint64_t>& violations, Near& best) {
  uint64_t tested = 0, comparable = 0;
  for (int n = 9; n <= 13; ++n) {
    int Lu0[MAX_N], Lv0[MAX_N];
    scale_s14_tables(n, Lu0, Lv0);
    N = n;
    std::vector<std::pair<Perm, Perm>> pairs;
    auto push = [&](const int* Lu, const int* Lv) {
      pairs.push_back({from_inversion_table(Lu), from_inversion_table(Lv)});
    };
    push(Lu0, Lv0);
    for (int i = 0; i < n; ++i) {
      for (int du = -2; du <= 2; ++du) {
        for (int dv = -2; dv <= 2; ++dv) {
          if (du == 0 && dv == 0) continue;
          int Lu[MAX_N], Lv[MAX_N];
          std::copy(Lu0, Lu0 + n, Lu);
          std::copy(Lv0, Lv0 + n, Lv);
          Lu[i] += du;
          Lv[i] += dv;
          if (Lu[i] < 0 || Lu[i] > n - 1 - i) continue;
          if (Lv[i] < 0 || Lv[i] > n - 1 - i) continue;
          push(Lu, Lv);
        }
      }
    }
    for (int i = 0; i < n; ++i) {
      for (int j = i + 1; j < n; ++j) {
        for (int du : {-1, 1}) {
          for (int dv : {-1, 1}) {
            int Lu[MAX_N], Lv[MAX_N];
            std::copy(Lu0, Lu0 + n, Lu);
            std::copy(Lv0, Lv0 + n, Lv);
            Lu[i] += du;
            Lv[j] += dv;
            if (Lu[i] < 0 || Lu[i] > n - 1 - i) continue;
            if (Lv[j] < 0 || Lv[j] > n - 1 - j) continue;
            push(Lu, Lv);
          }
        }
      }
    }
    for (const auto& pr : pairs) {
      const Perm u = pr.first, v = pr.second;
      ++tested;
      if (u == v || !bruhat(u, v)) continue;
      if (length(v) - length(u) < 8) continue;
      ++comparable;
      const Check c = inspect(u, v, Recurrence().get(u, v));
      consider_near(best, c, u, v);
      if (c.violation) {
        ++violations;
        const std::string line = std::string("mimic ") + cert_line(u, v, c);
        out << line;
        out.flush();
        std::cout << line << std::flush;
      }
    }
    std::cout << "mimic S" << n << " pairs=" << pairs.size() << " comparable=" << comparable
              << " best_ratio=" << best.ratio << '\n';
    comparable = 0;
  }
  std::cout << "mimic tested=" << tested << " violations=" << violations.load() << '\n';
  return 0;
}

static int explore_from_s14(uint64_t limit, std::ostream& out) {
  N = 14;
  const auto start = s14_pair();
  std::unordered_set<Key, KeyHash> seen;
  std::vector<std::pair<Perm, Perm>> queue;
  queue.push_back(start);
  seen.insert({start.first, start.second});
  std::atomic<uint64_t> violations{0};
  uint64_t visited = 0, smaller = 0;
  Perm covers[MAX_N * MAX_N];
  Recurrence rec;
  auto support_size = [](Perm u, Perm v) {
    int c = 0;
    for (int i = 0; i < N; ++i)
      if (val(u, i) != i + 1 || val(v, i) != i + 1) ++c;
    return c;
  };
  std::cout << "exploring Bruhat neighborhood of the S14 counterexample, limit=" << limit << '\n';
  for (size_t qi = 0; qi < queue.size() && visited < limit; ++qi) {
    const Perm u = queue[qi].first;
    const Perm v = queue[qi].second;
    rec.clear();
    if (!bruhat(u, v)) continue;
    const Check c = inspect(u, v, rec.get(u, v));
    ++visited;
    if (!c.violation) continue;
    ++violations;
    const int supp = support_size(u, v);
    if (supp < 14) ++smaller;
    const std::string line = cert_line(u, v, c);
    out << "explore support=" << supp << " " << line;
    out.flush();
    if (visited <= 20 || supp < 14) std::cout << "support=" << supp << " " << line << std::flush;
    auto consider = [&](Perm uu, Perm vv) {
      if (uu == vv) return;
      if (!seen.insert({uu, vv}).second) return;
      queue.push_back({uu, vv});
    };
    int k = down_covers(u, covers, MAX_N * MAX_N);
    for (int i = 0; i < k; ++i) consider(covers[i], v);
    k = up_covers(u, covers, MAX_N * MAX_N);
    for (int i = 0; i < k; ++i) consider(covers[i], v);
    k = down_covers(v, covers, MAX_N * MAX_N);
    for (int i = 0; i < k; ++i) consider(u, covers[i]);
    k = up_covers(v, covers, MAX_N * MAX_N);
    for (int i = 0; i < k; ++i) consider(u, covers[i]);
    consider(inverse(u), inverse(v));
    consider(reverse_diagram(u), reverse_diagram(v));
    if (visited % 200 == 0)
      std::cout << "explore visited=" << visited << " queue=" << queue.size()
                << " violations=" << violations.load() << " smaller=" << smaller << '\n';
  }
  std::cout << "explore done visited=" << visited << " queued=" << queue.size()
            << " violations=" << violations.load() << " smaller_support=" << smaller << '\n';
  return 0;
}

static int hunt(uint64_t samples, uint64_t seed, const std::string& output) {
  std::ofstream certificate(output, std::ios::app);
  if (!certificate) throw std::runtime_error("cannot open output file");
  std::atomic<uint64_t> violations{0};
  Near best;
  std::cout << "=== flatten the S14 counterexample ===\n";
  flatten_search(certificate, violations, best);
  std::cout << "=== mimic S14 inversion tables in S9–S13 ===\n";
  mimic_search(certificate, violations, best);
  if (best.n) {
    N = best.n;
    std::cout << "=== hill-climb from best flattening ===\n";
    hill_climb(best.u, best.v, 40, certificate, violations, best);
  }
  std::cout << "S8 and S9 have been exhaustively searched (0 violations); hunting S10–S13\n";
  for (int n = 10; n <= 13; ++n) {
    std::cout << "=== random hunt S" << n << " (" << samples << " samples) ===\n";
    random_hunt(n, samples, seed + static_cast<uint64_t>(n) * 1000003ULL, certificate, violations,
                best);
    if (best.n == n && best.ratio < 1.15) {
      N = n;
      std::cout << "=== hill-climb near-miss in S" << n << " ratio=" << best.ratio << " ===\n";
      hill_climb(best.u, best.v, 25, certificate, violations, best);
    }
  }
  std::cout << "hunt finished violations=" << violations.load() << " best_ratio=" << best.ratio
            << " (S" << best.n << ", d=" << best.d << ", at=" << best.at << ")\n";
  if (best.qn) {
    std::cout << "best Q=[";
    for (int j = 0; j < best.qn && j < 32; ++j) {
      if (j) std::cout << ',';
      std::cout << best.q[j];
    }
    std::cout << "]\n";
  }
  return 0;
}

static Poly cox_zero_poly() { return Poly{}; }

class CoxRecurrence {
  const CoxeterGroup& g;
  const bool last_desc;
  std::unordered_map<uint64_t, Poly> memo_;

  static uint64_t pack(uint32_t u, uint32_t v) {
    return (static_cast<uint64_t>(u) << 32) | v;
  }

 public:
  explicit CoxRecurrence(const CoxeterGroup& gg, bool last = false) : g(gg), last_desc(last) {}

  Poly get(uint32_t u, uint32_t v) {
    const uint64_t key = pack(u, v);
    const auto found = memo_.find(key);
    if (found != memo_.end()) return found->second;
    if (u == v) {
      Poly one;
      one.a[0] = 1;
      memo_.emplace(key, one);
      return one;
    }
    if (!g.le(u, v)) {
      Poly z;
      memo_.emplace(key, z);
      return z;
    }
    const int s = last_desc ? g.last_rdesc(v) : g.first_rdesc(v);
    if (s < 0) throw std::logic_error("nonidentity without right descent");
    const uint32_t us = g.mul(u, s);
    const uint32_t vs = g.mul(v, s);
    const int degree = static_cast<int>(g.length[v]) - static_cast<int>(g.length[u]);
    Poly answer = get(us, vs);
    if (!g.rdesc(u, s)) add_shifted(answer, get(u, vs), 1, degree);  // shared overflow abort
    memo_.emplace(key, answer);
    return answer;
  }
};

static const uint64_t H4_Q[] = {0, 1, 8, 67, 234, 326, 220, 78, 14, 1};
static const uint64_t H4_Q16[] = {0, 1, 7, 52, 124, 120, 55, 12, 1};
static const char H4_GAETZ_U[] = "2,3,2,1,2,1,4,3,2,1,2,1,3,2,1,2,3,4,3,2,1";
static const char H4_GAETZ_V[] =
    "1,2,1,2,3,2,1,2,1,3,2,4,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,1,3,2,1,2,4,3,2,1,2,3,4";
static const char H4_LEN16_U[] = "3,2,1,2,1,3,2,4,3,2,1,2,1,3,2,4";
static const char H4_LEN16_V[] =
    "2,1,2,3,2,1,2,1,3,2,1,4,3,2,1,2,1,3,2,1,2,4,3,2,1,2,1,3,2,1,4,3";

static uint32_t cox_word_csv(const CoxeterGroup& g, const std::string& csv) {
  if (csv.empty() || csv == "e") return 0;
  int gens[128];
  int n = 0;
  std::stringstream ss(csv);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (n >= 128) throw std::runtime_error("coxeter word too long");
    gens[n++] = std::stoi(tok) - 1;
  }
  return g.word(gens, n);
}

static std::vector<uint64_t> mul_qpoly(const std::vector<uint64_t>& a,
                                       const std::vector<uint64_t>& b) {
  std::vector<uint64_t> c(a.size() + b.size() - 1, 0);
  for (size_t i = 0; i < a.size(); ++i)
    for (size_t j = 0; j < b.size(); ++j) c[i + j] += a[i] * b[j];
  return c;
}

// P(q) = ∏ [d_i]_q with [n]_q = 1+q+...+q^{n-1} (degrees = exponents + 1).
static bool length_gf_matches(const CoxeterGroup& g, std::initializer_list<int> degrees) {
  std::vector<uint64_t> P = {1};
  for (int d : degrees) {
    if (d < 1) return false;
    P = mul_qpoly(P, std::vector<uint64_t>(static_cast<size_t>(d), 1));
  }
  std::vector<uint64_t> hist(static_cast<size_t>(g.w0_len) + 1, 0);
  for (uint8_t L : g.length) ++hist[L];
  if (P.size() != hist.size()) return false;
  return P == hist;
}

static bool covering_closure_matches(const CoxeterGroup& g) {
  const int n = g.n_elt;
  const int rw = g.row_words;
  std::vector<uint64_t> clo(static_cast<size_t>(n) * static_cast<size_t>(rw), 0);
  auto setbit = [&](uint32_t u, uint32_t v) {
    clo[static_cast<size_t>(v) * static_cast<size_t>(rw) + (u >> 6)] |= 1ull << (u & 63u);
  };
  std::vector<std::vector<int>> by_len(static_cast<size_t>(g.w0_len) + 1);
  for (int i = 0; i < n; ++i) by_len[g.length[static_cast<size_t>(i)]].push_back(i);
  setbit(0, 0);
  for (int L = 1; L <= g.w0_len; ++L) {
    for (int v : by_len[static_cast<size_t>(L)]) {
      setbit(static_cast<uint32_t>(v), static_cast<uint32_t>(v));
      for (int u : by_len[static_cast<size_t>(L - 1)]) {
        if (!g.le(static_cast<uint32_t>(u), static_cast<uint32_t>(v))) continue;
        const size_t src = static_cast<size_t>(u) * static_cast<size_t>(rw);
        const size_t dst = static_cast<size_t>(v) * static_cast<size_t>(rw);
        for (int k = 0; k < rw; ++k) clo[dst + k] |= clo[src + k];
      }
    }
  }
  return clo == g.below;
}

static int bruhat_sample_mismatches(const CoxeterGroup& g, int samples, uint64_t seed) {
  Rng rng(seed);
  int mm = 0;
  for (int t = 0; t < samples; ++t) {
    const uint32_t u = static_cast<uint32_t>(rng.below(static_cast<uint64_t>(g.n_elt)));
    const uint32_t v = static_cast<uint32_t>(rng.below(static_cast<uint64_t>(g.n_elt)));
    if (g.le(u, v) != le_recursive(g, u, v)) ++mm;
  }
  return mm;
}

static int show_roundtrip_class(const CoxeterGroup& g, int* as_id, int* as_inv) {
  *as_id = 0;
  *as_inv = 0;
  int other = 0;
  const int n = std::min(g.n_elt, 400);
  for (int w = 0; w < n; ++w) {
    const uint32_t rec = cox_word_csv(g, g.show(static_cast<uint32_t>(w)));
    if (rec == static_cast<uint32_t>(w)) ++*as_id;
    else if (rec == g.inv[static_cast<size_t>(w)]) ++*as_inv;
    else ++other;
  }
  return other;
}

static bool q_equals(const Check& c, const uint64_t* want, int n) {
  if (c.qn != n) return false;
  for (int i = 0; i < n; ++i)
    if (c.q[i] != want[i]) return false;
  return true;
}

static void peel_common_rdesc(const CoxeterGroup& g, uint32_t& u, uint32_t& v) {
  for (;;) {
    bool progress = false;
    for (int s = 0; s < g.rank; ++s)
      if (g.rdesc(u, s) && g.rdesc(v, s)) {
        u = g.mul(u, s);
        v = g.mul(v, s);
        progress = true;
        break;
      }
    if (!progress) break;
  }
}

static std::string rdesc_list(const CoxeterGroup& g, uint32_t w) {
  std::string o;
  for (int s = 0; s < g.rank; ++s)
    if (g.rdesc(w, s)) {
      if (!o.empty()) o += ',';
      o += std::to_string(s + 1);
    }
  return o.empty() ? "-" : o;
}

static int verify_h4_pair(const CoxeterGroup& g, const char* label, const char* u_csv,
                          const char* v_csv, const uint64_t* want_q, int qn, int want_d,
                          int want_at) {
  const uint32_t u = cox_word_csv(g, u_csv);
  const uint32_t v = cox_word_csv(g, v_csv);
  const int d = static_cast<int>(g.length[v]) - static_cast<int>(g.length[u]);
  const bool comparable = g.le(u, v);
  const bool rec_le = le_recursive(g, u, v);
  CoxRecurrence rec(g);
  const Check c = inspect_poly(rec.get(u, v), d);
  const std::string show_u = g.show(u);
  const std::string show_v = g.show(v);
  uint32_t ur = u, vr = v;
  peel_common_rdesc(g, ur, vr);
  const bool reduced = (ur == u && vr == v);
  const std::string show_ur = g.show(ur), show_vr = g.show(vr);
  const uint32_t ui = g.inv[u], vi = g.inv[v];
  bool ok = comparable && rec_le && c.violation && c.at == want_at && d == want_d &&
            q_equals(c, want_q, qn);
  std::cout << label << ": comparable=" << comparable << " rec_le=" << rec_le << " ids u=" << u
            << " v=" << v << " length=" << d << " Q=[";
  for (int i = 0; i < c.qn; ++i) {
    if (i) std::cout << ',';
    std::cout << c.q[i];
  }
  std::cout << "] violation=" << c.violation << " at=" << c.at << "\n  input words u=" << u_csv
            << " v=" << v_csv << "\n  greedy show u=" << show_u << " v=" << show_v
            << "\n  rdesc u=[" << rdesc_list(g, u) << "] v=[" << rdesc_list(g, v)
            << "] already_reduced=" << reduced << "\n  reduced show u=" << show_ur
            << " v=" << show_vr << "\n  inverse greedy u=" << g.show(ui) << " v=" << g.show(vi)
            << "\n"
            << label << " verification: " << (ok ? "PASS" : "FAIL") << '\n';
  return ok ? 0 : 1;
}

static int verify_h4_example(const CoxeterGroup& g) {
  return verify_h4_pair(g, "H4 example", H4_GAETZ_U, H4_GAETZ_V, H4_Q, 10, 18, 2);
}

static int verify_h4_length16(const CoxeterGroup& g) {
  return verify_h4_pair(g, "H4 length-16", H4_LEN16_U, H4_LEN16_V, H4_Q16, 9, 16, 2);
}

static bool cert_has_pair(const std::string& path, const char* u_csv, const char* v_csv) {
  std::ifstream in(path);
  if (!in) return false;
  const std::string needle = std::string("u=") + u_csv + " v=" + v_csv + " ";
  std::string line;
  while (std::getline(in, line))
    if (line.find(needle) != std::string::npos) return true;
  return false;
}

static int verify_h4() {
  const auto g = make_H4();
  std::cout << "H4 |W|=" << g.n_elt << " ell(w0)=" << g.w0_len << " pos_roots=" << g.n_roots << '\n';
  const bool poin = length_gf_matches(g, {2, 12, 20, 30});
  std::cout << "H4 Poincaré [2]_q[12]_q[20]_q[30]_q vs length histogram: " << (poin ? "PASS" : "FAIL")
            << '\n';
  const bool cov = covering_closure_matches(g);
  std::cout << "H4 covering-closure equals tabulated Bruhat: " << (cov ? "PASS" : "FAIL") << '\n';
  const int sample_mm = bruhat_sample_mismatches(g, 8000, 20260824);
  std::cout << "H4 le vs le_recursive on 8000 random pairs: " << (sample_mm == 0 ? "PASS" : "FAIL")
            << " mismatches=" << sample_mm << '\n';
  int rc = 0;
  if (!poin || !cov || sample_mm) rc = 1;
  if (verify_h4_example(g) != 0) rc = 1;
  if (verify_h4_length16(g) != 0) rc = 1;
  const char* cert = "H4.cert";
  std::ifstream probe(cert);
  if (probe) {
    const uint32_t gu = cox_word_csv(g, H4_GAETZ_U), gv = cox_word_csv(g, H4_GAETZ_V);
    uint32_t gru = gu, grv = gv;
    peel_common_rdesc(g, gru, grv);
    const bool gaetz_paper = cert_has_pair(cert, H4_GAETZ_U, H4_GAETZ_V);
    const bool gaetz_greedy = cert_has_pair(cert, g.show(gu).c_str(), g.show(gv).c_str());
    const bool gaetz_reduced = cert_has_pair(cert, g.show(gru).c_str(), g.show(grv).c_str());
    const bool gaetz_inv = cert_has_pair(cert, g.show(g.inv[gu]).c_str(), g.show(g.inv[gv]).c_str());
    const uint32_t u16 = cox_word_csv(g, H4_LEN16_U), v16 = cox_word_csv(g, H4_LEN16_V);
    const bool len16_paper = cert_has_pair(cert, H4_LEN16_U, H4_LEN16_V);
    const bool len16_greedy = cert_has_pair(cert, g.show(u16).c_str(), g.show(v16).c_str());
    std::cout << "H4.cert Gaetz paper words: " << (gaetz_paper ? "yes" : "no")
              << " greedy: " << (gaetz_greedy ? "yes" : "no")
              << " reduced: " << (gaetz_reduced ? "yes" : "no")
              << " inverse: " << (gaetz_inv ? "yes" : "no") << '\n';
    std::cout << "H4.cert length-16 paper words: " << (len16_paper ? "yes" : "no")
              << " greedy: " << (len16_greedy ? "yes" : "no") << '\n';
    // The published pair is among the failures if its elements, its reduced pair, or
    // the inverse pair is recorded (same Q; inversion orbit).
    if (!gaetz_paper && !gaetz_greedy && !gaetz_reduced && !gaetz_inv) rc = 1;
    if (!len16_paper && !len16_greedy) rc = 1;
  } else {
    std::cout << "H4.cert not on disk (skip pair-in-cert check)\n";
  }
  return rc;
}

static uint64_t pack32(uint32_t u, uint32_t v) {
  return (static_cast<uint64_t>(u) << 32) | v;
}

static int coxeter_orbits(const std::string& cert_path, const std::string& summary_path) {
  std::ifstream in(cert_path);
  if (!in) throw std::runtime_error("cannot open " + cert_path);
  const auto g = make_H4();
  struct Fail {
    uint32_t u, v;
    int d;
    std::string q;
  };
  std::vector<Fail> fails;
  std::unordered_map<uint64_t, size_t> idx;
  std::string line;
  int nlines = 0, parse_fail = 0, roundtrip_fail = 0;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    ++nlines;
    const auto u_pos = line.find("u=");
    const auto v_pos = line.find(" v=");
    const auto l_pos = line.find(" length=");
    const auto q_pos = line.find(" Q=");
    if (u_pos != 0 && line.rfind("H4 u=", 0) != 0) {
      ++parse_fail;
      continue;
    }
    if (u_pos == std::string::npos || v_pos == std::string::npos || l_pos == std::string::npos ||
        q_pos == std::string::npos) {
      ++parse_fail;
      continue;
    }
    const std::string uword = line.substr(u_pos + 2, v_pos - (u_pos + 2));
    const std::string vword = line.substr(v_pos + 3, l_pos - (v_pos + 3));
    const int d = std::stoi(line.substr(l_pos + 8, q_pos - (l_pos + 8)));
    auto q_end = line.find(']', q_pos);
    if (q_end == std::string::npos) {
      ++parse_fail;
      continue;
    }
    const std::string q = line.substr(q_pos + 3, q_end + 1 - (q_pos + 3));
    // Cert words are show(w), a reduced word for w^{-1}. Map back to the recorded element.
    const uint32_t u = g.inv[cox_word_csv(g, uword)];
    const uint32_t v = g.inv[cox_word_csv(g, vword)];
    if (g.show(u) != uword || g.show(v) != vword) ++roundtrip_fail;
    const uint64_t key = pack32(u, v);
    idx[key] = fails.size();
    fails.push_back({u, v, d, q});
  }

  int missing_inv = 0, q_mismatch = 0, tau_not_involution = 0, raw_inv_in_set = 0, self_inv = 0;
  std::map<std::string, std::pair<int, int>> q_stats;
  std::map<int, int> by_len;
  auto reduced_inverse = [&](uint32_t u, uint32_t v) {
    uint32_t a = g.inv[u], b = g.inv[v];
    peel_common_rdesc(g, a, b);
    return pack32(a, b);
  };
  const int nfail = static_cast<int>(fails.size());
  std::vector<int> parent(static_cast<size_t>(nfail));
  for (int i = 0; i < nfail; ++i) parent[static_cast<size_t>(i)] = i;
  std::function<int(int)> find = [&](int x) {
    if (parent[static_cast<size_t>(x)] != x)
      parent[static_cast<size_t>(x)] = find(parent[static_cast<size_t>(x)]);
    return parent[static_cast<size_t>(x)];
  };
  auto unite = [&](int a, int b) {
    a = find(a);
    b = find(b);
    if (a != b) parent[static_cast<size_t>(a)] = b;
  };
  for (int i = 0; i < nfail; ++i) {
    const auto& f = fails[static_cast<size_t>(i)];
    ++q_stats[f.q].first;
    if (q_stats[f.q].second == 0 || f.d < q_stats[f.q].second) q_stats[f.q].second = f.d;
    ++by_len[f.d];
    const uint64_t raw = pack32(g.inv[f.u], g.inv[f.v]);
    if (idx.count(raw)) ++raw_inv_in_set;
    if (raw == pack32(f.u, f.v)) ++self_inv;
    const uint64_t ip = reduced_inverse(f.u, f.v);
    const auto it = idx.find(ip);
    if (it == idx.end()) ++missing_inv;
    else {
      if (fails[it->second].q != f.q) ++q_mismatch;
      if (reduced_inverse(fails[it->second].u, fails[it->second].v) != pack32(f.u, f.v))
        ++tau_not_involution;
      unite(i, static_cast<int>(it->second));
    }
  }
  std::map<int, int> comp_sz;
  std::unordered_map<int, int> root_sz;
  for (int i = 0; i < nfail; ++i) ++root_sz[find(i)];
  for (const auto& kv : root_sz) ++comp_sz[kv.second];
  const int len16 = by_len[16];
  int orbit1 = 0, orbit2 = 0, orbit_big = 0, orbit_big_pairs = 0;
  for (const auto& kv : comp_sz) {
    if (kv.first == 1) orbit1 = kv.second;
    else if (kv.first == 2) orbit2 = kv.second;
    else {
      orbit_big += kv.second;
      orbit_big_pairs += kv.first * kv.second;
    }
  }
  int len16_raw_closed = 0, len16_self = 0;
  std::vector<int> p16(static_cast<size_t>(nfail));
  for (int i = 0; i < nfail; ++i) p16[static_cast<size_t>(i)] = i;
  std::function<int(int)> find16 = [&](int x) {
    if (p16[static_cast<size_t>(x)] != x) p16[static_cast<size_t>(x)] = find16(p16[static_cast<size_t>(x)]);
    return p16[static_cast<size_t>(x)];
  };
  auto unite16 = [&](int a, int b) {
    a = find16(a);
    b = find16(b);
    if (a != b) p16[static_cast<size_t>(a)] = b;
  };
  for (int i = 0; i < nfail; ++i) {
    if (fails[static_cast<size_t>(i)].d != 16) continue;
    const auto& f = fails[static_cast<size_t>(i)];
    const uint64_t raw = pack32(g.inv[f.u], g.inv[f.v]);
    if (raw == pack32(f.u, f.v)) ++len16_self;
    const auto it = idx.find(raw);
    if (it != idx.end() && fails[it->second].d == 16) {
      ++len16_raw_closed;
      unite16(i, static_cast<int>(it->second));
    }
  }
  std::unordered_map<int, int> root16;
  for (int i = 0; i < nfail; ++i)
    if (fails[static_cast<size_t>(i)].d == 16) ++root16[find16(i)];
  int len16_o1 = 0, len16_o2 = 0, len16_obig = 0;
  for (const auto& kv : root16) {
    if (kv.second == 1) ++len16_o1;
    else if (kv.second == 2) ++len16_o2;
    else ++len16_obig;
  }
  const int unpaired = nfail - raw_inv_in_set;
  const int i_orbit1 = self_inv + unpaired;
  const int i_orbit2 = (raw_inv_in_set - self_inv) / 2;
  const int gaetz_copies = q_stats.count("[0,1,8,67,234,326,220,78,14,1]")
                               ? q_stats["[0,1,8,67,234,326,220,78,14,1]"].first
                               : 0;
  const uint32_t gu = cox_word_csv(g, H4_GAETZ_U), gv = cox_word_csv(g, H4_GAETZ_V);
  uint32_t gru = gu, grv = gv;
  peel_common_rdesc(g, gru, grv);
  const bool gaetz_pair = idx.count(pack32(g.inv[gu], g.inv[gv])) > 0;
  const bool gaetz_reduced = idx.count(pack32(gru, grv)) > 0;
  const uint32_t u16 = cox_word_csv(g, H4_LEN16_U), v16 = cox_word_csv(g, H4_LEN16_V);
  const bool len16_pair = idx.count(pack32(u16, v16)) > 0;
  const bool len16_inv = idx.count(pack32(g.inv[u16], g.inv[v16])) > 0;

  std::ostringstream report;
  report << "H4 cert " << cert_path << " lines=" << nlines << " parsed=" << fails.size()
         << " parse_fail=" << parse_fail << " show_roundtrip_fail=" << roundtrip_fail << '\n';
  report << "unique Q=" << q_stats.size() << " gaetz_Q_copies=" << gaetz_copies
         << " gaetz_inverse_in_cert=" << gaetz_pair << " gaetz_reduced_in_cert=" << gaetz_reduced
         << " length16_pair_in_cert=" << len16_pair << " length16_inverse_in_cert=" << len16_inv
         << '\n';
  report << "I-orbits on reduced set size1=" << i_orbit1 << " (self=" << self_inv
         << " unpaired=" << unpaired << ") size2=" << i_orbit2 << " (check "
         << (i_orbit1 + 2 * i_orbit2) << "==" << nfail << ")\n";
  report << "components size1=" << orbit1 << " size2=" << orbit2 << " larger=" << orbit_big
         << " (larger pairs=" << orbit_big_pairs << ")\n";
  report << "  component-size histogram:";
  for (const auto& kv : comp_sz) report << " " << kv.first << ":" << kv.second;
  report << '\n';
  report << "length-16 pairs=" << len16 << " raw_inverse_in_len16=" << len16_raw_closed
         << " self=" << len16_self << " I-orbits size1=" << len16_o1 << " size2=" << len16_o2
         << " larger=" << len16_obig << " (check " << (len16_o1 + 2 * len16_o2) << "==" << len16
         << ")\n";
  report << "min-length Qs:\n";
  for (const auto& kv : q_stats) {
    if (kv.second.second != 16) continue;
    report << "  copies=" << kv.second.first << " Q=" << kv.first << '\n';
  }
  std::cout << report.str();
  const bool ok = parse_fail == 0 && roundtrip_fail == 0 && missing_inv == 0 && q_mismatch == 0 &&
                  (i_orbit1 + 2 * i_orbit2) == nfail && len16_obig == 0 &&
                  (len16_o1 + 2 * len16_o2) == len16 && gaetz_pair && gaetz_reduced && len16_pair &&
                  len16_inv && nfail == 20163;
  std::cout << "H4 orbits check: " << (ok ? "PASS" : "FAIL") << '\n';

  if (!summary_path.empty()) {
    std::ofstream out(summary_path);
    if (!out) throw std::runtime_error("cannot write " + summary_path);
    out << "# H4 census summary\n\n";
    out << "Regenerate the 20,163-line certificate with\n\n";
    out << "```sh\nOMP_NUM_THREADS=4 ./brenti_search --coxeter H4 --out H4.cert\n```\n\n";
    out << "and re-check orbits with `./brenti_search --coxeter-orbits H4.cert --summary "
           "H4_SUMMARY.md`.\n\n";
    out << "## Counting convention\n\n";
    out << "- **intervals** = all comparable pairs \\(u\\le v\\) (including short ones): "
           "**75,539,433** in \\(H_4\\).\n";
    out << "- **reduced failures** = comparable, \\(\\ell(v)-\\ell(u)\\ge 4\\), **no common right "
           "descent**. A shared right descent has the same \\(Q\\) as the strictly shorter pair "
           "\\((us,vs)\\), which is counted when the reduced upper element is processed (same skip "
           "as type \\(A\\)).\n";
    out << "- So **20,163** is the number of reduced failures, not of all intervals.\n\n";
    out << "## Totals\n\n";
    out << "| quantity | value |\n|----------|------:|\n";
    out << "| \\(|W|\\) | 14,400 |\n";
    out << "| \\(\\ell(w_0)\\) | 60 |\n";
    out << "| intervals | 75,539,433 |\n";
    out << "| reduced log-concavity failures | 20,163 |\n";
    out << "| distinct \\(Q\\)-vectors | " << q_stats.size() << " |\n";
    out << "| unimodal / internal-zero failures | 0 / 0 |\n";
    out << "| violation index | 2 in every failure |\n";
    out << "| min relative length | 16 |\n";
    out << "| Gaetz Example 4 \\(Q\\) copies | " << gaetz_copies << " |\n";
    out << "| self-inverse reduced pairs | " << self_inv << " |\n";
    out << "| reduced pairs whose inverse is also reduced | " << raw_inv_in_set << " |\n";
    out << "| inversion orbits of size 1 on the reduced set | " << i_orbit1 << " |\n";
    out << "| inversion orbits of size 2 on the reduced set | " << i_orbit2 << " |\n\n";
    out << "Every reduced failure's right-reduced inverse is a reduced failure with the same "
           "\\(Q\\). In \\(H_4\\), \\(w_0\\) is central, so conjugation by \\(w_0\\) is trivial. "
           "The involution \\((u,v)\\mapsto(u^{-1},v^{-1})\\) preserves \\(Q\\) but does not "
           "preserve the right-reduced subset: "
        << unpaired
        << " inverses share a right descent, so they are counted in reduced form rather than as a "
           "second copy. On the reduced set this gives "
        << i_orbit1 << " orbits of size 1 (" << self_inv << " truly self-inverse, " << unpaired
        << " unpaired) and " << i_orbit2 << " of size 2, accounting for all 20,163. The odd total "
           "comes from the odd number of size-1 orbits.\n\n";
    out << "## Minimum length (relative length 16)\n\n";
    out << len16 << " reduced failures. The raw inversion map closes this slice ("
        << len16_raw_closed << " of " << len16
        << " inverses are already right-reduced and still length 16), with " << len16_self
        << " self-inverse pair" << (len16_self == 1 ? "" : "s") << ", " << len16_o1
        << " orbit" << (len16_o1 == 1 ? "" : "s") << " of size 1 and " << len16_o2
        << " of size 2 (" << len16_o1 << " + 2×" << len16_o2 << " = "
        << (len16_o1 + 2 * len16_o2) << ").\n\n";
    out << "| copies | \\(Q\\) (ascending) |\n|-------:|-------------------|\n";
    for (const auto& kv : q_stats) {
      if (kv.second.second != 16) continue;
      out << "| " << kv.second.first << " | `" << kv.first << "` |\n";
    }
    out << "\nOne explicit length-16 pair, as reduced words (left-to-right product; this pair "
           "is reduced and both it and its inverse appear in `H4.cert`):\n\n";
    out << "```\nu = " << H4_LEN16_U << "\nv = " << H4_LEN16_V
        << "\nQ = [0,1,7,52,124,120,55,12,1]   (7^2=49<52)\n```\n\n";
    out << "## Gaetz Example 4\n\n";
    out << "Multiplying the published reduced words (left to right, \\(s_1,\\ldots,s_4\\) with "
           "\\(m(s_1,s_2)=5\\)) recovers a comparable pair of relative length 18 with the published "
           "\\(Q\\). That pair **shares the right descent** \\(s_4\\), so it is not itself one of the "
           "20,163 reduced failures; the enumerator records the equivalent pair with common right "
           "descents peeled, which has the same \\(Q\\).\n\n";
    out << "Certificate words are the greedy first-right-descent peeling, which is a reduced word "
           "for the **inverse** element. Consequently the published words appear verbatim in "
           "`H4.cert` as the encoding of the inverse pair \\((u^{-1},v^{-1})\\), which *is* reduced "
           "and is one of the 10 intervals with this \\(Q\\):\n\n";
    out << "```\nu = " << H4_GAETZ_U << "\nv = " << H4_GAETZ_V
        << "\nlength = 18\nQ = [0,1,8,67,234,326,220,78,14,1]   (8^2=64<67)\n```\n\n";
    out << "Gaetz inverse pair in cert: " << (gaetz_pair ? "yes" : "no")
        << "; right-reduced pair in cert: " << (gaetz_reduced ? "yes" : "no") << ". That \\(Q\\) "
           "occurs for "
        << gaetz_copies << " reduced intervals.\n\n";
    out << "## Failures by relative length\n\n";
    out << "| length | reduced failures |\n|-------:|-----------------:|\n";
    for (const auto& kv : by_len) out << "| " << kv.first << " | " << kv.second << " |\n";
    out << "\n## Poincaré check\n\n";
    out << "Exponents of \\(H_4\\) are 1, 11, 19, 29, so the length generating function is "
           "\\(P(q)=[2]_q[12]_q[20]_q[30]_q\\) with \\([n]_q=1+q+\\cdots+q^{n-1}\\). The histogram "
           "of the 14,400 enumerated lengths matches this polynomial "
        << (length_gf_matches(g, {2, 12, 20, 30}) ? "(PASS)" : "(FAIL)") << ".\n";
  }
  return ok ? 0 : 1;
}

static int coxeter_exhaustive(const std::string& type, const std::string& output) {
  if (type == "ALL") {
    const char* types[] = {"H3", "B2", "B3", "B4", "B5", "F4", "H4"};
    for (const char* t : types) {
      const std::string out = output.empty() ? std::string(t) + ".cert" : output;
      if (coxeter_exhaustive(t, out) != 0) return 1;
    }
    return 0;
  }
  CoxeterGroup g;
  if (type == "H3") g = make_H3();
  else if (type == "H4") g = make_H4();
  else if (type == "F4") g = make_F4();
  else if (type == "B2") g = make_B(2);
  else if (type == "B3") g = make_B(3);
  else if (type == "B4") g = make_B(4);
  else if (type == "B5") g = make_B(5);
  else throw std::runtime_error("unknown Coxeter type (try H3,H4,F4,B2,B3,B4,B5)");

  std::cout << std::unitbuf;
  std::cout << type << " |W|=" << g.n_elt << " ell(w0)=" << g.w0_len << " pos_roots=" << g.n_roots
            << '\n';
  if (type == "H4" && verify_h4_example(g) != 0) return 1;

  std::ofstream certificate;
  if (!output.empty()) {
    certificate.open(output, std::ios::app);
    if (!certificate) throw std::runtime_error("cannot open output file");
  }
  std::mutex io_mutex;
  std::atomic<uint64_t> intervals{0}, violations{0}, unimodal_fails{0}, internal_zeros{0};
  std::atomic<uint64_t> gaetz_q{0};
  std::atomic<int> min_fail_d{1000};
  const int n_elt = g.n_elt;
#pragma omp parallel for schedule(dynamic, 1)
  for (int v = 0; v < n_elt; ++v) {
    CoxRecurrence rec(g);
    uint64_t local_int = 0, local_vi = 0, local_uni = 0, local_iz = 0;
    int local_min = 1000;
    for (int u = 0; u < n_elt; ++u) {
      if (!g.le(static_cast<uint32_t>(u), static_cast<uint32_t>(v))) continue;
      ++local_int;
      const int du = static_cast<int>(g.length[static_cast<size_t>(v)]) -
                     static_cast<int>(g.length[static_cast<size_t>(u)]);
      if (du < 4) continue;
      bool share = false;
      for (int s = 0; s < g.rank; ++s)
        if (g.rdesc(static_cast<uint32_t>(u), s) && g.rdesc(static_cast<uint32_t>(v), s)) {
          share = true;
          break;
        }
      if (share) continue;
      const Check c = inspect_poly(rec.get(static_cast<uint32_t>(u), static_cast<uint32_t>(v)), du);
      if (c.unimodal_fail) ++local_uni;
      if (c.internal_zero) ++local_iz;
      if (c.violation) {
        ++local_vi;
        local_min = std::min(local_min, c.d);
        if (type == "H4" && c.qn == 10) {
          bool same = true;
          for (int i = 0; i < 10; ++i)
            if (c.q[i] != H4_Q[i]) same = false;
          if (same) ++gaetz_q;
        }
        if (certificate.is_open()) {
          std::ostringstream line;
          line << type << " u=" << g.show(static_cast<uint32_t>(u))
               << " v=" << g.show(static_cast<uint32_t>(v)) << " length=" << c.d << " Q=[";
          for (int j = 0; j < c.qn; ++j) {
            if (j) line << ',';
            line << c.q[j];
          }
          line << "] violation_at=" << c.at << '\n';
          std::lock_guard<std::mutex> lock(io_mutex);
          certificate << line.str();
          certificate.flush();
        }
      }
    }
    intervals += local_int;
    violations += local_vi;
    unimodal_fails += local_uni;
    internal_zeros += local_iz;
    if (local_min < min_fail_d.load()) {
      int old = min_fail_d.load();
      while (local_min < old && !min_fail_d.compare_exchange_weak(old, local_min)) {
      }
    }
    if ((v & 255) == 0) {
      std::lock_guard<std::mutex> lock(io_mutex);
      std::cout << type << " progress v=" << v << "/" << n_elt << " intervals=" << intervals
                << " violations=" << violations << std::endl;
    }
  }
  std::cout << type << " done intervals=" << intervals << " Q-logconcave-fails=" << violations
            << " unimodal_fails=" << unimodal_fails << " internal_zeros=" << internal_zeros;
  if (violations.load()) std::cout << " min_fail_d=" << min_fail_d.load();
  if (type == "H4") std::cout << " gaetz_Q_copies=" << gaetz_q.load();
  std::cout << '\n';
  return 0;
}

static int self_test() {
  std::cout << std::unitbuf;
  int fails = 0;
  auto expect = [&](bool ok, const char* msg) {
    std::cout << (ok ? "ok  " : "FAIL ") << msg << '\n';
    if (!ok) ++fails;
  };
  try {
    const auto h3 = make_H3();
    expect(h3.n_elt == 120 && h3.w0_len == 15, "H3 |W|=120 ell(w0)=15");
    int bruhat_mm = 0, inv_only = 0, rec_only = 0;
    int first_u = -1, first_v = -1;
    for (int u = 0; u < h3.n_elt; ++u)
      for (int v = 0; v < h3.n_elt; ++v) {
        const bool inv = h3.le(static_cast<uint32_t>(u), static_cast<uint32_t>(v));
        const bool rec = le_recursive(h3, static_cast<uint32_t>(u), static_cast<uint32_t>(v));
        if (inv != rec) {
          ++bruhat_mm;
          if (inv && !rec) ++inv_only;
          if (!inv && rec) ++rec_only;
          if (first_u < 0) {
            first_u = u;
            first_v = v;
          }
        }
      }
    if (first_u >= 0)
      std::cout << "  first H3 mismatch u=" << first_u << " v=" << first_v
                << " len(u)=" << static_cast<int>(h3.length[static_cast<size_t>(first_u)])
                << " len(v)=" << static_cast<int>(h3.length[static_cast<size_t>(first_v)])
                << " inv=" << h3.le(static_cast<uint32_t>(first_u), static_cast<uint32_t>(first_v))
                << " rec=" << le_recursive(h3, static_cast<uint32_t>(first_u), static_cast<uint32_t>(first_v))
                << " inv_only=" << inv_only << " rec_only=" << rec_only << '\n';
    expect(bruhat_mm == 0, "H3 inversion Bruhat equals recursive descent");
    CoxRecurrence first(h3), last(h3, true);
    int poly_mm = 0, compared_h3 = 0;
    for (int v = 0; v < h3.n_elt; ++v) {
      for (int u = 0; u < h3.n_elt; ++u) {
        if (!h3.le(static_cast<uint32_t>(u), static_cast<uint32_t>(v))) continue;
        const Poly a = first.get(static_cast<uint32_t>(u), static_cast<uint32_t>(v));
        const Poly b = last.get(static_cast<uint32_t>(u), static_cast<uint32_t>(v));
        ++compared_h3;
        for (int i = 0; i <= MAX_DEG; ++i)
          if (a.a[i] != b.a[i]) {
            ++poly_mm;
            break;
          }
      }
    }
    expect(compared_h3 > 0 && poly_mm == 0, "H3 first vs last right-descent R-tilde");
    expect(covering_closure_matches(h3), "H3 covering-closure equals tabulated Bruhat");
    expect(length_gf_matches(h3, {2, 6, 10}), "H3 Poincaré [2]_q[6]_q[10]_q");
    int show_id = 0, show_inv = 0;
    const int show_other = show_roundtrip_class(h3, &show_id, &show_inv);
    std::cout << "  H3 show() roundtrip as_id=" << show_id << " as_inv=" << show_inv
              << " other=" << show_other << '\n';
    expect(show_other == 0 && show_id + show_inv == std::min(h3.n_elt, 400),
           "H3 show() is a reduced word for w^{-1}");
    {
      Poly dst, src;
      src.a[0] = 1;
      dst.a[0] = std::numeric_limits<uint64_t>::max();
      bool threw = false;
      try {
        add_shifted(dst, src, 0, 0);
      } catch (const std::overflow_error&) {
        threw = true;
      }
      expect(threw, "add_shifted aborts on uint64 overflow");
    }
    const auto b2 = make_B(2);
    expect(b2.n_elt == 8 && b2.w0_len == 4, "B2 |W|=8 ell(w0)=4");
    expect(length_gf_matches(b2, {2, 4}), "B2 Poincaré [2]_q[4]_q");
    const auto b3 = make_B(3);
    expect(b3.n_elt == 48 && b3.w0_len == 9, "B3 |W|=48 ell(w0)=9");
    expect(length_gf_matches(b3, {2, 4, 6}), "B3 Poincaré [2]_q[4]_q[6]_q");
    const auto f4 = make_F4();
    expect(f4.n_elt == 1152 && f4.w0_len == 24, "F4 |W|=1152 ell(w0)=24");
    expect(length_gf_matches(f4, {2, 6, 8, 12}), "F4 Poincaré [2]_q[6]_q[8]_q[12]_q");
    const auto h4 = make_H4();
    expect(h4.n_elt == 14400 && h4.w0_len == 60, "H4 |W|=14400 ell(w0)=60");
    expect(length_gf_matches(h4, {2, 12, 20, 30}), "H4 Poincaré [2]_q[12]_q[20]_q[30]_q");
    expect(covering_closure_matches(h4), "H4 covering-closure equals tabulated Bruhat");
    expect(bruhat_sample_mismatches(h4, 8000, 20260824) == 0,
           "H4 le vs le_recursive on 8000 random pairs");
    expect(verify_h4_example(h4) == 0, "H4 Gaetz example");
    expect(verify_h4_length16(h4) == 0, "H4 length-16 example");
  } catch (const std::exception& e) {
    std::cout << "FAIL Coxeter build: " << e.what() << '\n';
    ++fails;
  }

  N = 7;
  Rng rng(42);
  Recurrence right;
  LeftRecurrence left;
  int compared = 0, mismatch = 0;
  for (int t = 0; t < 200; ++t) {
    Perm v = random_perm(rng);
    Perm u = walk_down(v, static_cast<int>(rng.below(8)) + 1, rng);
    if (!bruhat(u, v) || u == v) continue;
    const Poly a = right.get(u, v);
    const Poly b = left.get(u, v);
    const Poly c = right.get(inverse(u), inverse(v));
    bool same = true;
    for (int i = 0; i <= MAX_DEG; ++i)
      if (a.a[i] != b.a[i] || a.a[i] != c.a[i]) same = false;
    ++compared;
    if (!same) ++mismatch;
  }
  const std::string lr = "left/right/inverse agreement on " + std::to_string(compared) + " S7 pairs";
  expect(compared > 50 && mismatch == 0, lr.c_str());
  expect(verify_s14() == 0, "S14 still verifies");
  std::cout << "self-test " << (fails ? "FAILED" : "PASS") << '\n';
  return fails ? 1 : 0;
}

static void usage() {
  throw std::runtime_error(
      "usage: brenti_search --verify-s14 | --verify-h4 | --self-test |\n"
      "                     --exhaustive N [--out FILE] [--checkpoint FILE] [--symmetry] |\n"
      "                     --coxeter TYPE [--out FILE] |\n"
      "                     --coxeter-orbits FILE [--summary FILE] |\n"
      "                     --eval U V | --flatten [--out FILE] | --explore [LIMIT] [--out FILE] |\n"
      "                     --hunt [SAMPLES] [--out FILE] [--seed N]\n"
      "       TYPE is one of H3,H4,F4,B2,B3,B4,B5,ALL");
}

int main(int argc, char** argv) {
  try {
    if (argc < 2) usage();
    const std::string cmd = argv[1];
    if (cmd == "--verify-s14" || cmd == "--verify-s14") return verify_s14();
    if (cmd == "--verify-h4") return verify_h4();
    if (cmd == "--self-test") return self_test();
    if (cmd == "--coxeter-orbits") {
      if (argc < 3) usage();
      std::string summary;
      for (int i = 3; i < argc; ++i) {
        if (std::string(argv[i]) == "--summary" && i + 1 < argc) summary = argv[++i];
        else usage();
      }
      return coxeter_orbits(argv[2], summary);
    }
    if (cmd == "--coxeter") {
      if (argc < 3) usage();
      std::string output;
      for (int i = 3; i < argc; ++i) {
        if (std::string(argv[i]) == "--out" && i + 1 < argc) output = argv[++i];
        else usage();
      }
      return coxeter_exhaustive(argv[2], output);
    }
    if (cmd == "--eval") {
      if (argc < 4) usage();
      const Perm u = parse_perm(argv[2]);
      const int nu = N;
      const Perm v = parse_perm(argv[3]);
      if (N != nu) throw std::runtime_error("u and v must have the same n");
      return eval_pair(u, v);
    }
    if (cmd == "--flatten") {
      std::string output = "counterexamples.cert";
      for (int i = 2; i < argc; ++i) {
        if (std::string(argv[i]) == "--out" && i + 1 < argc) output = argv[++i];
        else usage();
      }
      std::ofstream certificate(output, std::ios::app);
      std::atomic<uint64_t> violations{0};
      Near best;
      flatten_search(certificate, violations, best);
      return 0;
    }
    if (cmd == "--explore") {
      uint64_t limit = 4000;
      std::string output = "counterexamples.cert";
      for (int i = 2; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--out" && i + 1 < argc) output = argv[++i];
        else if (argv[i][0] != '-') limit = std::strtoull(argv[i], nullptr, 10);
        else usage();
      }
      std::ofstream certificate(output, std::ios::app);
      if (!certificate) throw std::runtime_error("cannot open output file");
      return explore_from_s14(limit, certificate);
    }
    if (cmd == "--hunt") {
      uint64_t samples = 400000;
      uint64_t seed = 1;
      std::string output = "counterexamples.cert";
      int i = 2;
      if (i < argc && argv[i][0] != '-') samples = std::strtoull(argv[i++], nullptr, 10);
      for (; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--out" && i + 1 < argc) output = argv[++i];
        else if (option == "--seed" && i + 1 < argc) seed = std::strtoull(argv[++i], nullptr, 10);
        else usage();
      }
      return hunt(samples, seed, output);
    }
    if (cmd != "--exhaustive") usage();
    if (argc < 3) usage();
    const int n = std::atoi(argv[2]);
    if (n < 1 || n > MAX_N) throw std::runtime_error("N must be in 1..15");
    std::string output = "counterexamples.cert";
    std::string checkpoint;
    bool symmetry = false;
    for (int i = 3; i < argc; ++i) {
      const std::string option = argv[i];
      if (option == "--out" && i + 1 < argc) output = argv[++i];
      else if (option == "--checkpoint" && i + 1 < argc) checkpoint = argv[++i];
      else if (option == "--symmetry") symmetry = true;
      else throw std::runtime_error("bad option: " + option);
    }
    return exhaustive(n, output, checkpoint, symmetry);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
