#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

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

struct Check {
  bool violation = false;
  int at = -1;
  int d = 0;
  int qn = 0;
  uint64_t q[MAX_DEG / 2 + 2]{};
};

static Check inspect(Perm u, Perm v, const Poly& p) {
  Check c;
  c.d = length(v) - length(u);
  int k = 0;
  for (int exponent = c.d & 1; exponent <= c.d; exponent += 2) c.q[k++] = p.a[exponent];
  c.qn = k;
  for (int i = 1; i + 1 < k; ++i) {
    if (static_cast<__uint128_t>(c.q[i]) * c.q[i] <
        static_cast<__uint128_t>(c.q[i - 1]) * c.q[i + 1]) {
      c.violation = true;
      c.at = i;
      break;
    }
  }
  return c;
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

static int exhaustive(int n, const std::string& output, const std::string& checkpoint) {
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
  std::atomic<uint64_t> intervals{0}, violations{0}, completed{0};
  const uint64_t total = factorial(N);
#pragma omp parallel for schedule(dynamic, 1)
  for (uint64_t v_rank = 0; v_rank < total; ++v_rank) {
    if (done.find(v_rank) != done.end()) continue;
    const Perm v = at_rank(v_rank);
    Recurrence recurrence;
    uint64_t local_intervals = 0, local_violations = 0;
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
        std::ostringstream line;
        if (report(u, v, recurrence.get(u, v), line)) {
          ++local_violations;
          std::lock_guard<std::mutex> lock(io_mutex);
          certificate << line.str();
          certificate.flush();
        }
      }
      const int c = down_covers(u, covers, MAX_N * MAX_N);
      for (int i = 0; i < c; ++i)
        if (seen.insert(covers[i]).second) stack.push_back(covers[i]);
    }
    intervals += local_intervals;
    violations += local_violations;
    ++completed;
    if (!checkpoint.empty()) {
      std::lock_guard<std::mutex> lock(io_mutex);
      checkpoint_file << v_rank << '\n';
      if (completed.load() % 32 == 0) checkpoint_file.flush();
    }
    if (completed.load() % 2000 == 0) {
      std::lock_guard<std::mutex> lock(io_mutex);
      std::cout << "progress v=" << completed.load() << "/" << total << " intervals=" << intervals
                << " violations=" << violations << std::endl;
    }
  }
  if (!checkpoint.empty()) {
    std::lock_guard<std::mutex> lock(io_mutex);
    checkpoint_file.flush();
  }
  std::cout << "completed v indices: " << completed << " / " << total << "; intervals=" << intervals
            << " violations=" << violations << '\n';
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

static void usage() {
  throw std::runtime_error(
      "usage: brenti_search --verify-s14 | --exhaustive N [--out FILE] [--checkpoint FILE] | "
      "--eval U V | --flatten [--out FILE] | --explore [LIMIT] [--out FILE] | "
      "--hunt [SAMPLES] [--out FILE] [--seed N]");
}

int main(int argc, char** argv) {
  try {
    if (argc < 2) usage();
    const std::string cmd = argv[1];
    if (cmd == "--verify-s14" || cmd == "--verify-s14") return verify_s14();
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
    for (int i = 3; i < argc; ++i) {
      const std::string option = argv[i];
      if (option == "--out" && i + 1 < argc) output = argv[++i];
      else if (option == "--checkpoint" && i + 1 < argc) checkpoint = argv[++i];
      else throw std::runtime_error("bad option: " + option);
    }
    return exhaustive(n, output, checkpoint);
  } catch (const std::exception& error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
