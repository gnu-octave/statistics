/*
Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This file is part of the statistics package for GNU Octave.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

#include <octave/oct.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// Dissimilarity measures for nominal data, after the R package nomclust
// 2.8.1 (Sulc and Rezankova), whose formulas and order of arithmetic are
// followed.  Level frequencies come from a reference sample R, so every
// per-level quantity is computed once and each pair costs O(K).
//
// Levels are codes 1..L per variable, shared by X, Y and R.  A level of Y
// that R does not hold has a count of 0; each measure answers it with the
// limit of its formula as that count falls to 0, and returns NaN where no
// such limit exists (IOF, and OF and Burnaby against a level filling R).
// The callers refuse those cases before they get here.

static const double Inf = std::numeric_limits<double>::infinity ();
static const double NaN = std::numeric_limits<double>::quiet_NaN ();

enum MeasureCode
{
  M_SM = 0,
  M_ESKIN,
  M_ANDERBERG,
  M_BURNABY,
  M_GAMBARYAN,
  M_GOODALL1,
  M_GOODALL2,
  M_GOODALL3,
  M_GOODALL4,
  M_IOF,
  M_LIN,
  M_LIN1,
  M_OF,
  M_SMIRNOV,
  M_VE,
  M_VM
};

// What the reference sample says about each level.  A level is addressed by
// its flat index, off[k] plus its zero-based code in variable k.
struct Stats
{
  double r;
  octave_idx_type K;
  std::vector<octave_idx_type> off;
  std::vector<double> f;
  std::vector<double> nlev;
  std::vector<double> w;
  double W;
  double sumlev;
  // Per level: the match term, and a second quantity some measures need
  std::vector<double> m;
  std::vector<double> q;
  // Per variable: a constant some measures need
  std::vector<double> c;
  // Lin1: each level's range of positions among its variable's levels
  // sorted by count, and the running sum of counts over them
  std::vector<octave_idx_type> first;
  std::vector<octave_idx_type> last;
  std::vector<std::vector<double>> P;
};

// The levels of variable k that R holds, in increasing order of count.
static std::vector<octave_idx_type>
sortedLevels (const Stats& S, octave_idx_type k)
{
  std::vector<octave_idx_type> idx;
  for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
  {
    if (S.f[v] > 0)
    {
      idx.push_back (v);
    }
  }
  std::stable_sort (idx.begin (), idx.end (),
                    [&S] (octave_idx_type a, octave_idx_type b)
                    { return S.f[a] < S.f[b]; });
  return idx;
}

// Every measure accumulates two sums over the variables of a pair and
// finishes them into a dissimilarity.  X is the side that R covers.
struct SM
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double&) const
  {
    if (x == y)
    {
      a += S.w[k];
    }
  }
  inline double finish (double a, double) const
  {
    return 1.0 - (1.0 / S.W * a);
  }
};

struct Eskin
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double&) const
  {
    a += (x == y) ? S.w[k] : S.w[k] * S.c[k];
  }
  inline double finish (double a, double) const
  {
    return (1.0 / (1.0 / S.W * a)) - 1.0;
  }
};

struct Anderberg
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double& b) const
  {
    const double px = S.f[x] / S.r;
    const double py = S.f[y] / S.r;
    if (x == y)
    {
      a += (1.0 / px) * (1.0 / px) * S.c[k];
    }
    else if (py == 0)
    {
      b = Inf;
    }
    else
    {
      b += 1.0 / 2.0 / px / py * S.c[k];
    }
  }
  inline double finish (double a, double b) const
  {
    return 1.0 - (a / (a + b));
  }
};

struct Burnaby
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double&) const
  {
    if (x == y)
    {
      a += S.w[k];
      return;
    }
    const double px = S.f[x] / S.r;
    const double py = S.f[y] / S.r;
    if (py == 0)
    {
      a += (px == 1) ? NaN : 0.0;
      return;
    }
    const double N = S.c[k];
    a += S.w[k] * (N / (std::log (px * py / (1.0 - px) / (1.0 - py)) + N));
  }
  inline double finish (double a, double) const
  {
    return 1.0 - (a / S.W);
  }
};

// Match terms only, over a divisor fixed for the measure.
struct MatchOnly
{
  const Stats& S;
  double div;
  bool weighted;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double&) const
  {
    if (x == y)
    {
      a += weighted ? S.w[k] * S.m[x] : S.m[x];
    }
  }
  inline double finish (double a, double) const
  {
    return 1.0 - (1.0 / div * a);
  }
};

struct IOF
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double&) const
  {
    if (x == y)
    {
      a += S.w[k];
    }
    else if (S.f[y] == 0)
    {
      a += NaN;
    }
    else
    {
      a += S.w[k] / (1.0 + S.m[x] * S.m[y]);
    }
  }
  inline double finish (double a, double) const
  {
    return (1.0 / (1.0 / S.W * a)) - 1.0;
  }
};

struct OF
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double&) const
  {
    if (x == y)
    {
      a += S.w[k];
    }
    else if (S.f[y] == 0)
    {
      a += (S.f[x] == S.r) ? NaN : 0.0;
    }
    else
    {
      a += S.w[k] / (1.0 + S.m[x] * S.m[y]);
    }
  }
  inline double finish (double a, double) const
  {
    if (a == 0)
    {
      return Inf;
    }
    return (1.0 / (1.0 / S.W * a)) - 1.0;
  }
};

// Lin and Lin1 divide the shared information by the total, which an unseen
// level makes infinite, so the pair is infinitely far apart.  Frequencies
// are summed as counts, so a sum covering every level is exactly 1 and its
// log exactly 0, where nomclust's rounding leaves about 1e16.  A total of 0
// needs every variable to hold one level in R, and the limit there is that
// of identical rows.
struct Lin
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double& b) const
  {
    if (S.f[y] == 0)
    {
      b = -Inf;
      return;
    }
    const double px = S.f[x] / S.r;
    const double py = S.f[y] / S.r;
    if (x == y)
    {
      a += S.w[k] * (2.0 * std::log (px));
    }
    else
    {
      a += S.w[k] * (2.0 * std::log ((S.f[x] + S.f[y]) / S.r));
    }
    b += (std::log (px) + std::log (py)) * S.w[k];
  }
  inline double finish (double a, double b) const
  {
    if (std::isinf (b))
    {
      return Inf;
    }
    if (b == 0)
    {
      return 0.0;
    }
    const double s = 1.0 / b * a;
    return (s == 0) ? Inf : (1.0 / s) - 1.0;
  }
};

struct Lin1
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double& b) const
  {
    if (S.f[y] == 0)
    {
      b = -Inf;
      return;
    }
    const double px = S.f[x] / S.r;
    const double py = S.f[y] / S.r;
    if (x == y)
    {
      a += S.w[k] * std::log (S.m[x]);
    }
    else
    {
      octave_idx_type lo = (S.f[x] <= S.f[y]) ? x : y;
      octave_idx_type hi = (S.f[x] <= S.f[y]) ? y : x;
      const std::vector<double>& P = S.P[k];
      a += S.w[k] * (2.0 * std::log ((P[S.last[hi] + 1] - P[S.first[lo]])
                                     / S.r));
    }
    b += (std::log (px) + std::log (py)) * S.w[k];
  }
  // With every variable holding one level, each match term is half its
  // weight in the limit, so identical rows sit at 1 as they do elsewhere.
  inline double finish (double a, double b) const
  {
    if (std::isinf (b))
    {
      return Inf;
    }
    if (b == 0)
    {
      return 1.0;
    }
    const double s = 1.0 / b * a;
    return (s == 0) ? Inf : (1.0 / s) - 1.0;
  }
};

struct Smirnov
{
  const Stats& S;
  inline void add (octave_idx_type k, octave_idx_type x, octave_idx_type y,
                   double& a, double&) const
  {
    if (x == y)
    {
      a += 2.0 + (S.r - S.f[x]) / S.f[x] + (S.c[k] - S.q[x]);
    }
    else
    {
      a += std::max (0.0, S.c[k] - S.q[x] - S.q[y]);
    }
  }
  inline double finish (double a, double) const
  {
    return 1.0 / (1.0 / S.sumlev * a + 1.0);
  }
};

// The per-level and per-variable quantities each measure reads.
static void
prepare (Stats& S, int measure)
{
  const octave_idx_type K = S.K;
  const octave_idx_type nl = S.off[K];
  const double r = S.r;
  S.m.assign (nl, 0.0);
  S.q.assign (nl, 0.0);
  S.c.assign (K, 0.0);

  // The chance of drawing a level twice without replacement; a reference of
  // one row cannot draw twice
  auto p2 = [r] (double f)
  {
    return (r > 1) ? f * (f - 1.0) / r / (r - 1.0) : 0.0;
  };

  for (octave_idx_type k = 0; k < K; k++)
  {
    const double n = S.nlev[k];
    switch (measure)
    {
      case M_ESKIN:
      {
        S.c[k] = (n * n) / ((n * n) + 2.0);
        break;
      }
      case M_ANDERBERG:
      {
        S.c[k] = 2.0 / n / (n + 1.0);
        break;
      }
      case M_BURNABY:
      {
        double N = 0.0;
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          N = N + (2.0 * std::log (1.0 - S.f[v] / r));
        }
        S.c[k] = N;
        break;
      }
      case M_GAMBARYAN:
      {
        // A level held by every row has entropy 0, not 0 * log2 (0)
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          const double p = S.f[v] / r;
          if (p > 0 && p < 1)
          {
            S.m[v] = - (p * std::log2 (p) + (1.0 - p) * std::log2 (1.0 - p));
          }
        }
        break;
      }
      case M_GOODALL1:
      case M_GOODALL2:
      {
        // Levels at most (Goodall 1) or at least (Goodall 2) as frequent,
        // as a running sum over the levels sorted by count
        std::vector<octave_idx_type> idx = sortedLevels (S, k);
        const std::size_t L = idx.size ();
        std::vector<double> T (L + 1, 0.0);
        for (std::size_t i = 0; i < L; i++)
        {
          T[i+1] = T[i] + p2 (S.f[idx[i]]);
        }
        std::size_t i = 0;
        while (i < L)
        {
          std::size_t j = i;
          while (j < L && S.f[idx[j]] == S.f[idx[i]])
          {
            j++;
          }
          const double s = (measure == M_GOODALL1) ? T[j] : T[L] - T[i];
          for (std::size_t u = i; u < j; u++)
          {
            S.m[idx[u]] = 1.0 * (1.0 - s);
          }
          i = j;
        }
        break;
      }
      case M_GOODALL3:
      case M_GOODALL4:
      {
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          S.m[v] = (measure == M_GOODALL3) ? 1.0 - p2 (S.f[v]) : p2 (S.f[v]);
        }
        break;
      }
      case M_IOF:
      {
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          S.m[v] = std::log (S.f[v]);
        }
        break;
      }
      case M_OF:
      {
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          S.m[v] = std::log (r / S.f[v]);
        }
        break;
      }
      case M_LIN1:
      {
        std::vector<octave_idx_type> idx = sortedLevels (S, k);
        const std::size_t L = idx.size ();
        std::vector<double> P (L + 1, 0.0);
        for (std::size_t i = 0; i < L; i++)
        {
          P[i+1] = P[i] + S.f[idx[i]];
        }
        std::size_t i = 0;
        while (i < L)
        {
          std::size_t j = i;
          double t = 0.0;
          while (j < L && S.f[idx[j]] == S.f[idx[i]])
          {
            t = t + S.f[idx[j]];
            j++;
          }
          for (std::size_t u = i; u < j; u++)
          {
            S.first[idx[u]] = i;
            S.last[idx[u]] = j - 1;
            S.m[idx[u]] = t / r;
          }
          i = j;
        }
        S.P[k] = P;
        break;
      }
      case M_SMIRNOV:
      {
        // A level held by every row is the only one, so never among the
        // others any sum here runs over
        double T = 0.0;
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          if (S.f[v] < r)
          {
            S.q[v] = S.f[v] / (r - S.f[v]);
          }
          T = T + S.q[v];
        }
        S.c[k] = T;
        break;
      }
      case M_VE:
      {
        double e = 0.0;
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          const double p = S.f[v] / r;
          if (p > 0)
          {
            e = e + p * std::log (p);
          }
        }
        const double ne = (n > 1) ? - e / std::log (n) : 0.0;
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          S.m[v] = ne;
        }
        break;
      }
      case M_VM:
      {
        double g = 0.0;
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          const double p = S.f[v] / r;
          g = g + p * p;
        }
        const double ng = (n > 1) ? (1.0 - g) * n / (n - 1.0) : 0.0;
        for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
        {
          S.m[v] = ng;
        }
        break;
      }
      default:
        break;
    }
  }
}

// Every pair within X in pdist order, or every row of X against every row
// of Y.  Codes arrive as flat level indices, one row contiguous.
template <typename F>
static void
pairs (const F& fn, const std::vector<octave_idx_type>& xc,
       const std::vector<octave_idx_type>& yc, octave_idx_type n,
       octave_idx_type m, octave_idx_type K, bool within,
       const std::vector<octave_idx_type>& act, double *out)
{
  std::size_t at = 0;
  for (octave_idx_type i = 0; i < n; i++)
  {
    octave_quit ();
    const octave_idx_type *xp = &xc[(std::size_t) i * K];
    const octave_idx_type j0 = within ? i + 1 : 0;
    const octave_idx_type j1 = within ? n : m;
    for (octave_idx_type j = j0; j < j1; j++)
    {
      const octave_idx_type *yp = within ? &xc[(std::size_t) j * K]
                                         : &yc[(std::size_t) j * K];
      double a = 0.0;
      double b = 0.0;
      for (octave_idx_type k : act)
      {
        fn.add (k, xp[k], yp[k], a, b);
      }
      if (within)
      {
        out[at++] = fn.finish (a, b);
      }
      else
      {
        out[i + j * n] = fn.finish (a, b);
      }
    }
  }
}

// Codes as flat level indices, one row contiguous, checked on the way.
static std::vector<octave_idx_type>
pack (const Matrix& A, const Stats& S)
{
  const octave_idx_type n = A.rows ();
  const octave_idx_type K = S.K;
  std::vector<octave_idx_type> out ((std::size_t) n * K);
  for (octave_idx_type k = 0; k < K; k++)
  {
    for (octave_idx_type i = 0; i < n; i++)
    {
      const double v = A(i, k);
      out[(std::size_t) i * K + k] = S.off[k] + (octave_idx_type) v - 1;
    }
  }
  return out;
}

static bool
validCodes (const Matrix& A)
{
  const octave_idx_type N = A.numel ();
  for (octave_idx_type i = 0; i < N; i++)
  {
    const double v = A.xelem (i);
    if (! (v >= 1) || v != std::floor (v) || ! std::isfinite (v))
    {
      return false;
    }
  }
  return true;
}

DEFUN_DLD(__nomdist__, args, ,
"-*- texinfo -*-\n\
@deftypefn  {statistics} {@var{D} =} __nomdist__ (@var{X}, [], @var{R}, @var{measure}, @var{w})\n\
@deftypefnx {statistics} {@var{D} =} __nomdist__ (@var{X}, @var{Y}, @var{R}, @var{measure}, @var{w})\n\
\n\
Dissimilarity measures for nominal data.  Internal; called by @code{nomdist}\n\
and @code{nomdist2} and not meant to be used directly.\n\
\n\
@var{X}, @var{Y} and @var{R} hold level codes, positive integers shared\n\
across the three, one variable per column.  Level frequencies are counted\n\
over @var{R}, which must hold every level @var{X} holds, or over @var{X}\n\
itself where @var{R} is empty.  @var{measure} is one of\n\
@qcode{'anderberg'}, @qcode{'burnaby'}, @qcode{'eskin'},\n\
@qcode{'gambaryan'}, @qcode{'goodall1'}, @qcode{'goodall2'},\n\
@qcode{'goodall3'}, @qcode{'goodall4'}, @qcode{'iof'}, @qcode{'lin'},\n\
@qcode{'lin1'}, @qcode{'of'}, @qcode{'sm'}, @qcode{'smirnov'},\n\
@qcode{'ve'} or @qcode{'vm'}, and @var{w} holds one weight per variable,\n\
read by every measure but @qcode{'anderberg'}, @qcode{'gambaryan'} and\n\
@qcode{'smirnov'}.\n\
\n\
With @var{Y} empty, @var{D} is a row vector over every pair of rows of\n\
@var{X}, in the order @code{pdist} returns them.  Otherwise @var{D} is\n\
the @math{N*M} matrix of each row of @var{X} against each row of @var{Y}.\n\
\n\
@end deftypefn")
{
  if (args.length () != 5)
  {
    print_usage ();
  }

  for (int i = 0; i < 3; i++)
  {
    if (! args(i).isnumeric () || args(i).iscomplex ()
        || args(i).ndims () != 2)
    {
      error ("__nomdist__: X, Y and R must be real numeric matrices.");
    }
  }
  const Matrix X = args(0).matrix_value ();
  const Matrix Y = args(1).matrix_value ();
  const Matrix R = args(2).isempty () ? X : args(2).matrix_value ();
  const bool within = Y.isempty ();
  const octave_idx_type K = X.columns ();
  if (R.columns () != K || (! within && Y.columns () != K))
  {
    error ("__nomdist__: X, Y and R must have the same number of columns.");
  }
  if (! validCodes (X) || ! validCodes (Y) || ! validCodes (R))
  {
    error ("__nomdist__: X, Y and R must hold positive integer codes.");
  }
  if (! args(3).is_string ())
  {
    error ("__nomdist__: MEASURE must be a character vector.");
  }
  if (! args(4).isnumeric () || args(4).numel () != K)
  {
    error ("__nomdist__: W must hold one weight per variable.");
  }

  static const char *names[] = {"sm", "eskin", "anderberg", "burnaby",
                                "gambaryan", "goodall1", "goodall2",
                                "goodall3", "goodall4", "iof", "lin",
                                "lin1", "of", "smirnov", "ve", "vm"};
  const std::string name = args(3).string_value ();
  int measure = -1;
  for (int i = 0; i < 16; i++)
  {
    if (name == names[i])
    {
      measure = i;
    }
  }
  if (measure < 0)
  {
    error ("__nomdist__: unsupported MEASURE '%s'.", name.c_str ());
  }

  // Levels per variable, over all three sets, and their counts in R
  Stats S;
  S.K = K;
  S.r = (double) R.rows ();
  S.off.assign (K + 1, 0);
  for (octave_idx_type k = 0; k < K; k++)
  {
    double L = 0;
    for (octave_idx_type i = 0; i < X.rows (); i++)
    {
      L = std::max (L, X(i, k));
    }
    for (octave_idx_type i = 0; i < Y.rows (); i++)
    {
      L = std::max (L, Y(i, k));
    }
    for (octave_idx_type i = 0; i < R.rows (); i++)
    {
      L = std::max (L, R(i, k));
    }
    S.off[k+1] = S.off[k] + (octave_idx_type) L;
  }
  S.f.assign (S.off[K], 0.0);
  S.nlev.assign (K, 0.0);
  for (octave_idx_type k = 0; k < K; k++)
  {
    for (octave_idx_type i = 0; i < R.rows (); i++)
    {
      S.f[S.off[k] + (octave_idx_type) R(i, k) - 1] += 1.0;
    }
    for (octave_idx_type v = S.off[k]; v < S.off[k+1]; v++)
    {
      S.nlev[k] += (S.f[v] > 0) ? 1.0 : 0.0;
    }
  }
  S.sumlev = 0.0;
  for (octave_idx_type k = 0; k < K; k++)
  {
    S.sumlev += S.nlev[k];
  }

  const NDArray wv = args(4).array_value ();
  S.w.assign (K, 1.0);
  const bool weighted = (measure != M_ANDERBERG && measure != M_GAMBARYAN
                         && measure != M_SMIRNOV);
  if (weighted)
  {
    for (octave_idx_type k = 0; k < K; k++)
    {
      S.w[k] = wv(k);
    }
  }
  S.W = 0.0;
  for (octave_idx_type k = 0; k < K; k++)
  {
    S.W += S.w[k];
  }

  // A variable of weight 0 adds nothing, and skipping it keeps an unseen
  // level there from reaching a log
  std::vector<octave_idx_type> act;
  for (octave_idx_type k = 0; k < K; k++)
  {
    if (S.w[k] != 0)
    {
      act.push_back (k);
    }
  }

  const std::vector<octave_idx_type> xc = pack (X, S);
  const std::vector<octave_idx_type> yc = pack (Y, S);
  for (octave_idx_type v : xc)
  {
    if (S.f[v] == 0)
    {
      error ("__nomdist__: R must hold every level X holds.");
    }
  }

  if (measure == M_LIN1)
  {
    S.first.assign (S.off[K], 0);
    S.last.assign (S.off[K], 0);
    S.P.assign (K, std::vector<double> ());
  }
  prepare (S, measure);

  const octave_idx_type n = X.rows ();
  const octave_idx_type m = Y.rows ();
  NDArray D;
  if (within)
  {
    D = NDArray (dim_vector (1, n > 1 ? n * (n - 1) / 2 : 0));
  }
  else
  {
    D = NDArray (dim_vector (n, m));
  }
  double *out = D.fortran_vec ();

  switch (measure)
  {
    case M_SM:
      pairs (SM {S}, xc, yc, n, m, K, within, act, out);
      break;
    case M_ESKIN:
      pairs (Eskin {S}, xc, yc, n, m, K, within, act, out);
      break;
    case M_ANDERBERG:
      pairs (Anderberg {S}, xc, yc, n, m, K, within, act, out);
      break;
    case M_BURNABY:
      pairs (Burnaby {S}, xc, yc, n, m, K, within, act, out);
      break;
    case M_GAMBARYAN:
      pairs (MatchOnly {S, S.sumlev, false}, xc, yc, n, m, K, within, act,
             out);
      break;
    case M_GOODALL1:
    case M_GOODALL2:
    case M_GOODALL3:
    case M_GOODALL4:
    case M_VE:
    case M_VM:
      pairs (MatchOnly {S, S.W, true}, xc, yc, n, m, K, within, act, out);
      break;
    case M_IOF:
      pairs (IOF {S}, xc, yc, n, m, K, within, act, out);
      break;
    case M_LIN:
      pairs (Lin {S}, xc, yc, n, m, K, within, act, out);
      break;
    case M_LIN1:
      pairs (Lin1 {S}, xc, yc, n, m, K, within, act, out);
      break;
    case M_OF:
      pairs (OF {S}, xc, yc, n, m, K, within, act, out);
      break;
    case M_SMIRNOV:
      pairs (Smirnov {S}, xc, yc, n, m, K, within, act, out);
      break;
  }

  return ovl (D);
}

/*
%!shared X, W
%! X = [1, 1, 1; 1, 2, 1; 1, 1, 2; 2, 2, 2; 2, 1, 1; 3, 2, 2; 4, 1, 3];
%! W = [0.7, 1, 0.4];

## Expectations from nomclust 2.8.1 on the same data
%!test
%! D = [0.3333333333333334, 0.3333333333333334, 1, 0.3333333333333334, 1, ...
%!      0.6666666666666667, 0.6666666666666667, 0.6666666666666667, ...
%!      0.6666666666666667, 0.6666666666666667, 1, 0.6666666666666667, ...
%!      0.6666666666666667, 0.6666666666666667, 0.6666666666666667, ...
%!      0.6666666666666667, 0.3333333333333334, 1, 1, 0.6666666666666667, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'sm', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.125, 0.06451612903225801, 0.2638297872340427, 0.03846153846153855, ...
%!      0.2638297872340427, 0.1082089552238805, 0.2073170731707317, ...
%!      0.1082089552238805, 0.173913043478261, 0.1082089552238805, ...
%!      0.2638297872340427, 0.173913043478261, 0.1082089552238805, ...
%!      0.173913043478261, 0.1082089552238805, 0.2073170731707317, ...
%!      0.03846153846153855, 0.2638297872340427, 0.2638297872340427, ...
%!      0.1082089552238805, 0.2638297872340427];
%! assert_equal (__nomdist__ (X, [], X, 'eskin', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.3191489361702128, 0.2247191011235955, 1, 0.174757281553398, 1, ...
%!      0.6808510638297873, 0.6756756756756757, 0.3220338983050848, ...
%!      0.5454545454545454, 0.4117647058823529, 1, 0.5454545454545454, ...
%!      0.4578313253012049, 0.6226415094339623, 0.6808510638297873, ...
%!      0.4807692307692308, 0.3103448275862069, 1, 1, 0.7169811320754718, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'anderberg', ones (1, 3)), D, -1e-14);
%!test
%! D = [0, 0.06142861791239795, 0.1725141049787992, 0.1110854870664012, ...
%!      0.2158655957496493, 0.3042675674408121, 0.06142861791239795, ...
%!      0.1725141049787992, 0.1110854870664012, 0.2158655957496493, ...
%!      0.3042675674408121, 0.1110854870664012, 0.1725141049787992, ...
%!      0.1544369778372513, 0.3042675674408121, 0.06142861791239795, ...
%!      0.1764146125938005, 0.3262452021973615, 0.2378432305061988, ...
%!      0.3262452021973615, 0.3491708712867553];
%! assert_equal (__nomdist__ (X, [], X, 'burnaby', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.7810604142146108, 0.7810604142146108, 1, 0.7810604142146108, 1, ...
%!      0.8905302071073053, 0.8905302071073053, 0.8905302071073053, ...
%!      0.8905302071073053, 0.8905302071073053, 1, 0.8905302071073053, ...
%!      0.8905302071073053, 0.8905302071073053, 0.8905302071073053, ...
%!      0.9040977146037077, 0.7810604142146108, 1, 1, 0.8905302071073053, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'gambaryan', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.4920634920634921, 0.5396825396825398, 1, 0.5714285714285715, 1, ...
%!      0.8095238095238095, 0.7301587301587302, 0.7142857142857143, ...
%!      0.7619047619047619, 0.7142857142857143, 1, 0.7619047619047619, ...
%!      0.8095238095238095, 0.7619047619047619, 0.8095238095238095, ...
%!      0.6825396825396826, 0.4761904761904762, 1, 1, 0.8095238095238095, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'goodall1', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.4761904761904762, 0.4761904761904762, 1, 0.5238095238095238, 1, ...
%!      0.7619047619047619, 0.7142857142857143, 0.8095238095238095, ...
%!      0.7619047619047619, 0.8095238095238095, 1, 0.7619047619047619, ...
%!      0.7619047619047619, 0.7619047619047619, 0.7619047619047619, ...
%!      0.7301587301587302, 0.5714285714285715, 1, 1, 0.7619047619047619, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'goodall2', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.4285714285714286, 0.4761904761904762, 1, 0.4761904761904762, 1, ...
%!      0.7619047619047619, 0.7142857142857143, 0.7142857142857143, ...
%!      0.7142857142857143, 0.7142857142857143, 1, 0.7142857142857143, ...
%!      0.7619047619047619, 0.7142857142857143, 0.7619047619047619, ...
%!      0.6825396825396826, 0.4285714285714286, 1, 1, 0.7619047619047619, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'goodall3', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.9047619047619048, 0.8571428571428572, 1, 0.8571428571428572, 1, ...
%!      0.9047619047619048, 0.9523809523809523, 0.9523809523809523, ...
%!      0.9523809523809523, 0.9523809523809523, 1, 0.9523809523809523, ...
%!      0.9047619047619048, 0.9523809523809523, 0.9047619047619048, ...
%!      0.9841269841269842, 0.9047619047619048, 1, 1, 0.9047619047619048, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'goodall4', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.2519020857326395, 0.222935300643843, 1.116901264683761, ...
%!      0.1683617083596181, 0.6220882713463793, 0, 0.6220882713463793, ...
%!      0.4845515922241048, 0.5274548356775992, 0.222935300643843, ...
%!      0.2519020857326395, 0.5274548356775992, 0.4845515922241048, ...
%!      0.2519020857326395, 0, 0.6220882713463793, 0, 0.2519020857326395, ...
%!      0.6220882713463793, 0, 0.2519020857326395];
%! assert_equal (__nomdist__ (X, [], X, 'iof', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.4151177862291795, 0.4440221764476788, 4.300985770938718, ...
%!      0.4092944562546224, 3.129303942378162, 0.9970986461666451, ...
%!      1.394583893868549, 1.051411550470893, 1.197035645318646, ...
%!      0.9801872797639546, 2.124165640951619, 1.197035645318646, ...
%!      1.33941488968952, 1.094910865911154, 0.9970986461666451, ...
%!      0.9926721560959628, 0.2958576645228652, 1.629441680424714, ...
%!      2.145534809656081, 0.8080361706650088, 1.23240918245247];
%! assert_equal (__nomdist__ (X, [], X, 'lin', ones (1, 3)), D, -1e-14);
%!test
%! ## Pair 11 shares nothing, so its similarity is exactly 0, where
%! ## nomclust's rounding returns about 1.6e16.
%! D = [3.789167787737098, 1.628488554759508, 4.300985770938718, ...
%!      2.543556180472745, 18.11998502097161, 10.98259187699986, ...
%!      3.150318732214346, 2.001980368635569, 5.288962253828324, ...
%!      4.34995121471969, Inf, 5.288962253828324, 2.189052189432408, ...
%!      37.23997004194317, 10.98259187699986, 2.591802852053521, ...
%!      2.106486692231752, 5.61060225142327, 3.413335993771041, ...
%!      3.235711272013311, 2.229638071755347];
%! assert_equal (__nomdist__ (X, [], X, 'lin1', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.1200918256837458, 0.1618443665773319, 0.7186603816400106, ...
%!      0.2071986251911722, 0.8315156613470343, 0.7093347701634469, ...
%!      0.3271674504213862, 0.4512427967395944, 0.3866778238431001, ...
%!      0.5308967156237416, 1.092895882914944, 0.3866778238431001, ...
%!      0.4512427967395944, 0.4592247360215442, 0.7093347701634469, ...
%!      0.3271674504213862, 0.3095365880953216, 1.227546987265509, ...
%!      0.933812054699924, 0.7981072047001967, 1.371908574090516];
%! assert_equal (__nomdist__ (X, [], X, 'of', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.519730510105871, 0.5378486055776893, 0.9473684210526315, ...
%!      0.526829268292683, 0.9246575342465754, 0.6513872135102533, ...
%!      0.6801007556675063, 0.6352941176470589, 0.6625766871165645, 0.625, ...
%!      0.8723747980613893, 0.6625766871165645, 0.6923076923076923, ...
%!      0.6513872135102533, 0.6513872135102533, 0.6101694915254238, ...
%!      0.4778761061946903, 0.84375, 0.8925619834710745, 0.6352941176470589, ...
%!      0.8256880733944955];
%! assert_equal (__nomdist__ (X, [], X, 'smirnov', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.3882378704642936, 0.3645287891257314, 1, 0.3668903239823945, 1, ...
%!      0.6715906213219162, 0.6929381678038152, 0.6715906213219162, ...
%!      0.6952997026604784, 0.6715906213219162, 1, 0.6952997026604784, ...
%!      0.6715906213219162, 0.6952997026604784, 0.6715906213219162, ...
%!      0.6929381678038152, 0.3668903239823945, 1, 1, 0.6715906213219162, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 've', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.3854875283446713, 0.365079365079365, 1, 0.3673469387755102, 1, ...
%!      0.673469387755102, 0.691609977324263, 0.673469387755102, ...
%!      0.6938775510204082, 0.673469387755102, 1, 0.6938775510204082, ...
%!      0.673469387755102, 0.6938775510204082, 0.673469387755102, ...
%!      0.691609977324263, 0.3673469387755102, 1, 1, 0.673469387755102, 1];
%! assert_equal (__nomdist__ (X, [], X, 'vm', ones (1, 3)), D, -1e-14);
%!test
%! D = [0.4761904761904762, 0.1904761904761906, 1, 0.3333333333333334, 1, ...
%!      0.5238095238095238, 0.6666666666666667, 0.5238095238095238, ...
%!      0.8095238095238095, 0.5238095238095238, 1, 0.8095238095238095, ...
%!      0.5238095238095238, 0.8095238095238095, 0.5238095238095238, ...
%!      0.6666666666666667, 0.3333333333333334, 1, 1, 0.5238095238095238, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'sm', W), D, -1e-14);
%!test
%! D = [0.188679245283019, 0.03587443946188351, 0.2993750000000002, ...
%!      0.03846153846153855, 0.2993750000000002, 0.07720207253886002, ...
%!      0.2397137745974953, 0.07720207253886002, 0.2434210526315792, ...
%!      0.07720207253886002, 0.2993750000000002, 0.2434210526315792, ...
%!      0.07720207253886002, 0.2434210526315792, 0.07720207253886002, ...
%!      0.2397137745974953, 0.03846153846153855, 0.2993750000000002, ...
%!      0.2993750000000002, 0.07720207253886002, 0.2993750000000002];
%! assert_equal (__nomdist__ (X, [], X, 'eskin', W), D, -1e-14);
%!test
%! D = [2.220446049250313e-16, 0.03510206737851329, 0.1461875544449145, ...
%!      0.1110854870664013, 0.1895390452157646, 0.2400544576107146, ...
%!      0.03510206737851329, 0.1461875544449145, 0.1110854870664013, ...
%!      0.1895390452157645, 0.2400544576107146, 0.1110854870664013, ...
%!      0.1461875544449145, 0.1544369778372515, 0.2400544576107146, ...
%!      0.03510206737851329, 0.1764146125938008, 0.262032092367264, ...
%!      0.2115166799723139, 0.262032092367264, 0.2849577614566577];
%! assert_equal (__nomdist__ (X, [], X, 'burnaby', W), D, -1e-14);
%!test
%! D = [0.5941043083900226, 0.4580498866213153, 1, 0.5918367346938775, 1, ...
%!      0.7278911564625851, 0.7301587301587302, 0.5918367346938775, ...
%!      0.8639455782312925, 0.5918367346938775, 1, 0.8639455782312925, ...
%!      0.7278911564625851, 0.8639455782312925, 0.7278911564625851, ...
%!      0.6825396825396826, 0.45578231292517, 1, 1, 0.7278911564625851, 1];
%! assert_equal (__nomdist__ (X, [], X, 'goodall1', W), D, -1e-14);
%!test
%! D = [0.5782312925170068, 0.3741496598639457, 1, 0.5238095238095238, 1, ...
%!      0.6598639455782314, 0.7142857142857143, 0.7278911564625851, ...
%!      0.8639455782312925, 0.7278911564625851, 1, 0.8639455782312925, ...
%!      0.6598639455782314, 0.8639455782312925, 0.6598639455782314, ...
%!      0.7301587301587302, 0.5918367346938775, 1, 1, 0.6598639455782314, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'goodall2', W), D, -1e-14);
%!test
%! D = [0.5510204081632654, 0.3741496598639457, 1, 0.4965986394557823, 1, ...
%!      0.6598639455782314, 0.7142857142857143, 0.5918367346938775, ...
%!      0.8367346938775511, 0.5918367346938775, 1, 0.8367346938775511, ...
%!      0.6598639455782314, 0.8367346938775511, 0.6598639455782314, ...
%!      0.6825396825396826, 0.4285714285714285, 1, 1, 0.6598639455782314, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'goodall3', W), D, -1e-14);
%!test
%! D = [0.9251700680272109, 0.8163265306122449, 1, 0.8367346938775511, 1, ...
%!      0.8639455782312926, 0.9523809523809524, 0.9319727891156463, ...
%!      0.9727891156462585, 0.9319727891156463, 1, 0.9727891156462585, ...
%!      0.8639455782312926, 0.9727891156462585, 0.8639455782312926, ...
%!      0.9841269841269842, 0.9047619047619048, 1, 1, 0.8639455782312926, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 'goodall4', W), D, -1e-14);
%!test
%! D = [0.4034116524755296, 0.1162816237598046, 1.153873045642792, ...
%!      0.1683617083596181, 0.6437079284563265, 0, 0.6437079284563265, ...
%!      0.3302637748244865, 0.7591738998286246, 0.1162816237598046, ...
%!      0.4034116524755296, 0.7591738998286246, 0.3302637748244865, ...
%!      0.4034116524755296, 0, 0.6437079284563265, 0, 0.4034116524755296, ...
%!      0.6437079284563265, 0, 0.4034116524755296];
%! assert_equal (__nomdist__ (X, [], X, 'iof', W), D, -1e-14);
%!test
%! D = [0.7547596113257813, 0.2283122503928992, 4.98065967640019, ...
%!      0.4404425718584952, 3.455309722495803, 0.7834346505595069, ...
%!      1.497804622067723, 0.678686996286928, 2.094086789730795, ...
%!      0.663610211168179, 2.638407372343877, 2.094086789730795, ...
%!      0.9065634806963383, 1.764659180935495, 0.7834346505595069, ...
%!      1.044877050999849, 0.2958576645228654, 1.915257938793264, ...
%!      2.301793922738517, 0.6256289159595232, 1.383977355327267];
%! assert_equal (__nomdist__ (X, [], X, 'lin', W), D, -1e-14);
%!test
%! ## Pair 11 shares nothing, so its similarity is exactly 0, where
%! ## nomclust's rounding returns about 1.6e16.
%! D = [3.995609244135445, 1.337916967393231, 4.98065967640019, ...
%!      1.990937661841057, 31.76013401116178, 6.490425532349926, ...
%!      3.565656513831422, 1.66527788918046, 5.672895993089274, ...
%!      3.458682832253068, Inf, 5.672895993089274, 1.831126193046863, ...
%!      64.52026802232339, 6.490425532349926, 2.837642944938641, ...
%!      1.72479913371176, 5.079768268917555, 3.768312597687203, ...
%!      2.332335370699156, 1.992511747898871];
%! assert_equal (__nomdist__ (X, [], X, 'lin1', W), D, -1e-14);
%!test
%! D = [0.1808686870497029, 0.08648381587311493, 0.6789831874074297, ...
%!      0.2071986251911724, 0.7865248082245031, 0.4837971469738953, ...
%!      0.3033825055232831, 0.3355336337184947, 0.4810463690758979, ...
%!      0.4026982054804156, 0.9201937519825203, 0.4810463690758979, ...
%!      0.3355336337184947, 0.5640991370603945, 0.4837971469738953, ...
%!      0.3033825055232831, 0.3095365880953216, 1.032940807336712, ...
%!      0.8837237854764428, 0.5502337057652003, 1.152503424548615];
%! assert_equal (__nomdist__ (X, [], X, 'of', W), D, -1e-14);
%!test
%! D = [0.5188237121812314, 0.2237819125494098, 1, 0.3567292891230108, 1, ...
%!      0.5308437447455945, 0.6929381678038152, 0.5308437447455945, ...
%!      0.8258855443774162, 0.5308437447455945, 1, 0.8258855443774162, ...
%!      0.5308437447455945, 0.8258855443774162, 0.5308437447455945, ...
%!      0.6929381678038152, 0.3567292891230108, 1, 1, 0.5308437447455945, ...
%!      1];
%! assert_equal (__nomdist__ (X, [], X, 've', W), D, -1e-14);
%!test
%! D = [0.5166828636216392, 0.2251376741172659, 1, 0.358600583090379, 1, ...
%!      0.5335276967930028, 0.691609977324263, 0.5335276967930028, ...
%!      0.8250728862973761, 0.5335276967930028, 1, 0.8250728862973761, ...
%!      0.5335276967930028, 0.8250728862973761, 0.5335276967930028, ...
%!      0.691609977324263, 0.358600583090379, 1, 1, 0.5335276967930028, 1];
%! assert_equal (__nomdist__ (X, [], X, 'vm', W), D, -1e-14);
%!test
%! ## Row against row matches the pairs, bit for bit
%! D = squareform (__nomdist__ (X, [], X, 'goodall3', ones (1, 3)));
%! D2 = __nomdist__ (X, X, X, 'goodall3', ones (1, 3));
%! assert_equal (D2(! eye (7)), D(! eye (7)));
%!test
%! ## Identical rows need not be at 0: 1 - (3 - 24/42) / 3
%! D2 = __nomdist__ (X, X, X, 'goodall3', ones (1, 3));
%! assert_equal (D2(1,1), 4 / 21, -1e-14);
%!assert_equal (__nomdist__ ([1; 2], [1], [1; 2], 'goodall3', 1), [0; 1])
%!assert_equal (__nomdist__ ([1; 2], [1], [], 'goodall3', 1), [0; 1])
%!assert_equal (__nomdist__ (X, [], [], 'lin1', W), ...
%!              __nomdist__ (X, [], X, 'lin1', W))
%!assert_equal (__nomdist__ ([1; 2], [1], [1; 1; 1; 2], 'goodall3', 1), ...
%!              [0.5; 1])
%!assert_equal (__nomdist__ ([1 1 1], [], [1 1 1], 'sm', ones (1, 3)), ...
%!              zeros (1, 0))
%!assert_equal (__nomdist__ (X, [], X, 'sm', [1, 0, 1]), ...
%!              __nomdist__ (X(:,[1, 3]), [], X(:,[1, 3]), 'sm', [1, 1]))
## A level of Y that R does not hold
%!assert_equal (__nomdist__ ([1, 1; 2, 1; 2, 2], [3, 1; 3, 2], ...
%!                           [1, 1; 2, 1; 2, 2], 'goodall3', [1, 1]), ...
%!              [2/3, 1; 2/3, 1; 1, 0.5], -1e-14)
%!test
%! Z = [1, 1; 2, 1; 2, 2];
%! s = 1 / (1 + log (3) * log (3 / 2));
%! D = __nomdist__ (Z, [3, 1; 3, 2], Z, 'of', [1, 1]);
%! assert_equal (D, [1, 2/s-1; 1, 2/s-1; 2/s-1, 1], -1e-14);
%!assert_equal (__nomdist__ ([1, 1; 2, 1; 2, 2], [3, 1; 3, 2], ...
%!                           [1, 1; 2, 1; 2, 2], 'burnaby', [1, 1]), ...
%!              0.5 * ones (3, 2))
%!assert_equal (__nomdist__ ([1, 1; 2, 1; 2, 2], [3, 1; 3, 2], ...
%!                           [1, 1; 2, 1; 2, 2], 'anderberg', [1, 1]), ...
%!              ones (3, 2))
%!assert_equal (__nomdist__ ([1, 1; 2, 1; 2, 2], [3, 1; 3, 2], ...
%!                           [1, 1; 2, 1; 2, 2], 'lin', [1, 1]), ...
%!              Inf (3, 2))
%!assert_equal (__nomdist__ ([1, 1; 2, 1; 2, 2], [3, 1; 3, 2], ...
%!                           [1, 1; 2, 1; 2, 2], 'lin1', [1, 1]), ...
%!              Inf (3, 2))
%!assert_equal (__nomdist__ ([1, 1; 2, 1; 2, 2], [3, 1; 3, 2], ...
%!                           [1, 1; 2, 1; 2, 2], 'lin', [0, 1]), ...
%!              [0, Inf; 0, Inf; Inf, 0])
%!assert_equal (__nomdist__ ([1, 1; 2, 1; 2, 2], [3, 1; 3, 2], ...
%!                           [1, 1; 2, 1; 2, 2], 'iof', [1, 1]), NaN (3, 2))
%!assert_equal (__nomdist__ ([1; 1], 2, [1; 1], 'of', 1), [NaN; NaN])
%!assert_equal (__nomdist__ ([1; 1], 2, [1; 1], 'burnaby', 1), [NaN; NaN])
## Limits where nomclust returns NaN or a substitute
%!assert_equal (__nomdist__ ([1; 2], [], [1; 2], 'lin', 1), Inf)
%!assert_equal (__nomdist__ ([1, 1; 1, 1], [], [1, 1; 1, 1], 'lin', [1, 1]), 0)
%!assert_equal (__nomdist__ ([1, 1; 1, 1], [], [1, 1; 1, 1], 'lin1', ...
%!                           [1, 1]), 1)
%!test
%! h = - (2/3 * log2 (2/3) + 1/3 * log2 (1/3));
%! D = __nomdist__ ([1, 1; 1, 2; 1, 1], [], [1, 1; 1, 2; 1, 1], ...
%!                  'gambaryan', [1, 1]);
%! assert_equal (D, [1, 1 - h / 3, 1], -1e-14);
%!assert_equal (__nomdist__ (1, 1, 1, 'goodall3', 1), 0)
%!assert_equal (__nomdist__ (1, 1, 1, 'goodall4', 1), 1)

%!error<Invalid call to __nomdist__> __nomdist__ (1, [], 1, 'sm')
%!error<__nomdist__: X, Y and R must be real numeric matrices.> ...
%! __nomdist__ ({1}, [], 1, 'sm', 1)
%!error<__nomdist__: X, Y and R must have the same number of columns.> ...
%! __nomdist__ ([1, 1], [], 1, 'sm', [1, 1])
%!error<__nomdist__: X, Y and R must hold positive integer codes.> ...
%! __nomdist__ (1.5, [], 1, 'sm', 1)
%!error<__nomdist__: MEASURE must be a character vector.> ...
%! __nomdist__ (1, [], 1, 1, 1)
%!error<__nomdist__: W must hold one weight per variable.> ...
%! __nomdist__ (1, [], 1, 'sm', [1, 1])
%!error<__nomdist__: unsupported MEASURE 'morlini'.> ...
%! __nomdist__ (1, [], 1, 'morlini', 1)
%!error<__nomdist__: R must hold every level X holds.> ...
%! __nomdist__ (2, [], 1, 'sm', 1)
*/
