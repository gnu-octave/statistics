/*
Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>

This file is part of the statistics package for GNU Octave.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with
this program; if not, see <http://www.gnu.org/licenses/>.
*/

// Limited-memory BFGS with a strong Wolfe line search.
//
// The engine is deliberately free of any Octave type so that it can serve two
// callers with opposite needs: a compiled learner whose objective is C++ and
// which must never cross the interpreter boundary, and the __lbfgs__ oct-file,
// whose objective is an Octave function handle.  The objective is therefore a
// template parameter, and anything modelling
//
//   double operator () (const std::vector<double> &x, std::vector<double> &g)
//
// which returns the value at x and writes the gradient into g, will do.
//
// Stopping behaviour follows MATLAB's, as measured on R2024a with the
// Statistics and Machine Learning Toolbox 24.1.  Three findings are baked in
// here because they are not what the documentation's wording suggests, and
// a fourth that it does not state at all:
//
//   1. All three tolerances are ABSOLUTE comparisons against the scalar the
//      iteration records, not relative to that scalar's starting value.
//      MATLAB reports "Relative gradient tolerance reached." but stops at the
//      first iteration whose recorded gradient falls to or below the
//      tolerance; scaling the objective by 1e6 moved the recorded gradient up
//      rather than leaving it fixed, so no division by the objective value is
//      involved either.
//   2. LossTolerance tests the objective VALUE, not the change in it.  Given a
//      tolerance of 0.5 against a first-iteration loss of 0.3551, MATLAB
//      stopped at iteration 1.
//   3. When several tolerances are met at once the reported criterion follows
//      a fixed order: gradient, then step, then loss.
//   4. The gradient and step are reduced by their infinity norm, which is
//      what MATLAB reports.  Letting R2024a's fitrnet converge and computing
//      both norms at the weights it returned, its reported gradient matched
//      the infinity norm to twelve significant digits and the two-norm not at
//      all, the latter being 9.5 per cent larger on that fit.

#ifndef OCTAVE_STATISTICS_LBFGS_H
#define OCTAVE_STATISTICS_LBFGS_H

#include <cmath>
#include <cstddef>
#include <vector>

using namespace std;

namespace lbfgs
{

  // Why the iteration stopped.  Callers map these onto MATLAB's own strings.
  enum criterion
  {
    CRIT_GRADIENT = 0,
    CRIT_STEP,
    CRIT_LOSS,
    CRIT_ITERATION_LIMIT,
    CRIT_LINE_SEARCH
  };

  struct options
  {
    int iteration_limit;
    double gradient_tolerance;
    double loss_tolerance;
    double step_tolerance;
    int history_size;
    double initial_step_size;   // non-positive selects the automatic value
    int max_line_search;
    double wolfe_c1;
    double wolfe_c2;

    // Defaults are MATLAB's, read off fitcnet's ModelParameters on R2024a.
    // HistorySize is not exposed by fitcnet; 10 is what fmincon and the deep
    // learning solver both use, and 15 is what fitclinear uses.
    options ()
      : iteration_limit (1000), gradient_tolerance (1e-6),
        loss_tolerance (1e-6), step_tolerance (1e-6), history_size (10),
        initial_step_size (0.0), max_line_search (25), wolfe_c1 (1e-4),
        wolfe_c2 (0.9)
    { }
  };

  struct record
  {
    double fval;
    double gradient;
    double step;
  };

  struct result
  {
    int crit;
    int iterations;
    int funcount;
    double fval;
    double gradient;
    double step;
    vector<record> history;
  };

  // MATLAB's own wording, kept verbatim so a caller rebuilding fitcnet's
  // ConvergenceInfo can pass the string straight through.  The gradient test
  // it names is absolute rather than relative; see the note above.
  inline const char * criterion_message (int crit)
  {
    switch (crit)
    {
      case CRIT_GRADIENT:
        return "Relative gradient tolerance reached.";
      case CRIT_STEP:
        return "Step size tolerance reached.";
      case CRIT_LOSS:
        return "Loss tolerance reached.";
      case CRIT_ITERATION_LIMIT:
        return "Iteration limit reached.";
      default:
        return "Line search could not improve the objective.";
    }
  }

  inline double dot (const vector<double> &a, const vector<double> &b)
  {
    double s = 0.0;
    for (size_t i = 0; i < a.size (); i++)
    {
      s += a[i] * b[i];
    }
    return s;
  }

  inline double inf_norm (const vector<double> &a)
  {
    double m = 0.0;
    for (size_t i = 0; i < a.size (); i++)
    {
      double v = fabs (a[i]);
      if (v > m)
      {
        m = v;
      }
    }
    return m;
  }

  inline bool is_finite (double v)
  {
    return (v == v && v < HUGE_VAL && v > -HUGE_VAL);
  }

  // Safeguarded quadratic interpolation between the two bracket ends, falling
  // back to bisection whenever the interpolant is not usable.  Keeping the
  // trial away from the ends by a tenth of the bracket is what stops the zoom
  // from stalling on a flat interpolant.
  inline double interpolate (double a_lo, double f_lo, double d_lo,
                             double a_hi, double f_hi)
  {
    double lo = (a_lo < a_hi ? a_lo : a_hi);
    double hi = (a_lo < a_hi ? a_hi : a_lo);
    double edge = 0.1 * (hi - lo);
    double bisect = 0.5 * (a_lo + a_hi);

    double width = a_hi - a_lo;
    double denom = 2.0 * (f_hi - f_lo - d_lo * width);
    if (denom == 0.0 || ! is_finite (denom))
    {
      return bisect;
    }

    double a = a_lo - d_lo * width * width / denom;
    if (! is_finite (a) || a < lo + edge || a > hi - edge)
    {
      return bisect;
    }
    return a;
  }

  // Strong Wolfe line search, Nocedal and Wright algorithms 3.5 and 3.6.  On
  // success the accepted point, its gradient and its value are left in x, g
  // and f, so an accepted search costs the caller no further objective call.
  template <class Objective>
  bool line_search (Objective &fun, const vector<double> &x0,
                    const vector<double> &d, double f0, double dphi0,
                    double alpha0, const options &opt, vector<double> &x,
                    vector<double> &g, double &f, double &alpha, int &nfev)
  {
    const size_t n = x0.size ();
    const double c1 = opt.wolfe_c1;
    const double c2 = opt.wolfe_c2;
    double dphi = 0.0;

    // phi (a) together with phi' (a), evaluated in place: the trial point is
    // left in x, its gradient in g and its value in f, so an accepted step
    // needs no repeat call.
    auto evaluate = [&] (double a) -> void
    {
      for (size_t i = 0; i < n; i++)
      {
        x[i] = x0[i] + a * d[i];
      }
      f = fun (x, g);
      dphi = dot (g, d);
      nfev++;
    };

    double a_lo = 0.0, f_lo = f0, d_lo = dphi0;
    double a_hi = 0.0, f_hi = 0.0;
    bool bracketed = false;

    double a_prev = 0.0, f_prev = f0, d_prev = dphi0;
    double a = alpha0;

    // Bracketing.  The first trial that breaks sufficient decrease, or that
    // rises above its predecessor, or that turns the slope upwards, encloses
    // an acceptable step together with the trial before it.
    for (int i = 1; i <= opt.max_line_search; i++)
    {
      evaluate (a);

      if (! is_finite (f) || f > f0 + c1 * a * dphi0 || (i > 1 && f >= f_prev))
      {
        a_lo = a_prev; f_lo = f_prev; d_lo = d_prev;
        a_hi = a; f_hi = f;
        bracketed = true;
        break;
      }

      if (fabs (dphi) <= -c2 * dphi0)
      {
        alpha = a;
        return true;
      }

      if (dphi >= 0.0)
      {
        a_lo = a; f_lo = f; d_lo = dphi;
        a_hi = a_prev; f_hi = f_prev;
        bracketed = true;
        break;
      }

      a_prev = a; f_prev = f; d_prev = dphi;
      a *= 2.0;
    }

    if (! bracketed)
    {
      // Every trial improved and none met the curvature condition, which
      // means the last one is still a descent step worth taking.
      if (a_prev > 0.0)
      {
        evaluate (a_prev);
        alpha = a_prev;
        return true;
      }
      return false;
    }

    // Zoom.  The bracket shrinks until a trial satisfies both Wolfe
    // conditions, and a_lo always holds the best point seen so far.
    for (int j = 1; j <= opt.max_line_search; j++)
    {
      a = interpolate (a_lo, f_lo, d_lo, a_hi, f_hi);
      evaluate (a);

      if (! is_finite (f) || f > f0 + c1 * a * dphi0 || f >= f_lo)
      {
        a_hi = a; f_hi = f;
      }
      else
      {
        if (fabs (dphi) <= -c2 * dphi0)
        {
          alpha = a;
          return true;
        }
        if (dphi * (a_hi - a_lo) >= 0.0)
        {
          a_hi = a_lo; f_hi = f_lo;
        }
        a_lo = a; f_lo = f; d_lo = dphi;
      }

      if (fabs (a_hi - a_lo) <= 1e-16 * (1.0 + fabs (a_lo)))
      {
        break;
      }
    }

    // The budget ran out.  Accepting the best bracket end is better than
    // failing outright: it still decreases the objective, and the curvature
    // safeguard on the (s, y) pair will discard the update if it has to.
    if (a_lo > 0.0)
    {
      evaluate (a_lo);
      alpha = a_lo;
      return true;
    }
    return false;
  }

  // Minimise fun starting from x, which is overwritten with the final point.
  template <class Objective>
  result minimize (Objective &fun, vector<double> &x, const options &opt)
  {
    const size_t n = x.size ();
    const int m = (opt.history_size > 0 ? opt.history_size : 1);

    vector<double> g (n, 0.0), d (n, 0.0), q (n, 0.0);
    vector<double> x_new (n, 0.0), g_new (n, 0.0);
    vector<vector<double> > S (m, vector<double> (n, 0.0));
    vector<vector<double> > Y (m, vector<double> (n, 0.0));
    vector<double> rho (m, 0.0), coef (m, 0.0);
    int stored = 0, newest = -1;
    double gamma = 1.0;

    result out;
    out.crit = CRIT_ITERATION_LIMIT;
    out.iterations = 0;
    out.funcount = 1;
    out.fval = fun (x, g);
    out.gradient = inf_norm (g);
    out.step = 0.0;

    for (int k = 1; k <= opt.iteration_limit; k++)
    {
      // Two-loop recursion, d = -H * g.
      q = g;
      for (int j = 0; j < stored; j++)
      {
        int idx = (newest - j + m) % m;
        coef[idx] = rho[idx] * dot (S[idx], q);
        for (size_t i = 0; i < n; i++)
        {
          q[i] -= coef[idx] * Y[idx][i];
        }
      }
      for (size_t i = 0; i < n; i++)
      {
        q[i] *= gamma;
      }
      for (int j = stored - 1; j >= 0; j--)
      {
        int idx = (newest - j + m) % m;
        double beta = rho[idx] * dot (Y[idx], q);
        for (size_t i = 0; i < n; i++)
        {
          q[i] += (coef[idx] - beta) * S[idx][i];
        }
      }
      for (size_t i = 0; i < n; i++)
      {
        d[i] = -q[i];
      }

      // Rounding can leave the recursion with a direction that is not one of
      // descent; dropping the history and stepping down the gradient is the
      // standard recovery.
      double dphi0 = dot (g, d);
      if (! (dphi0 < 0.0))
      {
        for (size_t i = 0; i < n; i++)
        {
          d[i] = -g[i];
        }
        dphi0 = -dot (g, g);
        stored = 0;
        newest = -1;
        gamma = 1.0;
      }

      // A zero slope down the steepest descent direction means the gradient
      // itself is zero, so there is nothing left to search along.
      if (dphi0 == 0.0)
      {
        out.crit = CRIT_GRADIENT;
        break;
      }

      // A quasi-Newton step is unit-scaled once curvature information exists.
      // The first one has none, so it is scaled by the gradient instead.
      double alpha0 = 1.0;
      if (stored == 0)
      {
        if (opt.initial_step_size > 0.0)
        {
          alpha0 = opt.initial_step_size;
        }
        else
        {
          double gn = inf_norm (g);
          alpha0 = (gn > 1.0 ? 1.0 / gn : 1.0);
        }
      }

      // f_new is kept separate from out.fval: a failed search leaves its own
      // scratch value behind, and out.fval must go on describing the last
      // point actually accepted.
      double f_new = out.fval;
      double alpha = 0.0;
      bool ok = line_search (fun, x, d, out.fval, dphi0, alpha0, opt, x_new,
                             g_new, f_new, alpha, out.funcount);

      // A failed search usually means the stored curvature has gone stale
      // rather than that the point is stationary, so the memory is dropped
      // and the step retried down the gradient before giving up.  Only a
      // second failure, with nothing left to blame, ends the iteration.
      if (! ok && stored > 0)
      {
        stored = 0;
        newest = -1;
        gamma = 1.0;
        for (size_t i = 0; i < n; i++)
        {
          d[i] = -g[i];
        }
        dphi0 = -dot (g, g);
        double gn = inf_norm (g);
        alpha0 = (gn > 1.0 ? 1.0 / gn : 1.0);
        f_new = out.fval;
        ok = line_search (fun, x, d, out.fval, dphi0, alpha0, opt, x_new,
                          g_new, f_new, alpha, out.funcount);
      }

      if (! ok)
      {
        out.crit = CRIT_LINE_SEARCH;
        break;
      }

      // The curvature pair goes straight into the slot it would occupy, so
      // that the step norm can be read off it before the pair is judged.  A
      // rejected pair simply leaves the slot to be overwritten next time.
      int idx = (newest + 1) % m;
      double sy = 0.0, yy = 0.0;
      for (size_t i = 0; i < n; i++)
      {
        S[idx][i] = x_new[i] - x[i];
        Y[idx][i] = g_new[i] - g[i];
        sy += S[idx][i] * Y[idx][i];
        yy += Y[idx][i] * Y[idx][i];
      }
      double stepnorm = inf_norm (S[idx]);

      // Keeping only pairs with positive curvature is what holds the implicit
      // inverse Hessian positive definite.
      if (yy > 0.0 && sy > 1e-10 * yy)
      {
        rho[idx] = 1.0 / sy;
        gamma = sy / yy;
        newest = idx;
        if (stored < m)
        {
          stored++;
        }
      }

      x.swap (x_new);
      g.swap (g_new);
      out.fval = f_new;
      out.gradient = inf_norm (g);
      out.step = stepnorm;
      out.iterations = k;

      record r;
      r.fval = out.fval;
      r.gradient = out.gradient;
      r.step = out.step;
      out.history.push_back (r);

      // MATLAB's order, verified on R2024a: gradient, then step, then loss.
      if (out.gradient <= opt.gradient_tolerance)
      {
        out.crit = CRIT_GRADIENT;
        break;
      }
      if (out.step <= opt.step_tolerance)
      {
        out.crit = CRIT_STEP;
        break;
      }
      if (out.fval <= opt.loss_tolerance)
      {
        out.crit = CRIT_LOSS;
        break;
      }
    }

    return out;
  }

}

#endif
