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

#include <octave/oct.h>
#include <octave/parse.h>

#include <cmath>
#include <string>
#include <vector>

#include "lbfgs.h"

using namespace std;

// Objective adaptor for an Octave function handle.  One call is one
// interpreter round trip returning [f, g], which is why the compiled learners
// include lbfgs.h directly instead of coming through here.
class octave_objective
{
public:

  octave_objective (const octave_value& fcn, octave_idx_type n)
    : m_fcn (fcn), m_n (n)
  { }

  double operator () (const vector<double>& x, vector<double>& g)
  {
    octave_quit ();

    ColumnVector xv (m_n);
    for (octave_idx_type i = 0; i < m_n; i++)
    {
      xv(i) = x[i];
    }

    octave_value_list in (1);
    in(0) = xv;
    octave_value_list out = octave::feval (m_fcn, in, 2);

    if (out.length () < 2)
    {
      error ("__lbfgs__: FCN must return both a value and a gradient.");
    }
    if (! out(0).isnumeric () || ! out(0).is_scalar_type ()
        || out(0).iscomplex ())
    {
      error ("__lbfgs__: FCN must return a real scalar as its first output.");
    }
    if (! out(1).isnumeric () || out(1).iscomplex ()
        || out(1).numel () != m_n)
    {
      error ("__lbfgs__: FCN must return a real gradient with as many "
             "elements as X0.");
    }

    NDArray gv = out(1).array_value ();
    for (octave_idx_type i = 0; i < m_n; i++)
    {
      g[i] = gv(i);
    }

    return out(0).double_value ();
  }

private:

  octave_value m_fcn;
  octave_idx_type m_n;
};

static bool
is_real_scalar (const octave_value& v)
{
  if (! v.isnumeric () || ! v.is_scalar_type () || v.iscomplex ())
  {
    return false;
  }
  double d = v.double_value ();
  return (d == d);
}

static bool
is_whole (double d)
{
  return (d < HUGE_VAL && d > -HUGE_VAL && d == floor (d));
}

DEFUN_DLD (__lbfgs__, args, nargout,
           "-*- texinfo -*-\n\
@deftypefn  {statistics} {@var{x} =} __lbfgs__ (@var{fcn}, @var{x0})\n\
@deftypefnx {statistics} {@var{x} =} __lbfgs__ (@var{fcn}, @var{x0}, @var{options})\n\
@deftypefnx {statistics} {[@var{x}, @var{info}] =} __lbfgs__ (@dots{})\n\
\n\
Minimize @var{fcn} by limited-memory BFGS with a strong Wolfe line search.\n\
\n\
@var{fcn} is a function handle called as @code{[f, g] = fcn (x)}, returning\n\
the objective value at @var{x} and its gradient.  @var{x0} is the real\n\
numeric vector the iteration starts from.  @var{x} is returned as a column\n\
vector whatever the orientation of @var{x0}.\n\
\n\
@var{options} is a scalar struct.  Every field is optional and the defaults\n\
are MATLAB's, as measured on R2024a.\n\
\n\
@multitable @columnfractions 0.25 0.15 0.60\n\
@headitem Field @tab Default @tab Description\n\
\n\
@item @qcode{IterationLimit} @tab @qcode{1000} @tab Largest number of\n\
iterations taken.\n\
\n\
@item @qcode{GradientTolerance} @tab @qcode{1e-6} @tab Stop once the\n\
gradient's infinity norm falls to or below this value.\n\
\n\
@item @qcode{StepTolerance} @tab @qcode{1e-6} @tab Stop once the step's\n\
infinity norm falls to or below this value.\n\
\n\
@item @qcode{LossTolerance} @tab @qcode{1e-6} @tab Stop once the objective\n\
VALUE falls to or below this value.  This is not a test on the change in the\n\
objective, which is what the name suggests; it matches MATLAB.  Pass\n\
@code{-Inf} to disable it, which any objective that can go negative, or whose\n\
minimum is zero, will want.\n\
\n\
@item @qcode{BetaTolerance} @tab @qcode{0} @tab Stop once the RELATIVE change\n\
in the iterate falls to or below this value, measured as the two-norm of the\n\
step over the two-norm of the point it landed on, with the denominator\n\
unguarded.  A value of @qcode{0} does not test against zero: it switches the\n\
test off and leaves @qcode{RelativeChangeInBeta} at @code{NaN}, as MATLAB\n\
does.  Note that this disables with @qcode{0} where @qcode{LossTolerance}\n\
disables with @code{-Inf}; each follows MATLAB for its own quantity.\n\
\n\
@item @qcode{HistorySize} @tab @qcode{10} @tab Number of curvature pairs the\n\
inverse Hessian approximation is built from.\n\
\n\
@item @qcode{InitialStepSize} @tab @qcode{[]} @tab Step length tried first on\n\
the opening iteration.  Scaled from the gradient when left empty.\n\
@end multitable\n\
\n\
When more than one tolerance is met at the same iteration the reported\n\
criterion follows MATLAB's order: gradient, then the relative change in the\n\
iterate, then step, then loss.\n\
\n\
@var{info} is a struct with fields @qcode{Iterations}, @qcode{FuncCount},\n\
@qcode{Fval}, @qcode{Gradient}, @qcode{Step}, @qcode{RelativeChangeInBeta},\n\
@qcode{Criterion} and @qcode{History}, the last being a struct of\n\
per-iteration column vectors of the same names.\n\
\n\
@qcode{Criterion} is one of @qcode{'gradient'}, @qcode{'beta'},\n\
@qcode{'step'}, @qcode{'loss'}, @qcode{'iteration'} or\n\
@qcode{'linesearch'}.  It is a token and not a sentence on purpose: MATLAB\n\
words the same criterion differently per function, so each caller maps it to\n\
the wording of the function it mirrors.\n\
\n\
This is an internal function and is not meant to be called directly.\n\
@end deftypefn")
{
  octave_value_list retval;

  if (args.length () < 2 || args.length () > 3)
  {
    print_usage ();
  }

  if (! args(0).is_function_handle ())
  {
    error ("__lbfgs__: FCN must be a function handle.");
  }

  if (! args(1).isnumeric () || args(1).iscomplex () || args(1).isempty ()
      || (args(1).rows () != 1 && args(1).columns () != 1))
  {
    error ("__lbfgs__: X0 must be a real numeric vector.");
  }

  lbfgs::options opt;

  if (args.length () == 3)
  {
    if (! args(2).isstruct () || args(2).numel () != 1)
    {
      error ("__lbfgs__: OPTIONS must be a scalar struct.");
    }

    octave_scalar_map map = args(2).scalar_map_value ();
    string_vector fields = map.fieldnames ();

    for (octave_idx_type i = 0; i < fields.numel (); i++)
    {
      std::string name = fields(i);
      octave_value val = map.contents (name);

      if (name == "IterationLimit")
      {
        if (! is_real_scalar (val) || val.double_value () < 0
            || ! is_whole (val.double_value ()))
        {
          error ("__lbfgs__: 'IterationLimit' must be a nonnegative "
                 "integer scalar.");
        }
        opt.iteration_limit = val.int_value ();
      }
      else if (name == "GradientTolerance" || name == "StepTolerance"
               || name == "BetaTolerance")
      {
        if (! is_real_scalar (val) || val.double_value () < 0)
        {
          error ("__lbfgs__: '%s' must be a nonnegative scalar.",
                 name.c_str ());
        }
        if (name == "GradientTolerance")
        {
          opt.gradient_tolerance = val.double_value ();
        }
        else if (name == "StepTolerance")
        {
          opt.step_tolerance = val.double_value ();
        }
        else
        {
          opt.beta_tolerance = val.double_value ();
        }
      }
      else if (name == "LossTolerance")
      {
        // This one is a value rather than a norm, so it is allowed to go
        // negative.  An objective that can take negative values needs -Inf
        // here, or the test fires on the first step that crosses zero.
        if (! is_real_scalar (val))
        {
          error ("__lbfgs__: 'LossTolerance' must be a real scalar.");
        }
        opt.loss_tolerance = val.double_value ();
      }
      else if (name == "HistorySize")
      {
        if (! is_real_scalar (val) || val.double_value () < 1
            || ! is_whole (val.double_value ()))
        {
          error ("__lbfgs__: 'HistorySize' must be a positive integer "
                 "scalar.");
        }
        opt.history_size = val.int_value ();
      }
      else if (name == "InitialStepSize")
      {
        // An empty value is MATLAB's way of asking for the automatic step.
        if (! val.isempty ())
        {
          if (! is_real_scalar (val) || val.double_value () <= 0)
          {
            error ("__lbfgs__: 'InitialStepSize' must be a positive "
                   "scalar.");
          }
          opt.initial_step_size = val.double_value ();
        }
      }
      else
      {
        error ("__lbfgs__: '%s' is not a valid option.", name.c_str ());
      }
    }
  }

  octave_idx_type n = args(1).numel ();
  NDArray x0 = args(1).array_value ();
  vector<double> x (n);
  for (octave_idx_type i = 0; i < n; i++)
  {
    x[i] = x0(i);
  }

  octave_objective fun (args(0), n);
  lbfgs::result res = lbfgs::minimize (fun, x, opt);

  ColumnVector xout (n);
  for (octave_idx_type i = 0; i < n; i++)
  {
    xout(i) = x[i];
  }
  retval(0) = xout;

  if (nargout > 1)
  {
    octave_idx_type nh = res.history.size ();
    ColumnVector h_iter (nh), h_fval (nh), h_grad (nh), h_step (nh),
                 h_beta (nh);
    for (octave_idx_type i = 0; i < nh; i++)
    {
      h_iter(i) = i + 1;
      h_fval(i) = res.history[i].fval;
      h_grad(i) = res.history[i].gradient;
      h_step(i) = res.history[i].step;
      h_beta(i) = res.history[i].rel_beta;
    }

    octave_scalar_map history;
    history.assign ("Iteration", h_iter);
    history.assign ("Fval", h_fval);
    history.assign ("Gradient", h_grad);
    history.assign ("Step", h_step);
    history.assign ("RelativeChangeInBeta", h_beta);

    octave_scalar_map info;
    info.assign ("Iterations", static_cast<double> (res.iterations));
    info.assign ("FuncCount", static_cast<double> (res.funcount));
    info.assign ("Fval", res.fval);
    info.assign ("Gradient", res.gradient);
    info.assign ("Step", res.step);
    info.assign ("RelativeChangeInBeta", res.rel_beta);
    info.assign ("Criterion", lbfgs::criterion_token (res.crit));
    info.assign ("History", history);

    retval(1) = info;
  }

  return retval;
}

/*
%!function [f, g] = __rosen__ (x)
%!  t = x(2) - x(1)^2;
%!  f = 100 * t^2 + (1 - x(1))^2;
%!  g = zeros (2, 1);
%!  g(1) = -400 * x(1) * t - 2 * (1 - x(1));
%!  g(2) = 200 * t;
%!endfunction

%!function [f, g] = __exrosen__ (x)
%!  n = numel (x);
%!  o = (1:2:n-1)';
%!  e = (2:2:n)';
%!  t1 = x(e) - x(o).^2;
%!  t2 = 1 - x(o);
%!  f = sum (100 * t1.^2 + t2.^2);
%!  g = zeros (n, 1);
%!  g(o) = -400 * x(o) .* t1 - 2 * t2;
%!  g(e) = 200 * t1;
%!endfunction

%!function [f, g] = __beale__ (x)
%!  y = [1.5; 2.25; 2.625];
%!  k = (1:3)';
%!  t = 1 - x(2).^k;
%!  r = y - x(1) * t;
%!  f = sum (r.^2);
%!  g = zeros (2, 1);
%!  g(1) = -2 * sum (r .* t);
%!  g(2) = 2 * x(1) * sum (r .* k .* x(2).^(k - 1));
%!endfunction

%!function [f, g] = __wood__ (x)
%!  t1 = 100 * (x(2) - x(1)^2)^2 + (1 - x(1))^2;
%!  t2 = 90 * (x(4) - x(3)^2)^2 + (1 - x(3))^2;
%!  t3 = 10.1 * ((1 - x(2))^2 + (1 - x(4))^2);
%!  t4 = 19.8 * (1 - x(2)) * (1 - x(4));
%!  f = t1 + t2 + t3 + t4;
%!  g = zeros (4, 1);
%!  g(1) = -400 * x(1) * (x(2) - x(1)^2) - 2 * (1 - x(1));
%!  g(2) = 200 * (x(2) - x(1)^2) - 20.2 * (1 - x(2)) - 19.8 * (1 - x(4));
%!  g(3) = -360 * x(3) * (x(4) - x(3)^2) - 2 * (1 - x(3));
%!  g(4) = 180 * (x(4) - x(3)^2) - 20.2 * (1 - x(4)) - 19.8 * (1 - x(2));
%!endfunction

%!function [f, g] = __powellsq__ (x)
%!  a = x(1) + 10 * x(2);
%!  b = x(3) - x(4);
%!  c = x(2) - 2 * x(3);
%!  d = x(1) - x(4);
%!  f = a^2 + 5 * b^2 + c^4 + 10 * d^4;
%!  g = zeros (4, 1);
%!  g(1) = 2 * a + 40 * d^3;
%!  g(2) = 20 * a + 4 * c^3;
%!  g(3) = 10 * b - 8 * c^3;
%!  g(4) = -10 * b - 40 * d^3;
%!endfunction

%!function [f, g] = __quadtri__ (x)
%!  n = numel (x);
%!  A = 2 * eye (n) - diag (ones (n - 1, 1), 1) - diag (ones (n - 1, 1), -1);
%!  xs = ((1:n)' / n).^2;
%!  r = x - xs;
%!  f = 0.5 * r' * A * r;
%!  g = A * r;
%!endfunction

%!function [f, g] = __badgrad__ (x)
%!  f = sum (x.^2);
%!  g = 1;
%!endfunction

## Every Moré, Garbow and Hillstrom problem below has a zero minimum, so the
## loss test is switched off and convergence is judged by the gradient.  Left
## on, it stops the iteration first; that is the point of the next test.

## Problem 1, the standard start.
%!test
%! opt = struct ("LossTolerance", -Inf);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! crit = info.Criterion;
%! assert_equal (x, [1; 1], 1e-5);
%! assert_equal (info.Fval < 1e-8, true);
%! assert_equal (crit, "gradient");

## LossTolerance tests the objective VALUE, not the change in it, so on a
## problem whose minimum is zero the default 1e-6 pre-empts the gradient test
## while the gradient is still six orders of magnitude above its own
## tolerance.  Measured on MATLAB R2024a, where a LossTolerance of 0.5 against
## a first-iteration loss of 0.3551 stopped fitcnet at iteration 1.
%!test
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1]);
%! crit = info.Criterion;
%! assert_equal (crit, "loss");
%! assert_equal (info.Fval <= 1e-6, true);
%! assert_equal (info.Gradient > 1e-6, true);

## Problem 21 at n = 20, where the parameters outnumber the stored curvature
## pairs and the limited-memory recursion is what is actually being tested.
%!test
%! n = 20;
%! x0 = zeros (n, 1);
%! x0(1:2:n-1) = -1.2;
%! x0(2:2:n) = 1;
%! opt = struct ("LossTolerance", -Inf);
%! [x, info] = __lbfgs__ (@__exrosen__, x0, opt);
%! assert_equal (x, ones (n, 1), 1e-4);

## Problem 5.
%!test
%! opt = struct ("LossTolerance", -Inf);
%! [x, info] = __lbfgs__ (@__beale__, [1; 1], opt);
%! assert_equal (x, [3; 0.5], 1e-4);

## Problem 14.
%!test
%! opt = struct ("LossTolerance", -Inf);
%! [x, info] = __lbfgs__ (@__wood__, [-3; -1; -3; -1], opt);
%! assert_equal (x, ones (4, 1), 1e-4);

## Problem 13.  The Hessian is singular at the solution, so convergence is
## linear and the objective is the only thing worth asserting on.
%!test
%! opt = struct ("LossTolerance", -Inf);
%! [x, info] = __lbfgs__ (@__powellsq__, [3; -1; 0; 1], opt);
%! assert_equal (info.Fval < 1e-8, true);

## A quadratic in fifty variables, well past the ten curvature pairs the
## history holds, against its known minimiser.  StepTolerance is switched off
## along with the loss test: the steps go below 1e-6 while the gradient is
## still far from its own tolerance, and stopping there is correct but is not
## what this test is measuring.
%!test
%! n = 50;
%! xs = ((1:n)' / n).^2;
%! opt = struct ("GradientTolerance", 1e-10, "LossTolerance", -Inf, ...
%!               "StepTolerance", 0);
%! [x, info] = __lbfgs__ (@__quadtri__, zeros (n, 1), opt);
%! crit = info.Criterion;
%! assert_equal (x, xs, 1e-6);
%! assert_equal (crit, "gradient");

## A single stored pair still converges, it just takes longer.
%!test
%! opt = struct ("HistorySize", 1, "LossTolerance", -Inf);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! assert_equal (x, [1; 1], 1e-4);

## X0 may be a row, X never is.
%!test
%! opt = struct ("LossTolerance", -Inf);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2, 1], opt);
%! assert_equal (x, [1; 1], 1e-5);

## An empty InitialStepSize asks for the automatic one, as it does in MATLAB.
%!test
%! opt = struct ("InitialStepSize", 1e-3, "LossTolerance", -Inf);
%! [x1, i1] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! opt = struct ("InitialStepSize", [], "LossTolerance", -Inf);
%! [x2, i2] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! assert_equal (x1, [1; 1], 1e-5);
%! assert_equal (x2, [1; 1], 1e-5);

## A loss tolerance well above the starting value stops at once.
%!test
%! opt = struct ("LossTolerance", 1);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! crit = info.Criterion;
%! assert_equal (crit, "loss");
%! assert_equal (info.Fval <= 1, true);

## Reported criterion when several are met at once: gradient outranks step.
%!test
%! opt = struct ("GradientTolerance", 1e3, "StepTolerance", 1e3);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! crit = info.Criterion;
%! assert_equal (crit, "gradient");
%! assert_equal (info.Iterations, 1);

## And step outranks loss.
%!test
%! opt = struct ("StepTolerance", 1e3, "LossTolerance", 1e3);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! crit = info.Criterion;
%! assert_equal (crit, "step");
%! assert_equal (info.Iterations, 1);

## The history holds one row per iteration taken, and the line search costs
## objective calls over and above them.
%!test
%! opt = struct ("IterationLimit", 4);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! crit = info.Criterion;
%! assert_equal (crit, "iteration");
%! assert_equal (info.Iterations, 4);
%! assert_equal (info.History.Iteration, (1:4)');
%! assert_equal (numel (info.History.Fval), 4);
%! assert_equal (info.FuncCount >= info.Iterations, true);

## A zero limit evaluates the start point and stops there.
%!test
%! opt = struct ("IterationLimit", 0);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! assert_equal (x, [-1.2; 1]);
%! assert_equal (info.Iterations, 0);
%! assert_equal (isempty (info.History.Fval), true);

## BetaTolerance is off by default, and off means the quantity is not
## computed rather than tested against zero: MATLAB reports NaN for it, and
## so does this.
%!test
%! opt = struct ("LossTolerance", -Inf, "IterationLimit", 5);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! assert_equal (info.RelativeChangeInBeta, NaN);
%! assert_equal (info.History.RelativeChangeInBeta, NaN (5, 1));

## The relative change is the two-norm of the step over the two-norm of the
## point it landed on, unguarded.  From a zero start the first iteration is
## therefore exactly 1, which is what R2024a's fitclinear reports and what
## rules out every max (1, .) guarded candidate.
%!test
%! n = 50;
%! opt = struct ("LossTolerance", -Inf, "BetaTolerance", 1e-12, ...
%!               "IterationLimit", 4);
%! [x, info] = __lbfgs__ (@__quadtri__, zeros (n, 1), opt);
%! assert_equal (info.History.RelativeChangeInBeta(1), 1);

## It stops on the relative change, and reports which test did it.
%!test
%! opt = struct ("LossTolerance", -Inf, "BetaTolerance", 1e3);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! assert_equal (info.Criterion, "beta");
%! assert_equal (info.Iterations, 1);

## Gradient outranks it, and it outranks step, which is MATLAB's order.
%!test
%! opt = struct ("GradientTolerance", 1e3, "BetaTolerance", 1e3, ...
%!               "StepTolerance", 1e3);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! assert_equal (info.Criterion, "gradient");
%!test
%! opt = struct ("BetaTolerance", 1e3, "StepTolerance", 1e3, ...
%!               "LossTolerance", 1e3);
%! [x, info] = __lbfgs__ (@__rosen__, [-1.2; 1], opt);
%! assert_equal (info.Criterion, "beta");

%!error <Invalid call to __lbfgs__.> __lbfgs__ (@__rosen__)
%!error <__lbfgs__: FCN must be a function handle.> ...
%! __lbfgs__ (1, [1; 1])
%!error <__lbfgs__: X0 must be a real numeric vector.> ...
%! __lbfgs__ (@__rosen__, {1})
%!error <__lbfgs__: X0 must be a real numeric vector.> ...
%! __lbfgs__ (@__rosen__, ones (2, 2))
%!error <__lbfgs__: X0 must be a real numeric vector.> ...
%! __lbfgs__ (@__rosen__, [])
%!error <__lbfgs__: X0 must be a real numeric vector.> ...
%! __lbfgs__ (@__rosen__, complex ([1; 1]))
%!error <__lbfgs__: OPTIONS must be a scalar struct.> ...
%! __lbfgs__ (@__rosen__, [1; 1], 5)
%!error <__lbfgs__: 'IterationLimit' must be a nonnegative integer scalar.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("IterationLimit", -1))
%!error <__lbfgs__: 'IterationLimit' must be a nonnegative integer scalar.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("IterationLimit", 2.5))
%!error <__lbfgs__: 'GradientTolerance' must be a nonnegative scalar.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("GradientTolerance", -1))
%!error <__lbfgs__: 'StepTolerance' must be a nonnegative scalar.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("StepTolerance", -1))
%!error <__lbfgs__: 'LossTolerance' must be a real scalar.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("LossTolerance", NaN))
%!error <__lbfgs__: 'BetaTolerance' must be a nonnegative scalar.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("BetaTolerance", -1))
%!error <__lbfgs__: 'HistorySize' must be a positive integer scalar.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("HistorySize", 0))
%!error <__lbfgs__: 'InitialStepSize' must be a positive scalar.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("InitialStepSize", -1))
%!error <__lbfgs__: 'Bogus' is not a valid option.> ...
%! __lbfgs__ (@__rosen__, [1; 1], struct ("Bogus", 1))
%!error <__lbfgs__: FCN must return a real gradient with as many elements as X0.> ...
%! __lbfgs__ (@__badgrad__, [1; 1])
*/
