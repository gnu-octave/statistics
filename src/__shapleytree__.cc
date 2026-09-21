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

#include <cmath>
#include <vector>
#include <octave/oct.h>
#include <octave/Cell.h>

// One demand a leaf makes: the value of predictor VAR must reach the child
// the path took.  A cut on a threshold sends a value below it left, and a cut
// on levels sends the levels it lists either way and anything else neither,
// which is what a missing value does too.
struct Demand
{
  octave_idx_type var;
  bool left;
  bool categorical;
  double cut;
  std::vector<double> levels;   // the levels of the side the path took
};

// A leaf: what it answers, and everything the way to it demands, gathered by
// predictor so that a predictor cut more than once counts as one.
struct Leaf
{
  octave_idx_type node;
  std::vector<octave_idx_type> vars;          // distinct, ascending
  std::vector<std::vector<Demand> > demands;  // one list per entry of vars
};

// A tree as this needs it: the node table, and what it answers at each node.
struct Tree
{
  Matrix children;
  ColumnVector cutvar;
  ColumnVector cutpoint;
  Cell catleft;
  Cell catright;
  Matrix leafval;
  bool hascat;
  std::vector<Leaf> leaves;
};

static bool
demand_met (const Demand& d, double v)
{
  if (d.categorical)
    {
      for (size_t i = 0; i < d.levels.size (); i++)
        if (v == d.levels[i])
          return true;
      return false;
    }
  if (octave::math::isnan (v))
    return false;
  return d.left ? (v < d.cut) : (v >= d.cut);
}

// Walk the tree once and record every leaf with the demands on the way to it.
static void
collect_leaves (Tree& T)
{
  const octave_idx_type nn = T.children.rows ();
  if (nn <= 0)
    return;

  std::vector<octave_idx_type> stack_node;
  std::vector<std::vector<Demand> > stack_path;
  stack_node.push_back (0);
  stack_path.push_back (std::vector<Demand> ());

  while (! stack_node.empty ())
    {
      const octave_idx_type n = stack_node.back ();
      std::vector<Demand> path = stack_path.back ();
      stack_node.pop_back ();
      stack_path.pop_back ();

      const octave_idx_type var =
        static_cast<octave_idx_type> (T.cutvar(n));
      if (var == 0)
        {
          Leaf lf;
          lf.node = n;
          // Gather the demands by predictor, keeping the predictors ascending
          for (size_t i = 0; i < path.size (); i++)
            {
              size_t k = 0;
              while (k < lf.vars.size () && lf.vars[k] < path[i].var)
                k++;
              if (k == lf.vars.size () || lf.vars[k] != path[i].var)
                {
                  lf.vars.insert (lf.vars.begin () + k, path[i].var);
                  lf.demands.insert (lf.demands.begin () + k,
                                     std::vector<Demand> ());
                }
              lf.demands[k].push_back (path[i]);
            }
          T.leaves.push_back (lf);
          continue;
        }

      for (int kid = 0; kid < 2; kid++)
        {
          const octave_idx_type child =
            static_cast<octave_idx_type> (T.children(n, kid)) - 1;
          if (child < 0 || child >= nn)
            continue;
          Demand d;
          d.var = var - 1;
          d.left = (kid == 0);
          d.categorical = false;
          d.cut = T.cutpoint(n);
          if (T.hascat)
            {
              const octave_value cv = (kid == 0 ? T.catleft(n) : T.catright(n));
              if (! cv.isempty ())
                {
                  d.categorical = true;
                  const NDArray lv = cv.array_value ();
                  for (octave_idx_type j = 0; j < lv.numel (); j++)
                    d.levels.push_back (lv(j));
                }
            }
          std::vector<Demand> next = path;
          next.push_back (d);
          stack_node.push_back (child);
          stack_path.push_back (next);
        }
    }
}

// The sum, over every way of adding T of the M-A-B predictors a leaf does not
// demand, of the weight the Shapley definition gives a subset of that size.
// Taken in logs, so that neither the binomial coefficient nor the factorials
// overflow however wide the model is.
static double
kernel_sum (octave_idx_type k, octave_idx_type m, octave_idx_type M,
            const std::vector<double>& lg)
{
  double tot = 0.0;
  for (octave_idx_type t = 0; t <= m; t++)
    {
      const octave_idx_type s = k + t;
      const double lt = lg[m] - lg[t] - lg[m - t]
                        + lg[s] + lg[M - s - 1] - lg[M];
      tot += std::exp (lt);
    }
  return tot;
}

// The Shapley values of one tree at one query point, added into PHI.
//
// A leaf answers for the rows meeting every demand on the way to it.  Hold
// the predictors of a subset at the query point and leave the rest as an
// observation has them, and that row reaches the leaf exactly when every
// predictor the query fails is outside the subset and every predictor the
// observation fails is inside it.  Writing A for the predictors the query
// passes and the observation fails and B for the other way about, the leaf is
// reached exactly for the subsets holding all of A and none of B, so only
// those predictors can move a value, and by how much turns on nothing but how
// many there are.
static void
tree_values (const Tree& T, const Matrix& X, const RowVector& q,
             octave_idx_type M, octave_idx_type K,
             const std::vector<double>& lg, Matrix& phi)
{
  const octave_idx_type n = X.rows ();
  std::vector<char> passX;
  std::vector<char> passZ;
  std::vector<octave_idx_type> acount;
  std::vector<char> keep;

  for (size_t li = 0; li < T.leaves.size (); li++)
    {
      const Leaf& lf = T.leaves[li];
      const size_t nf = lf.vars.size ();
      if (nf == 0)
        continue;

      // Which of the leaf's predictors the query passes
      passX.assign (nf, 1);
      for (size_t j = 0; j < nf; j++)
        for (size_t d = 0; d < lf.demands[j].size (); d++)
          if (! demand_met (lf.demands[j][d], q(lf.vars[j])))
            {
              passX[j] = 0;
              break;
            }

      // And which each observation passes
      passZ.assign (static_cast<size_t> (n) * nf, 1);
      for (size_t j = 0; j < nf; j++)
        {
          const octave_idx_type v = lf.vars[j];
          for (octave_idx_type i = 0; i < n; i++)
            {
              char ok = 1;
              for (size_t d = 0; d < lf.demands[j].size (); d++)
                if (! demand_met (lf.demands[j][d], X(i, v)))
                  {
                    ok = 0;
                    break;
                  }
              passZ[j * static_cast<size_t> (n) + i] = ok;
            }
        }

      // An observation failing a predictor the query fails too puts the leaf
      // out of reach whatever the subset is
      octave_idx_type b = 0;
      for (size_t j = 0; j < nf; j++)
        if (! passX[j])
          b++;

      keep.assign (n, 1);
      for (size_t j = 0; j < nf; j++)
        if (! passX[j])
          for (octave_idx_type i = 0; i < n; i++)
            if (! passZ[j * static_cast<size_t> (n) + i])
              keep[i] = 0;

      acount.assign (n, 0);
      octave_idx_type nkeep = 0;
      octave_idx_type amax = 0;
      for (octave_idx_type i = 0; i < n; i++)
        {
          if (! keep[i])
            continue;
          nkeep++;
          octave_idx_type a = 0;
          for (size_t j = 0; j < nf; j++)
            if (passX[j] && ! passZ[j * static_cast<size_t> (n) + i])
              a++;
          acount[i] = a;
          if (a > amax)
            amax = a;
        }
      if (nkeep == 0)
        continue;

      // The weight turns on an observation only through how many predictors
      // fall in A, so it is worked out once for each count that occurs
      std::vector<double> gA (amax + 1, 0.0), gB (amax + 1, 0.0);
      std::vector<octave_idx_type> nat (amax + 1, 0);
      for (octave_idx_type a = 0; a <= amax; a++)
        {
          const octave_idx_type mv = M - a - b;
          if (mv < 0)
            continue;
          gB[a] = kernel_sum (a, mv, M, lg);
          if (a >= 1)
            gA[a] = kernel_sum (a - 1, mv, M, lg);
        }
      for (octave_idx_type i = 0; i < n; i++)
        if (keep[i])
          nat[acount[i]]++;

      for (octave_idx_type kk = 0; kk < K; kk++)
        {
          const double vL = T.leafval(lf.node, kk);
          if (vL == 0.0)
            continue;

          if (b > 0)
            {
              double wsum = 0.0;
              for (octave_idx_type a = 0; a <= amax; a++)
                wsum += static_cast<double> (nat[a]) * gB[a];
              for (size_t j = 0; j < nf; j++)
                if (! passX[j])
                  phi(lf.vars[j], kk) -= wsum * vL;
            }

          for (size_t j = 0; j < nf; j++)
            {
              if (! passX[j])
                continue;
              double wsum = 0.0;
              for (octave_idx_type i = 0; i < n; i++)
                if (keep[i] && ! passZ[j * static_cast<size_t> (n) + i])
                  wsum += gA[acount[i]];
              phi(lf.vars[j], kk) += wsum * vL;
            }
        }
    }
}

DEFUN_DLD (__shapleytree__, args, ,
           "-*- texinfo -*-\n\
@deftypefn {statistics} {@var{phi} =} \
__shapleytree__ (@var{trees}, @var{X}, @var{Q})\n\
\n\
Shapley values of an ensemble of decision trees, taken leaf by leaf.\n\
Internal; called by @code{shapley} and not meant to be used directly.\n\
\n\
@var{trees} is a cell array of structures, one per tree, each holding\n\
@qcode{Children}, the @math{Nx2} table of child node numbers and zero on a\n\
leaf; @qcode{CutVar}, the predictor each node cuts on and zero on a leaf;\n\
@qcode{CutPoint}, the value it cuts at; @qcode{CatLeft} and @qcode{CatRight},\n\
the levels each side of a cut on a categorical predictor takes and empty\n\
where the cut is on a threshold; and @qcode{Leaf}, what the tree answers at\n\
each node, one column per class and one column for a response, already\n\
carrying whatever weight the model gives the tree.\n\
\n\
@var{X} holds the observations averaged over, one per row, and @var{Q} the\n\
query points, one per row.  Both have one column per predictor.\n\
\n\
@var{phi} is @math{MxKxnq}: one row per predictor, one column per column of\n\
@qcode{Leaf}, and one page per query point.  The values of a query point sum\n\
to the deviation of the prediction from the average prediction.\n\
\n\
@end deftypefn")
{
  if (args.length () != 3)
    print_usage ();

  if (! args(0).iscell ())
    error ("__shapleytree__: TREES must be a cell array of structures.");
  if (! args(1).isnumeric () || args(1).iscomplex () || args(1).isempty ())
    error ("__shapleytree__: X must be a real numeric matrix.");
  if (! args(2).isnumeric () || args(2).iscomplex () || args(2).isempty ())
    error ("__shapleytree__: Q must be a real numeric matrix.");

  const Cell tc = args(0).cell_value ();
  const Matrix X = args(1).matrix_value ();
  const Matrix Q = args(2).matrix_value ();
  const octave_idx_type M = X.columns ();
  const octave_idx_type n = X.rows ();
  const octave_idx_type nq = Q.rows ();

  if (Q.columns () != M)
    error ("__shapleytree__: X and Q must have the same number of columns.");
  if (tc.numel () == 0)
    error ("__shapleytree__: TREES must hold at least one tree.");

  std::vector<Tree> trees (tc.numel ());
  octave_idx_type K = -1;
  for (octave_idx_type t = 0; t < tc.numel (); t++)
    {
      if (! tc(t).isstruct ())
        error ("__shapleytree__: TREES must be a cell array of structures.");
      const octave_scalar_map sm = tc(t).scalar_map_value ();
      Tree& T = trees[t];
      T.children = sm.contents ("Children").matrix_value ();
      T.cutvar = sm.contents ("CutVar").column_vector_value ();
      T.cutpoint = sm.contents ("CutPoint").column_vector_value ();
      T.leafval = sm.contents ("Leaf").matrix_value ();
      T.hascat = sm.isfield ("CatLeft")
                 && ! sm.contents ("CatLeft").isempty ();
      if (T.hascat)
        {
          T.catleft = sm.contents ("CatLeft").cell_value ();
          T.catright = sm.contents ("CatRight").cell_value ();
        }
      const octave_idx_type nn = T.children.rows ();
      if (T.cutvar.numel () != nn || T.cutpoint.numel () != nn
          || T.leafval.rows () != nn)
        error ("__shapleytree__: tree %d is not described consistently.",
               static_cast<int> (t + 1));
      if (K < 0)
        K = T.leafval.columns ();
      else if (T.leafval.columns () != K)
        error ("__shapleytree__: every tree must answer with as many columns.");
      collect_leaves (T);
    }

  // log-gamma of 1 .. M+1, so that lg[i] is log (i!)
  std::vector<double> lg (M + 2, 0.0);
  for (octave_idx_type i = 0; i <= M + 1; i++)
    lg[i] = std::lgamma (static_cast<double> (i) + 1.0);

  dim_vector dv (M, K, nq);
  NDArray out (dv, 0.0);
  const double invn = 1.0 / static_cast<double> (n);
  for (octave_idx_type p = 0; p < nq; p++)
    {
      Matrix phi (M, K, 0.0);
      RowVector q (M);
      for (octave_idx_type j = 0; j < M; j++)
        q(j) = Q(p, j);
      for (size_t t = 0; t < trees.size (); t++)
        tree_values (trees[t], X, q, M, K, lg, phi);
      for (octave_idx_type i = 0; i < M; i++)
        for (octave_idx_type k = 0; k < K; k++)
          out(i, k, p) = phi(i, k) * invn;
    }

  return ovl (out);
}

/*
%!error <Invalid call> __shapleytree__ ()
%!error <TREES must be a cell array of structures.> __shapleytree__ (1, 2, 3)
%!error <X must be a real numeric matrix.> __shapleytree__ ({1}, {}, 3)
%!error <Q must be a real numeric matrix.> __shapleytree__ ({1}, [1 2], {})
%!error <TREES must hold at least one tree.> __shapleytree__ ({}, [1 2], [1 2])
%!error <X and Q must have the same number of columns.> __shapleytree__ ({1}, [1 2], [1 2 3])

## A stump on the first of two predictors: holding it at the query point
## moves the answer by the whole deviation and the other predictor by nothing
%!test
%! T.Children = [2, 3; 0, 0; 0, 0];
%! T.CutVar = [1; 0; 0];
%! T.CutPoint = [0.5; NaN; NaN];
%! T.Leaf = [0; 10; 20];
%! X = [0, 7; 0, 8; 1, 9; 1, 6];
%! phi = __shapleytree__ ({T}, X, [1, 7]);
%! assert (size (phi), [2, 1]);
%! assert (phi(1), 5, 1e-12);
%! assert (phi(2), 0, 1e-12);

## Two stumps added together are the sum of what each gives on its own
%!test
%! T1.Children = [2, 3; 0, 0; 0, 0];
%! T1.CutVar = [1; 0; 0];
%! T1.CutPoint = [0.5; NaN; NaN];
%! T1.Leaf = [0; 10; 20];
%! T2.Children = [2, 3; 0, 0; 0, 0];
%! T2.CutVar = [2; 0; 0];
%! T2.CutPoint = [7.5; NaN; NaN];
%! T2.Leaf = [0; 1; 3];
%! X = [0, 7; 0, 8; 1, 9; 1, 6];
%! a = __shapleytree__ ({T1}, X, [1, 7]);
%! b = __shapleytree__ ({T2}, X, [1, 7]);
%! both = __shapleytree__ ({T1, T2}, X, [1, 7]);
%! assert (both, a + b, 1e-12);

## A cut on the levels of a categorical predictor sends the levels it lists
%!test
%! T.Children = [2, 3; 0, 0; 0, 0];
%! T.CutVar = [1; 0; 0];
%! T.CutPoint = [NaN; NaN; NaN];
%! T.CatLeft = {[1, 2]; []; []};
%! T.CatRight = {[3]; []; []};
%! T.Leaf = [0; 4; 10];
%! X = [1; 2; 3; 3];
%! phi = __shapleytree__ ({T}, X, 3);
%! assert (phi, 3, 1e-12);

## One page of values per query point
%!test
%! T.Children = [2, 3; 0, 0; 0, 0];
%! T.CutVar = [1; 0; 0];
%! T.CutPoint = [0.5; NaN; NaN];
%! T.Leaf = [0; 10; 20];
%! X = [0, 7; 0, 8; 1, 9; 1, 6];
%! phi = __shapleytree__ ({T}, X, [1, 7; 0, 7]);
%! assert (size (phi), [2, 1, 2]);
%! assert (phi(1,1,1), 5, 1e-12);
%! assert (phi(1,1,2), -5, 1e-12);
*/
