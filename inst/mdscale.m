## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{Y} =} mdscale (@var{D}, @var{p})
## @deftypefnx {statistics} {[@var{Y}, @var{stress}] =} mdscale (@var{D}, @var{p})
## @deftypefnx {statistics} {[@var{Y}, @var{stress}, @var{disparities}] =} mdscale (@var{D}, @var{p})
## @deftypefnx {statistics} {[@dots{}] =} mdscale (@dots{}, @var{Name}, @var{Value})
##
## Nonclassical (metric and nonmetric) multidimensional scaling.
##
## @code{@var{Y} = mdscale (@var{D}, @var{p})} takes a matrix of dissimilarities
## @var{D} and returns a configuration @var{Y} of @var{n} points in @var{p}
## dimensions (an @math{@var{n} * @var{p} } matrix) whose interpoint distances
## approximate @var{D}, by minimizing a stress criterion.  @var{D} may be given
## either as a full @math{@var{n} * @var{n}} symmetric matrix with zero diagonal,
## or as the vector of the @math{@var{n} (@var{n} - 1) / 2} upper-triangle
## dissimilarities returned by @code{pdist}.
##
## @code{[@var{Y}, @var{stress}, @var{disparities}] = mdscale (@dots{})} also
## returns the final value of the @var{stress} criterion and the
## @var{disparities} (the transformed dissimilarities the distances are fitted
## to).  For the nonmetric criteria the disparities are the monotone
## (isotonic) regression of the dissimilarities onto the distances; for the
## metric criteria they are the dissimilarities themselves.
##
## Name/Value pairs:
##
## @table @asis
## @item @qcode{'Criterion'}
## The goodness-of-fit criterion to minimize, one of:
##
## @table @asis
## @item @qcode{'stress'} (default)
## Kruskal's normalized stress-1, @math{sqrt (sum ((d - dhat)^2) / sum (d^2))},
## computed from disparities @var{dhat} (nonmetric).
##
## @item @qcode{'sstress'}
## Squared stress, @math{sqrt (sum ((d^2 - dhat^2)^2) / sum (d^4))} (nonmetric).
##
## @item @qcode{'metricstress'}
## Metric stress, @math{sqrt (sum ((d - delta)^2) / sum (delta^2))}, fitting the
## dissimilarities @var{delta} directly.
##
## @item @qcode{'metricsstress'}
## Metric squared stress, @math{sqrt (sum ((d^2 - delta^2)^2) / sum (delta^4))}.
##
## @item @qcode{'sammon'}
## Sammon's nonlinear mapping criterion,
## @math{(1 / sum (delta)) sum ((d - delta)^2 / delta)}.
##
## @item @qcode{'strain'}
## The classical scaling criterion; equivalent to @code{cmdscale}.
## @end table
##
## @item @qcode{'Weights'}
## A matrix or vector of nonnegative weights, the same size as @var{D}, weighting
## each dissimilarity in the criterion.
##
## @item @qcode{'Start'}
## The initial configuration: @qcode{'cmdscale'} (default, classical scaling),
## @qcode{'random'}, or an explicit @math{@var{n} * @var{p} } matrix.
##
## @item @qcode{'Replicates'}
## The number of times to repeat the minimization from different starting points,
## keeping the best (lowest-stress) result.  The default is 1.
##
## @item @qcode{'Options'}
## A structure of algorithm options (as returned by @code{statset}) whose
## @qcode{MaxIter}, @qcode{TolFun}, and @qcode{TolX} fields control the iterative
## minimization.
## @end table
##
## @subheading Non-uniqueness of the solution
##
## A stress-minimizing configuration is defined only up to a translation,
## rotation, and reflection, because these leave all interpoint distances (and
## hence the stress) unchanged.  @code{mdscale} removes this ambiguity by
## returning @var{Y} centred at the origin and rotated to its principal axes
## (with the largest-magnitude coordinate on each axis made positive), matching
## the convention used by MATLAB.
##
## Beyond that rigid ambiguity, the @strong{nonmetric} criteria (@qcode{'stress'}
## and @qcode{'sstress'}) are non-convex and typically have several local minima;
## the one reached depends on the starting configuration and the details of the
## optimizer.  As a result the returned configuration for these criteria may
## differ from the one another program (including MATLAB) reports even when the
## @strong{stress value} agrees, and different runs may find configurations with
## slightly different stress.  Use @qcode{'Replicates'} with a @qcode{'random'}
## start to search for a lower-stress solution.  The metric criteria and
## @qcode{'strain'} have an essentially unique solution and are reproducible up to
## the rigid ambiguity above.
##
## @seealso{cmdscale, pdist, squareform, procrustes, statset}
## @end deftypefn

function [Y, stress, disparities] = mdscale (D, p, varargin)

  ## Input validation
  if (nargin < 2)
    print_usage ();
  endif

  ## Resolve D into the vector of dissimilarities DELTA and the count N
  if (isvector (D))
    delta = D(:);
    m = numel (delta);
    n = (1 + sqrt (1 + 8 * m)) / 2;
    if (n != fix (n))
      error (strcat ("mdscale: the length of a dissimilarity vector D must", ...
                     " be N*(N-1)/2 for some integer N."));
    endif
  elseif (ismatrix (D) && ndims (D) == 2 && rows (D) == columns (D))
    n = rows (D);
    delta = D(tril (true (n), -1));
    m = numel (delta);
  else
    error ("mdscale: D must be a dissimilarity vector or a square matrix.");
  endif
  if (! (isnumeric (delta) && isreal (delta)) || any (delta < 0))
    error ("mdscale: D must contain real, nonnegative dissimilarities.");
  endif

  if (! (isscalar (p) && isnumeric (p) && p >= 1 && p == fix (p) && p < n))
    error ("mdscale: P must be a positive integer smaller than N.");
  endif

  ## Defaults and Name/Value parsing
  Criterion = "stress";
  Weights = [];
  Start = "cmdscale";
  Replicates = 1;
  opts = struct ("MaxIter", 1000, "TolFun", 1e-6, "TolX", 1e-6);
  args = varargin;
  if (mod (numel (args), 2) != 0)
    error ("mdscale: Name/Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (args)
    name = args{k};
    val = args{k+1};
    if (! (ischar (name) && isrow (name)))
      error ("mdscale: parameter names must be character vectors.");
    endif
    switch (lower (name))
      case 'criterion'
        Criterion = lower (val);
      case 'weights'
        Weights = val;
      case 'start'
        Start = val;
      case 'replicates'
        Replicates = val;
      case 'options'
        opts = merge_statset (opts, val);
      otherwise
        error ("mdscale: unknown parameter name '%s'.", name);
    endswitch
  endfor

  critlist = {'stress', 'sstress', 'metricstress', 'metricsstress', ...
              'sammon', 'strain'};
  if (! (ischar (Criterion) && any (strcmp (Criterion, critlist))))
    error ("mdscale: unknown 'Criterion' '%s'.", Criterion);
  endif
  nonmetric = any (strcmp (Criterion, {'stress', 'sstress'}));

  ## Weights vector aligned with DELTA
  if (isempty (Weights))
    w = ones (m, 1);
  elseif (isvector (Weights) && numel (Weights) == m)
    w = Weights(:);
  elseif (ismatrix (Weights) && isequal (size (Weights), [n, n]))
    w = Weights(tril (true (n), -1));
  else
    error ("mdscale: 'Weights' must match the size of D.");
  endif
  if (! (isnumeric (w) && isreal (w)) || any (w < 0))
    error ("mdscale: 'Weights' must be real and nonnegative.");
  endif

  if (! (isscalar (Replicates) && Replicates >= 1 && Replicates == fix (Replicates)))
    error ("mdscale: 'Replicates' must be a positive integer.");
  endif

  ## The 'strain' criterion is classical scaling: a direct eigen-solution.
  if (strcmp (Criterion, "strain"))
    Y = classical_scaling (delta, n, p);
    Y = orient (Y);
    stress = strain_value (delta, n, Y);
    disparities = delta;
    return;
  endif

  ## Pair index list, in pdist (upper-triangle) order
  P = zeros (m, 2);
  r = 1;
  for i = 1:n-1
    for j = i+1:n
      P(r,:) = [i, j];
      r++;
    endfor
  endfor

  ## Map the criterion name to the id used by crit_grad/eval_stress
  critmap = {'metricstress', 'metricsstress', 'sammon', 'stress', 'sstress'};
  critid = find (strcmp (Criterion, critmap));

  ## Run the minimization once per replicate, keeping the best result
  bestY = [];
  bestf = Inf;
  for rep = 1:Replicates
    Y0 = initial_config (Start, delta, n, p, rep);
    [Yr, fr] = minimize_stress (Y0, critid, delta, w, n, p, P, nonmetric, opts);
    if (fr < bestf)
      bestf = fr;
      bestY = Yr;
    endif
  endfor

  Y = orient (bestY);
  [stress, disparities] = eval_stress (Y, critid, delta, w, n, p, P, nonmetric);

endfunction

## ---------------------------------------------------------------------------
## Build the starting configuration.
function Y0 = initial_config (Start, delta, n, p, rep)

  if (ischar (Start))
    switch (lower (Start))
      case 'cmdscale'
        if (rep == 1)
          Y0 = classical_scaling (delta, n, p);
        else
          Y0 = randn (n, p);
        endif
      case 'random'
        Y0 = randn (n, p);
      otherwise
        error ("mdscale: unknown 'Start' '%s'.", Start);
    endswitch
  elseif (isnumeric (Start) && isequal (size (Start), [n, p]))
    Y0 = Start;
    if (rep > 1)
      Y0 += randn (n, p);
    endif
  else
    error ("mdscale: 'Start' must be 'cmdscale', 'random', or an N-by-P matrix.");
  endif

endfunction

## ---------------------------------------------------------------------------
## Classical scaling (the 'strain' solution and default start).
function Y = classical_scaling (delta, n, p)

  Dsq = squareform (delta) .^ 2;
  J = eye (n) - ones (n) / n;
  B = -0.5 * J * Dsq * J;
  [V, L] = eig ((B + B') / 2);
  [ev, idx] = sort (diag (L), "descend");
  V = V(:,idx);
  ev = max (ev(1:p), 0);
  Y = V(:,1:p) .* sqrt (ev)';

endfunction

## ---------------------------------------------------------------------------
## The 'strain' criterion value.
function s = strain_value (delta, n, Y)

  Dsq = squareform (delta) .^ 2;
  s = sqrt (sum (sumsq (-0.5 * Dsq - Y * Y')));

endfunction

## ---------------------------------------------------------------------------
## Centre at the origin and rotate to principal axes with a sign convention.
function Y = orient (Y)

  Y = Y - mean (Y);
  [~, ~, V] = svd (Y, 0);
  Y = Y * V;
  for j = 1:columns (Y)
    [~, im] = max (abs (Y(:,j)));
    if (Y(im,j) < 0)
      Y(:,j) = -Y(:,j);
    endif
  endfor

endfunction

## ---------------------------------------------------------------------------
## Isotonic (monotone) regression of X in the order of DELTA, with primary
## tie handling (ties in DELTA broken by X), by pool-adjacent-violators.
function dhat = pav (x, delta)

  [~, ord] = sortrows ([delta, x]);
  xo = x(ord);
  m = numel (xo);
  vals = xo;
  wts = ones (m, 1);
  cnt = ones (m, 1);
  K = m;
  k = 1;
  while (k < K)
    if (vals(k) > vals(k+1))
      vals(k) = (wts(k) * vals(k) + wts(k+1) * vals(k+1)) / (wts(k) + wts(k+1));
      wts(k) += wts(k+1);
      cnt(k) += cnt(k+1);
      vals(k+1) = [];
      wts(k+1) = [];
      cnt(k+1) = [];
      K--;
      if (k > 1)
        k--;
      endif
    else
      k++;
    endif
  endwhile
  xexp = zeros (m, 1);
  pos = 1;
  for b = 1:K
    xexp(pos:pos+cnt(b)-1) = vals(b);
    pos += cnt(b);
  endfor
  dhat = zeros (m, 1);
  dhat(ord) = xexp;

endfunction

## ---------------------------------------------------------------------------
## Criterion value and gradient at a configuration.  For the nonmetric criteria
## the disparities are recomputed and held fixed for the gradient (Kruskal).
function [f, g] = crit_grad (Yv, critid, delta, w, n, p, P, nonmetric)

  Y = reshape (Yv, n, p);
  dif = Y(P(:,1),:) - Y(P(:,2),:);
  d = sqrt (sum (dif .^ 2, 2));
  d(d == 0) = eps;

  switch (critid)
    case 1   ## metricstress
      T = sum (w .* delta .^ 2);
      res = d - delta;
      N = sum (w .* res .^ 2);
      f = sqrt (N / T);
      dfdd = (w .* res) / (f * T);
    case 2   ## metricsstress
      T = sum (w .* delta .^ 4);
      res = d .^ 2 - delta .^ 2;
      N = sum (w .* res .^ 2);
      f = sqrt (N / T);
      dfdd = (2 * w .* d .* res) / (f * T);
    case 3   ## sammon
      S = sum (w .* delta);
      res = d - delta;
      f = sum (w .* res .^ 2 ./ delta) / S;
      dfdd = (2 * w .* res ./ delta) / S;
    case 4   ## stress (nonmetric)
      dh = pav (d, delta);
      res = d - dh;
      T = sum (w .* d .^ 2);
      N = sum (w .* res .^ 2);
      f = sqrt (N / T);
      dfdd = (w .* res * T - N * (w .* d)) / (f * T ^ 2);
    case 5   ## sstress (nonmetric)
      dh = pav (d, delta);
      res = d .^ 2 - dh .^ 2;
      T = sum (w .* d .^ 4);
      N = sum (w .* res .^ 2);
      f = sqrt (N / T);
      dfdd = (2 * w .* d .* res * T - N * (2 * w .* d .^ 3)) / (f * T ^ 2);
  endswitch

  gg = (dfdd ./ d) .* dif;
  g = zeros (n, p);
  for k = 1:size (P, 1)
    g(P(k,1),:) += gg(k,:);
    g(P(k,2),:) -= gg(k,:);
  endfor
  g = g(:);

endfunction

## ---------------------------------------------------------------------------
## Gradient descent with backtracking (Armijo) line search.
function [Y, f] = minimize_stress (Y0, critid, delta, w, n, p, P, nonmetric, opts)

  Yv = Y0(:);
  [f, g] = crit_grad (Yv, critid, delta, w, n, p, P, nonmetric);
  al = 1;
  for it = 1:opts.MaxIter
    al = min (al * 2, 1e3);
    gg = g' * g;
    fn = f;
    Yn = Yv;
    for ls = 1:60
      Yn = Yv - al * g;
      fn = crit_grad (Yn, critid, delta, w, n, p, P, nonmetric);
      if (fn < f - 1e-4 * al * gg)
        break;
      endif
      al /= 2;
    endfor
    step = max (abs (Yv - Yn));
    Yv = Yn;
    [f, g] = crit_grad (Yv, critid, delta, w, n, p, P, nonmetric);
    if (abs (fn - f) <= opts.TolFun * max (1, f) && step <= opts.TolX)
      break;
    endif
  endfor
  Y = reshape (Yv, n, p);

endfunction

## ---------------------------------------------------------------------------
## Final stress value and disparities at the returned configuration.
function [stress, disparities] = eval_stress (Y, critid, delta, w, n, p, P, nonmetric)

  dif = Y(P(:,1),:) - Y(P(:,2),:);
  d = sqrt (sum (dif .^ 2, 2));
  if (nonmetric)
    disparities = pav (d, delta);
  else
    disparities = delta;
  endif
  dd = disparities;
  switch (critid)
    case {1}
      stress = sqrt (sum (w .* (d - dd) .^ 2) / sum (w .* dd .^ 2));
    case {2}
      stress = sqrt (sum (w .* (d.^2 - dd.^2) .^ 2) / sum (w .* dd .^ 4));
    case 3
      stress = sum (w .* (d - dd) .^ 2 ./ dd) / sum (w .* dd);
    case 4
      stress = sqrt (sum (w .* (d - dd) .^ 2) / sum (w .* d .^ 2));
    case 5
      stress = sqrt (sum (w .* (d.^2 - dd.^2) .^ 2) / sum (w .* d .^ 4));
  endswitch

endfunction

## ---------------------------------------------------------------------------
## Copy the recognised statset fields from S into OPTS.
function opts = merge_statset (opts, s)

  if (! isstruct (s))
    error ("mdscale: 'Options' must be a structure.");
  endif
  fn = fieldnames (s);
  for k = 1:numel (fn)
    val = s.(fn{k});
    if (isempty (val))
      continue;
    endif
    switch (lower (fn{k}))
      case 'maxiter'
        opts.MaxIter = val;
      case 'tolfun'
        opts.TolFun = val;
      case 'tolx'
        opts.TolX = val;
    endswitch
  endfor

endfunction

%!demo
%! ## Recover a 2-D map of 8 points from their pairwise distances.
%! rng = [0 0; 5 0; 5 4; 0 4; 8 2; -3 2; 2 7; 2 -3];
%! D = pdist (rng);
%! [Y, stress] = mdscale (D, 2, 'Criterion', 'metricstress');
%! stress
%! ## Y reproduces the pairwise distances closely.
%! max (abs (pdist (Y)' - D'))

## Reference values below are from MATLAB R2023b, for
## C = [0 0 0; 5 1 2; 2 6 1; 1 2 7; 7 3 4; 3 8 5; 6 4 9; 4 9 3];
## D = pdist (C);  Y0 = [1 -1;-1 1;2 2;-2 -2;3 1;1 3;-3 -1;-1 -3];
## The metric criteria and 'strain' reproduce MATLAB's configuration; the
## nonmetric criteria reproduce its stress value (see the non-uniqueness note).

%!shared C, D, Y0, opt
%! C = [0 0 0; 5 1 2; 2 6 1; 1 2 7; 7 3 4; 3 8 5; 6 4 9; 4 9 3];
%! D = pdist (C);
%! Y0 = [1 -1; -1 1; 2 2; -2 -2; 3 1; 1 3; -3 -1; -1 -3];
%! opt = struct ("MaxIter", 1000, "TolFun", 1e-10, "TolX", 1e-10);

%!test
%! [Y, s] = mdscale (D, 2, "Criterion", "metricstress", "Start", Y0, "Options", opt);
%! Yref = [ 6.810767239317778, -0.733620649959338; ...
%!          3.071301239764592,  1.094132882710339; ...
%!          1.065313739387755, -3.760130869234625; ...
%!          0.623478915333739,  4.891508042741564; ...
%!         -0.674391648191895,  1.443596654184263; ...
%!         -3.493671563292435, -2.348293046078109; ...
%!         -4.559656269176627,  3.674896334304512; ...
%!         -2.843141653142908, -4.262089348668606];
%! assert_equal (Y, Yref, 1e-4);
%! assert_equal (s, 0.144654108154698, 1e-8);

%!test
%! [Y, s] = mdscale (D, 2, "Criterion", "metricsstress", "Start", Y0, "Options", opt);
%! Yref = [ 6.635327214692762,  0.792856204849063; ...
%!          3.228563095536443, -1.334104353307981; ...
%!          1.003358563030293,  3.760704781755182; ...
%!          0.699645585641474, -3.784555012990314; ...
%!         -1.079052104978016, -2.376723256450389; ...
%!         -3.223960839101638,  2.563657848313004; ...
%!         -4.268551283289842, -3.733948071746559; ...
%!         -2.995330231531478,  4.112111859577994];
%! assert_equal (Y, Yref, 1e-4);
%! assert_equal (s, 0.190937534497622, 1e-8);

%!test
%! [Y, s] = mdscale (D, 2, "Criterion", "sammon", "Start", Y0, "Options", opt);
%! Yref = [ 6.933279805420950, -0.419508440326684; ...
%!          3.008163425596839,  1.373357387532562; ...
%!          1.273453983125238, -3.707332917059666; ...
%!          0.101742369038229,  5.027798678525714; ...
%!         -0.450738847427944,  1.372181450811530; ...
%!         -3.390706238135306, -2.442042299704407; ...
%!         -4.908058626481512,  3.317873521951261; ...
%!         -2.567135871136496, -4.522327381730312];
%! assert_equal (Y, Yref, 1e-4);
%! assert_equal (s, 0.023119816166763, 1e-8);

%!test
%! ## 'strain' is classical scaling: configuration and value are reproducible.
%! [Y, s] = mdscale (D, 2, "Criterion", "strain", "Start", Y0, "Options", opt);
%! Yref = [ 6.572163580404339, -0.496338025399502; ...
%!          2.719997306646783,  1.496763413554223; ...
%!          1.036596577326041, -3.598999625596409; ...
%!          0.538361534877868,  2.777882896586096; ...
%!         -0.658048085951181,  1.901764558503386; ...
%!         -3.170613455894927, -2.220122324409574; ...
%!         -4.005398518313616,  4.052005597122571; ...
%!         -3.033058939095310, -3.912956490360791];
%! assert_equal (Y, Yref, 1e-3);
%! assert_equal (s, 1.930074461959919e+02, 1e-6);

%!test
%! ## Nonmetric 'stress': the stress value matches MATLAB (configuration is
%! ## only defined up to the local minimum -- see the non-uniqueness note).
%! [Y, s, dsp] = mdscale (D, 2, "Criterion", "stress", "Start", Y0, "Options", opt);
%! assert_equal (s, 0.078131057459269, 1e-6);
%! assert_equal (size (Y), [8, 2]);
%! assert_equal (numel (dsp), 28);
%! ## disparities are monotone in the dissimilarities
%! [~, ord] = sort (D(:));
%! assert_equal (all (diff (dsp(ord)) >= -1e-9), true);

%!test
%! ## Nonmetric 'sstress' stress value.
%! [Y, s] = mdscale (D, 2, "Criterion", "sstress", "Start", Y0, "Options", opt);
%! assert_equal (s <= 0.089488341028430 + 1e-6, true);

%!test
%! ## Accepts a square dissimilarity matrix as well as a vector.
%! [Y1, s1] = mdscale (D, 2, "Criterion", "metricstress", "Start", Y0, "Options", opt);
%! [Y2, s2] = mdscale (squareform (D), 2, "Criterion", "metricstress", ...
%!                     "Start", Y0, "Options", opt);
%! assert_equal (Y1, Y2, 1e-12);

%!test
%! ## Returned configuration is centred and on principal axes.
%! Y = mdscale (D, 2, "Criterion", "metricstress", "Start", Y0, "Options", opt);
%! assert_equal (mean (Y), [0, 0], 1e-10);
%! c = cov (Y);
%! assert_equal (c(1,2), 0, 1e-8);

## Test input validation
%!error<Invalid call to mdscale> mdscale (D)
%!error<mdscale: the length of a dissimilarity vector D must be> mdscale (1:4, 2)
%!error<mdscale: D must be a dissimilarity vector or a square matrix.> ...
%! mdscale (ones (3, 4), 2)
%!error<mdscale: D must contain real, nonnegative dissimilarities.> ...
%! mdscale ([1, -2, 3], 1)
%!error<mdscale: P must be a positive integer smaller than N.> mdscale (D, 8)
%!error<mdscale: P must be a positive integer smaller than N.> mdscale (D, 0)
%!error<mdscale: unknown 'Criterion' 'foo'.> mdscale (D, 2, "Criterion", "foo")
%!error<mdscale: unknown parameter name 'bogus'.> mdscale (D, 2, "bogus", 1)
%!error<mdscale: 'Weights' must match the size of D.> ...
%! mdscale (D, 2, "Weights", [1, 2, 3])
%!error<mdscale: 'Replicates' must be a positive integer.> ...
%! mdscale (D, 2, "Replicates", 0)
%!error<mdscale: unknown 'Start' 'bad'.> mdscale (D, 2, "Start", "bad")
%!error<mdscale: 'Options' must be a structure.> mdscale (D, 2, "Options", 5)
