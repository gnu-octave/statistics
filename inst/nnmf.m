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
## @deftypefn  {statistics} {[@var{W}, @var{H}] =} nnmf (@var{A}, @var{K})
## @deftypefnx {statistics} {[@var{W}, @var{H}, @var{D}] =} nnmf (@var{A}, @var{K})
## @deftypefnx {statistics} {[@dots{}] =} nnmf (@dots{}, @var{Name}, @var{Value})
##
## Nonnegative matrix factorization.
##
## @code{[@var{W}, @var{H}] = nnmf (@var{A}, @var{K})} factors the nonnegative
## @math{N * M} matrix @var{A} into nonnegative factors @var{W} (@math{N *
## @var{K}}) and @var{H} (@math{@var{K} * M}) whose product approximates
## @var{A}, by minimizing the root-mean-square residual between @var{A} and
## @code{@var{W} * @var{H}}.  @var{K}, the number of factors, is typically smaller
## than @math{N} and @math{M}.
##
## @code{[@var{W}, @var{H}, @var{D}] = nnmf (@dots{})} also returns the
## root-mean-square residual @var{D}, that is @code{norm (@var{A} - @var{W} *
## @var{H}, "fro") / sqrt (N * M)}.
##
## The factorization is not unique: the returned factors are normalized so that
## the rows of @var{H} have unit length, and the columns of @var{W} (and the
## corresponding rows of @var{H}) are ordered by decreasing length of the
## columns of @var{W}.  Because the objective is not convex, the iteration
## converges to a local minimum that depends on the starting point; use
## @qcode{'Replicates'} to try several random starts and keep the best.
##
## Name/Value pairs:
##
## @table @asis
## @item @qcode{'Algorithm'}
## @qcode{'als'} (default) for alternating least squares, or @qcode{'mult'} for
## multiplicative updates.  Alternating least squares usually converges faster
## and more reliably; multiplicative updates are more sensitive to the starting
## point.
##
## @item @qcode{'W0'}
## An @math{N * @var{K}} initial value for @var{W}.
##
## @item @qcode{'H0'}
## A @math{@var{K} * M} initial value for @var{H}.
##
## @item @qcode{'Replicates'}
## The number of times to repeat the factorization from new random starting
## points, keeping the result with the smallest residual.  The default is 1.
## Ignored for a starting point fixed by both @qcode{'W0'} and @qcode{'H0'}.
##
## @item @qcode{'Options'}
## A structure of algorithm options (as returned by @code{statset}) whose
## @qcode{MaxIter}, @qcode{TolFun}, and @qcode{TolX} fields control the iteration.
## @end table
##
## @seealso{pca, statset}
## @end deftypefn

function [W, H, D] = nnmf (A, K, varargin)

  ## Input validation
  if (nargin < 2)
    print_usage ();
  endif
  if (! (isnumeric (A) && ismatrix (A) && ndims (A) == 2))
    error ("nnmf: A must be a numeric matrix.");
  endif
  if (! (isreal (A) && all (isfinite (A(:)))))
    error ("nnmf: A must be real and finite.");
  endif

  [n, m] = size (A);

  if (! (isscalar (K) && isnumeric (K) && K >= 1 && K == fix (K)))
    error ("nnmf: K must be a positive integer.");
  endif

  ## Defaults and Name/Value parsing
  algorithm = "als";
  W0 = [];
  H0 = [];
  Replicates = 1;
  opts = struct ("MaxIter", 100, "TolFun", 1e-4, "TolX", 1e-4);
  args = varargin;
  if (mod (numel (args), 2) != 0)
    error ("nnmf: Name/Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (args)
    name = args{k};
    val = args{k+1};
    if (! (ischar (name) && isrow (name)))
      error ("nnmf: parameter names must be character vectors.");
    endif
    switch (lower (name))
      case 'algorithm'
        algorithm = lower (val);
      case 'w0'
        W0 = val;
      case 'h0'
        H0 = val;
      case 'replicates'
        Replicates = val;
      case 'options'
        opts = merge_statset (opts, val);
      otherwise
        error ("nnmf: unknown parameter name '%s'.", name);
    endswitch
  endfor

  if (! any (strcmp (algorithm, {'als', 'mult'})))
    error ("nnmf: 'Algorithm' must be 'als' or 'mult'.");
  endif
  if (! isempty (W0) && ! isequal (size (W0), [n, K]))
    error ("nnmf: 'W0' must be an %d-by-%d matrix.", n, K);
  endif
  if (! isempty (H0) && ! isequal (size (H0), [K, m]))
    error ("nnmf: 'H0' must be an %d-by-%d matrix.", K, m);
  endif
  if (! (isscalar (Replicates) && Replicates >= 1 && Replicates == fix (Replicates)))
    error ("nnmf: 'Replicates' must be a positive integer.");
  endif

  ## With both factors fixed there is a single starting point to run
  fixedStart = ! isempty (W0) && ! isempty (H0);
  if (fixedStart)
    Replicates = 1;
  endif

  ## Run the factorization once per replicate, keeping the best
  bestD = Inf;
  bestW = [];
  bestH = [];
  scal = sqrt (mean (A(:)) / K);
  for rep = 1:Replicates
    if (isempty (W0))
      W = rand (n, K) * scal;
    else
      W = W0;
    endif
    if (isempty (H0))
      H = rand (K, m) * scal;
    else
      H = H0;
    endif
    [W, H, d] = factorize (A, W, H, algorithm, opts, n, m);
    if (d < bestD)
      bestD = d;
      bestW = W;
      bestH = H;
    endif
  endfor

  ## Normalize: unit-length rows of H, columns ordered by decreasing W length
  hlen = sqrt (sum (bestH .^ 2, 2));
  hlen(hlen == 0) = 1;
  H = bestH ./ hlen;
  W = bestW .* hlen';
  [~, idx] = sort (sqrt (sum (W .^ 2, 1)), "descend");
  W = W(:,idx);
  H = H(idx,:);
  D = bestD;

endfunction

## ---------------------------------------------------------------------------
## Iterate one algorithm from a starting point to convergence.
function [W, H, dnorm] = factorize (A, W, H, algorithm, opts, n, m)

  sqrteps = sqrt (eps);
  dnorm0 = Inf;
  for iter = 1:opts.MaxIter
    Wp = W;
    Hp = H;
    if (strcmp (algorithm, "als"))
      H = max (0, W \ A);
      W = max (0, A / H);
    else
      H = H .* (W' * A) ./ ((W' * W) * H + eps);
      W = W .* (A * H') ./ (W * (H * H') + eps);
    endif

    dnorm = norm (A - W * H, "fro") / sqrt (n * m);
    if (iter > 1)
      dw = max (max (abs (W - Wp))) / (sqrteps + max (max (abs (Wp))));
      dh = max (max (abs (H - Hp))) / (sqrteps + max (max (abs (Hp))));
      if (max (dw, dh) <= opts.TolX)
        break;
      endif
      if (dnorm0 - dnorm <= opts.TolFun * max (1, dnorm0))
        break;
      endif
    endif
    dnorm0 = dnorm;
  endfor

endfunction

## ---------------------------------------------------------------------------
## Copy the recognised statset fields from S into OPTS.
function opts = merge_statset (opts, s)

  if (! isstruct (s))
    error ("nnmf: 'Options' must be a structure.");
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
%! ## Factor a nonnegative matrix into two rank-2 nonnegative factors.
%! A = [1, 2, 3; 2, 4, 6; 3, 5, 7; 4, 8, 12];
%! [W, H, D] = nnmf (A, 2);
%! D
%! ## W * H approximates A.
%! W * H

## Reference values below are from MATLAB R2023b for
## A = reshape (mod ((1:20)*7, 13), 5, 4) + 1, with fixed starting factors.

%!shared A, W0, H0, opt
%! A = reshape (mod ((1:20)*7, 13), 5, 4) + 1;
%! W0 = reshape (mod ((1:10)*3, 7), 5, 2) + 1;
%! H0 = reshape (mod ((1:8)*3, 7), 2, 4) + 1;
%! opt = struct ("MaxIter", 500, "TolFun", 1e-12, "TolX", 1e-12);

%!test
%! [W, H, D] = nnmf (A, 2, "Algorithm", "als", "W0", W0, "H0", H0, "Options", opt);
%! Wref = [11.808229123411982,  9.049993647167147; ...
%!          3.212266148150507, 12.092875848976247; ...
%!         13.853561873197318,  1.116057457895120; ...
%!          4.618299368228214, 13.176478328920602; ...
%!         15.259595093275024,  2.199659937839474];
%! Href = [0.657913952127586, 0.187993469755446, 0.140048330155150, 0.715677407877266; ...
%!         0,                 0.730839581652543, 0.680241675322826, 0.056078240378338];
%! assert_equal (W, Wref, 1e-4);
%! assert_equal (H, Href, 1e-5);
%! assert_equal (D, 1.882172868745145, 1e-8);

%!test
%! [W, H, D] = nnmf (A, 2, "Algorithm", "mult", "W0", W0, "H0", H0, "Options", opt);
%! Wref = [11.221742977568461,  8.713876463682844; ...
%!          2.176215667360621, 12.316408085448018; ...
%!         14.030193052240865,  0.411019822854590; ...
%!          3.511868529923804, 13.360184914158390; ...
%!         15.365861312580710,  1.454761030110911];
%! Href = [0.647958006689051, 0.222457058550383, 0.172654918585345, 0.707710080299095; ...
%!         0.057102252901681, 0.727300231504821, 0.673915159707622, 0.116670748188381];
%! assert_equal (W, Wref, 1e-3);
%! assert_equal (H, Href, 1e-4);
%! assert_equal (D, 1.882172868753092, 1e-8);

%!test
%! ## Factors are nonnegative, correctly sized, and H rows are unit length.
%! [W, H] = nnmf (A, 2, "W0", W0, "H0", H0, "Options", opt);
%! assert_equal (size (W), [5, 2]);
%! assert_equal (size (H), [2, 4]);
%! assert_equal (all (W(:) >= 0) && all (H(:) >= 0), true);
%! assert_equal (sqrt (sum (H .^ 2, 2)), [1; 1], 1e-10);

%!test
%! ## Columns of W are ordered by decreasing length.
%! [W, H] = nnmf (A, 2, "W0", W0, "H0", H0, "Options", opt);
%! assert_equal (norm (W(:,1)) >= norm (W(:,2)), true);

%!test
%! ## The residual D matches the direct computation.
%! [W, H, D] = nnmf (A, 2, "W0", W0, "H0", H0, "Options", opt);
%! assert_equal (D, norm (A - W * H, "fro") / sqrt (numel (A)), 1e-12);

## Test input validation
%!error<Invalid call to nnmf> nnmf (ones (4, 3))
%!error<nnmf: A must be a numeric matrix.> nnmf ({1, 2}, 1)
%!error<nnmf: A must be real and finite.> nnmf ([1, Inf; 2, 3], 1)
%!error<nnmf: A must be real and finite.> nnmf ([1+2i, 3; 4, 5], 1)
%!error<nnmf: K must be a positive integer.> nnmf (ones (4, 3), 0)
%!error<nnmf: K must be a positive integer.> nnmf (ones (4, 3), 1.5)
%!error<nnmf: 'Algorithm' must be 'als' or 'mult'.> ...
%! nnmf (ones (4, 3), 2, "Algorithm", "foo")
%!error<nnmf: 'W0' must be an 4-by-2 matrix.> ...
%! nnmf (ones (4, 3), 2, "W0", ones (3, 3))
%!error<nnmf: 'H0' must be an 2-by-3 matrix.> ...
%! nnmf (ones (4, 3), 2, "H0", ones (2, 2))
%!error<nnmf: 'Replicates' must be a positive integer.> ...
%! nnmf (ones (4, 3), 2, "Replicates", 0)
%!error<nnmf: unknown parameter name 'bogus'.> nnmf (ones (4, 3), 2, "bogus", 1)
%!error<nnmf: 'Options' must be a structure.> nnmf (ones (4, 3), 2, "Options", 5)
