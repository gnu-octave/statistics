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
## @deftypefn  {statistics} {@var{Y} =} tsne (@var{X})
## @deftypefnx {statistics} {[@var{Y}, @var{loss}] =} tsne (@var{X})
## @deftypefnx {statistics} {[@dots{}] =} tsne (@dots{}, @var{Name}, @var{Value})
##
## t-distributed stochastic neighbor embedding (t-SNE).
##
## @code{@var{Y} = tsne (@var{X})} embeds the @math{N * P} data matrix @var{X}
## (rows are observations) into a low-dimensional space and returns the
## @math{N * @var{NumDimensions}} matrix @var{Y} of embedded points, whose
## pairwise (Student-t) affinities approximate the Gaussian affinities of the
## rows of @var{X}.
##
## @code{[@var{Y}, @var{loss}] = tsne (@dots{})} also returns the
## Kullback-Leibler divergence @var{loss} between the two affinity
## distributions at the returned embedding.
##
## Name/Value pairs:
##
## @table @asis
## @item @qcode{'Algorithm'}
## Only @qcode{'exact'} is supported (MATLAB's approximate @qcode{'barneshut'}
## algorithm is not available).
##
## @item @qcode{'Distance'}
## The distance metric used for the high-dimensional affinities, as accepted by
## @code{pdist} (default @qcode{'euclidean'}).
##
## @item @qcode{'NumDimensions'}
## The dimension of the embedding @var{Y} (default @code{min (P, 2)}).
##
## @item @qcode{'NumPCAComponents'}
## If positive, reduce @var{X} to this many principal components before embedding
## (default 0, no reduction).
##
## @item @qcode{'Standardize'}
## Logical; center and scale each column of @var{X} before embedding (default
## @qcode{false}).
##
## @item @qcode{'Perplexity'}
## The effective number of local neighbors (default 30).  It must be smaller
## than @math{N}.
##
## @item @qcode{'Exaggeration'}
## Tightness factor applied to the high-dimensional affinities for the first 100
## iterations (default 4, no less than 1).
##
## @item @qcode{'LearnRate'}
## The learning rate of the optimization (default 500).
##
## @item @qcode{'InitialY'}
## An @math{N * @var{NumDimensions}} initial embedding (default
## @code{1e-4 * randn}).
##
## @item @qcode{'Options'}
## A structure (as returned by @code{statset}) whose @qcode{MaxIter} (default
## 1000) and @qcode{TolFun} (default @code{1e-10}) fields control the
## optimization.
## @end table
##
## The embedding is not unique: it depends on the initial configuration and the
## random state.  Set @qcode{'InitialY'} (or the random seed) for a reproducible
## result.
##
## @seealso{pca, pdist, statset}
## @end deftypefn

function [Y, loss] = tsne (X, varargin)

  if (nargin < 1)
    print_usage ();
  endif
  if (! (isnumeric (X) && ismatrix (X) && ndims (X) == 2))
    error ("tsne: X must be a numeric matrix.");
  endif
  if (! (isreal (X) && all (! isinf (X(:)))))
    error ("tsne: X must be real and finite.");
  endif

  [N, p] = size (X);

  ## Defaults and Name/Value parsing
  Algorithm = "exact";
  Distance = "euclidean";
  ydims = [];
  numPCA = 0;
  Standardize = false;
  Perplexity = 30;
  Exaggeration = 4;
  LearnRate = 500;
  InitialY = [];
  opts = struct ("MaxIter", 1000, "TolFun", 1e-10);
  if (mod (numel (varargin), 2) != 0)
    error ("tsne: Name/Value arguments must come in pairs.");
  endif
  for k = 1:2:numel (varargin)
    name = varargin{k};
    val = varargin{k+1};
    switch (lower (name))
      case 'algorithm'
        Algorithm = lower (val);
      case 'distance'
        Distance = val;
      case 'numdimensions'
        ydims = val;
      case 'numpcacomponents'
        numPCA = val;
      case 'standardize'
        Standardize = val;
      case 'perplexity'
        Perplexity = val;
      case 'exaggeration'
        Exaggeration = val;
      case 'learnrate'
        LearnRate = val;
      case 'initialy'
        InitialY = val;
      case 'options'
        opts = merge_statset (opts, val);
      otherwise
        error ("tsne: unknown parameter name '%s'.", name);
    endswitch
  endfor

  if (! strcmp (Algorithm, "exact"))
    error ("tsne: only the 'exact' Algorithm is supported.");
  endif
  if (! (isscalar (Perplexity) && Perplexity > 0 && Perplexity < N))
    error ("tsne: 'Perplexity' must be a positive scalar smaller than N.");
  endif
  if (! (isscalar (Exaggeration) && Exaggeration >= 1))
    error ("tsne: 'Exaggeration' must be a scalar not less than 1.");
  endif
  if (! (isscalar (LearnRate) && LearnRate > 0))
    error ("tsne: 'LearnRate' must be a positive scalar.");
  endif

  ## Initial embedding
  if (isempty (InitialY))
    if (isempty (ydims))
      ydims = min (p, 2);
    endif
    Y = 1e-4 * randn (N, ydims);
  else
    if (! isequal (rows (InitialY), N))
      error ("tsne: 'InitialY' must have N rows.");
    endif
    ydims = columns (InitialY);
    Y = InitialY;
  endif

  ## Standardize and optionally reduce with PCA
  if (Standardize)
    sig = std (X, 0, 1);
    sig(range (X, 1) == 0) = 1;
    X = (X - mean (X, 1)) ./ sig;
  endif
  if (numPCA > 0)
    [~, X] = pca (X, "Centered", false, "NumComponents", numPCA);
  endif

  ## High-dimensional affinities P
  D = squareform (pdist (X, Distance)) .^ 2;
  P = binary_search_variance (D, Perplexity, N);
  P(1:N+1:end) = 0;
  P = (P + P') / (2 * N);
  P = max (P, realmin);

  ## Optimize the embedding
  [Y, loss] = tsne_embedding (Y, P, Exaggeration, LearnRate, opts, N, ydims);

endfunction

## ---------------------------------------------------------------------------
## Binary search for the per-point variance giving the target perplexity.
function condP = binary_search_variance (D, perplexity, N)

  condP = zeros (N);
  beta = ones (N, 1);
  H = log (perplexity);
  for i = 1:N
    a = -Inf;
    c = Inf;
    for iter = 1:100
      Pi = exp (-D(i,:) * beta(i));
      Pi(i) = 0;
      si = max (sum (Pi), realmin);
      Pi = Pi / si;
      Hi = log (si) + beta(i) * sum (D(i,:) .* Pi);
      fval = Hi - H;
      if (abs (fval) < 1e-5)
        break;
      endif
      if (fval > 0)
        a = beta(i);
        if (isinf (c))
          beta(i) = 2 * beta(i);
        else
          beta(i) = 0.5 * (beta(i) + c);
        endif
      else
        c = beta(i);
        if (isinf (a))
          beta(i) = 0.5 * beta(i);
        else
          beta(i) = 0.5 * (a + beta(i));
        endif
      endif
    endfor
    condP(i,:) = Pi;
  endfor

endfunction

## ---------------------------------------------------------------------------
## Gradient of the t-SNE cost and the low-dimensional affinities Q.
function [grad, Q] = tsne_gradient (P, Y, N)

  sumY = sum (Y .^ 2, 2);
  num = 1 ./ (1 + sumY + sumY' - 2 * (Y * Y'));
  num(1:N+1:end) = 0;
  Q = max (num ./ sum (num(:)), realmin);
  L = num .* (P - Q);
  grad = 4 * (diag (sum (L, 1)) - L) * Y;

endfunction

## ---------------------------------------------------------------------------
## Gradient descent with adaptive learning rate (Jacobi) and momentum.
function [Y, loss] = tsne_embedding (Y, P, exaggeration, learnrate, opts, N, ydims)

  adp = ones (N, ydims);
  minRate = 0.01;
  momentums = [0.5, 0.8];
  momentumChange = 250;
  exaggerationStop = 100;
  kk = 0.15;
  phi = 0.85;

  P = exaggeration * P;
  Ychange = zeros (N, ydims);
  Q = [];
  for iter = 1:opts.MaxIter
    if (iter == exaggerationStop)
      P = P / exaggeration;
      exaggeration = 1;
    endif

    [grad, Q] = tsne_gradient (P, Y, N);

    ops = sign (grad) != sign (Ychange);
    adp(ops) += kk;
    adp(! ops) *= phi;
    adpRate = learnrate * max (minRate, adp);

    if (iter < momentumChange)
      Ychange = momentums(1) * Ychange - adpRate .* grad;
    else
      Ychange = momentums(2) * Ychange - adpRate .* grad;
    endif
    Y = Y + Ychange;

    if (norm (grad, Inf) < opts.TolFun)
      break;
    endif
  endfor

  if (nargout > 1)
    loss = P(:)' * log (P(:)) - P(:)' * log (Q(:));
  endif

endfunction

## ---------------------------------------------------------------------------
## Copy the recognised statset fields from S into OPTS.
function opts = merge_statset (opts, s)

  if (! isstruct (s))
    error ("tsne: 'Options' must be a structure.");
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
    endswitch
  endfor

endfunction

%!demo
%! ## Embed a small five-dimensional data set into two dimensions.
%! X = [randn(20, 5); randn(20, 5) + 5];
%! Y = tsne (X, "Perplexity", 10);
%! plot (Y(1:20,1), Y(1:20,2), "bo", Y(21:end,1), Y(21:end,2), "rx");
%! title ("t-SNE embedding");

## Verified against MATLAB R2023b (exact algorithm) for
## X = reshape (mod ((1:60)*7, 13), 12, 5) - 6 with a fixed initial embedding.

%!test
%! X = reshape (mod ((1:60)*7, 13), 12, 5) - 6;
%! Y0 = reshape (mod ((1:24)*3, 7), 12, 2) - 3;
%! opt = struct ("MaxIter", 20, "TolFun", 0);
%! [Y, loss] = tsne (X, "Algorithm", "exact", "Distance", "euclidean", ...
%!                   "NumDimensions", 2, "Perplexity", 3, "Exaggeration", 1, ...
%!                   "LearnRate", 100, "Standardize", false, "InitialY", Y0, ...
%!                   "Options", opt);
%! Yref = [-11.087464047262014,  26.776776669616648; ...
%!           4.015876079561814,   4.314538936119215; ...
%!         -10.791976288780955, -28.890539694879799; ...
%!          -6.197553982245900,  11.136429372335288; ...
%!          -4.791764519072724, -29.162528702965172; ...
%!          -0.506807464921708,  13.209994463821561; ...
%!          -1.097171761767460, -24.062824347518902; ...
%!           0.108706666197103,  21.076444386535410; ...
%!           1.270590904095503, -17.709227609503525; ...
%!           7.840024872376860,  19.937127571156172; ...
%!           6.043642805354078, -19.531582453041814; ...
%!          13.753383923495731,  22.678754729129128];
%! assert_equal (Y, Yref, 1e-8);
%! assert_equal (loss, 0.298069770937047, 1e-10);

%!test
%! ## Output sizes and reproducibility with a fixed initial embedding.
%! X = reshape (mod ((1:60)*7, 13), 12, 5) - 6;
%! Y0 = reshape (mod ((1:24)*3, 7), 12, 2) - 3;
%! Y1 = tsne (X, "Perplexity", 3, "InitialY", Y0, ...
%!            "Options", struct ("MaxIter", 50));
%! Y2 = tsne (X, "Perplexity", 3, "InitialY", Y0, ...
%!            "Options", struct ("MaxIter", 50));
%! assert_equal (size (Y1), [12, 2]);
%! assert_equal (Y1, Y2, 1e-12);

%!test
%! ## NumDimensions controls the embedding dimension.
%! X = reshape (mod ((1:60)*7, 13), 12, 5) - 6;
%! Y = tsne (X, "NumDimensions", 3, "Perplexity", 3, ...
%!           "Options", struct ("MaxIter", 20));
%! assert_equal (size (Y), [12, 3]);

## Test input validation
%!error<Invalid call to tsne> tsne ()
%!error<tsne: X must be a numeric matrix.> tsne ({1, 2})
%!error<tsne: X must be real and finite.> tsne ([1, Inf; 2, 3])
%!error<tsne: only the 'exact' Algorithm is supported.> ...
%! tsne (ones (5, 3), "Algorithm", "barneshut")
%!error<tsne: 'Perplexity' must be a positive scalar smaller than N.> ...
%! tsne (ones (5, 3), "Perplexity", 5)
%!error<tsne: 'Exaggeration' must be a scalar not less than 1.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "Exaggeration", 0)
%!error<tsne: 'LearnRate' must be a positive scalar.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "LearnRate", -1)
%!error<tsne: unknown parameter name 'bogus'.> tsne (ones (5, 3), "bogus", 1)
