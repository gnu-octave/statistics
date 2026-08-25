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
## @qcode{'exact'} (default) forms the affinities and the gradient over every
## pair of points, which costs @math{O(N^2)} in time and memory at each
## iteration.  @qcode{'barneshut'} approximates both: the high-dimensional
## affinities are kept only over each point's @math{3 * Perplexity} nearest
## neighbours, and the repulsive part of the gradient is summed over a
## space-partitioning tree of the embedding, giving @math{O(N log N)}.  Use it
## when @math{N} is large enough that the exact algorithm is slow or cannot
## allocate; on this machine the two cost the same at about @math{N = 500} and
## @qcode{'barneshut'} is eight times faster at @math{N = 2000}.
##
## The two do not return the same embedding, and their @var{loss} values are
## not comparable either: the divergence is summed over the pairs that carry an
## affinity, and @qcode{'barneshut'} keeps far fewer of them.
##
## @item @qcode{'Theta'}
## The tree opening criterion for @qcode{'barneshut'}, a non-negative scalar
## (default 0.5).  A cell of the tree is collapsed to its centre of mass when
## its width is smaller than @var{Theta} times its distance from the point
## being pushed, so a larger value is faster and coarser.  @qcode{0} collapses
## nothing and makes the repulsion exact, at @math{O(N^2)}; note that this
## still leaves the affinities sparse, so it does not reproduce
## @qcode{'exact'}.  Ignored by @qcode{'exact'}.
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
  Theta = 0.5;
  opts = statset ("tsne");
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
      case 'theta'
        Theta = val;
      case 'options'
        if (! isstruct (val))
          error ("tsne: 'Options' must be a structure.");
        endif
        opts = statset (opts, val);
      otherwise
        error ("tsne: unknown parameter name '%s'.", name);
    endswitch
  endfor

  if (! any (strcmp (Algorithm, {'exact', 'barneshut'})))
    error ("tsne: 'Algorithm' must be 'exact' or 'barneshut'.");
  endif
  if (! (isscalar (Theta) && isnumeric (Theta) && isreal (Theta)
         && isfinite (Theta) && Theta >= 0))
    error ("tsne: 'Theta' must be a non-negative scalar.");
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
  ## An empty value falls back to the default below.  Anything else has to be a
  ## dimension the embedding can actually have: a non-positive one used to
  ## return an N-by-0 embedding rather than complain, and a fractional or
  ## non-numeric one leaked an internal conversion error.
  if (! isempty (ydims))
    if (! (isnumeric (ydims) && isscalar (ydims) && isreal (ydims)
           && isfinite (ydims) && ydims > 0 && fix (ydims) == ydims))
      error ("tsne: 'NumDimensions' must be a positive integer.");
    elseif (ydims > p)
      error (strcat ("tsne: 'NumDimensions' must not be greater than the", ...
                     " number of columns of X."));
    endif
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

  ## The Barnes-Hut tree is a 2^D-ary subdivision of the embedding, which is
  ## only built for the dimensions t-SNE is used in.
  if (strcmp (Algorithm, "barneshut") && ydims > 3)
    error (strcat ("tsne: 'NumDimensions' must not be greater than 3 for", ...
                   " the 'barneshut' Algorithm."));
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

  ## High-dimensional affinities P.  The exact algorithm forms them for every
  ## pair; Barnes-Hut keeps only each point's nearest neighbours, which is the
  ## half of the approximation that has nothing to do with the tree.
  if (strcmp (Algorithm, "exact"))
    D = squareform (pdist (X, Distance)) .^ 2;
    P = binary_search_variance (D, Perplexity, N);
    P(1:N+1:end) = 0;
    P = (P + P') / (2 * N);
    P = max (P, realmin);
  else
    P = neighbour_affinities (X, Distance, Perplexity, N);
  endif

  ## Optimize the embedding
  [Y, loss] = tsne_embedding (Y, P, Exaggeration, LearnRate, opts, N, ydims, ...
                              Algorithm, Theta);

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
## Sparse high-dimensional affinities over each point's nearest neighbours.
## Barnes-Hut t-SNE keeps 3 * Perplexity of them, the count at which the
## conditional distribution has essentially no mass left, and runs the same
## per-point variance search over that row instead of over all N.
function P = neighbour_affinities (X, Distance, perplexity, N)

  K = min (N - 1, max (1, round (3 * perplexity)));

  ## Neighbours a block of rows at a time.  The exact algorithm forms the whole
  ## N-by-N distance matrix, which is the memory wall this path exists to get
  ## around, so the block is sized to a fixed budget rather than to N.  A block
  ## against all of X is what pdist2 is for; knnsearch would serve too and is
  ## two orders of magnitude slower at this K.
  blk = max (1, min (N, floor (1e7 / N)));
  idx = zeros (N, K);
  D2 = zeros (N, K);
  for lo = 1:blk:N
    hi = min (lo + blk - 1, N);
    Db = pdist2 (X(lo:hi,:), X, Distance) .^ 2;
    ## A point is its own nearest neighbour.  Removing it by index rather than
    ## by distance keeps a duplicated row, which is a neighbour at distance
    ## zero and belongs in the list.
    for r = lo:hi
      Db(r - lo + 1, r) = Inf;
    endfor
    [sv, so] = sort (Db, 2);
    idx(lo:hi,:) = so(:,1:K);
    D2(lo:hi,:) = sv(:,1:K);
  endfor

  ## The same per-point variance search the exact algorithm runs, over the
  ## neighbour row instead of over all N.
  condP = zeros (N, K);
  H = log (perplexity);
  for i = 1:N
    beta = 1;
    a = -Inf;
    c = Inf;
    for iter = 1:100
      Pi = exp (-D2(i,:) * beta);
      si = max (sum (Pi), realmin);
      Pi = Pi / si;
      Hi = log (si) + beta * sum (D2(i,:) .* Pi);
      fval = Hi - H;
      if (abs (fval) < 1e-5)
        break;
      endif
      if (fval > 0)
        a = beta;
        if (isinf (c))
          beta = 2 * beta;
        else
          beta = 0.5 * (beta + c);
        endif
      else
        c = beta;
        if (isinf (a))
          beta = 0.5 * beta;
        else
          beta = 0.5 * (a + beta);
        endif
      endif
    endfor
    condP(i,:) = Pi;
  endfor

  rowi = repmat ((1:N)', 1, K);
  P = sparse (rowi(:), idx(:), condP(:), N, N);
  P = (P + P') / (2 * N);

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
## Gradient of the t-SNE cost by Barnes-Hut summation.  The attraction runs
## over the stored pairs alone and is formed here; the repulsion runs over
## every pair and is summed by the tree in __bhtsne__.  The loss comes with it
## because both its terms are already to hand.
function [grad, loss] = tsne_gradient_bh (P, Y, N, ydims, theta)

  [ii, jj, pv] = find (P);
  dY = Y(ii,:) - Y(jj,:);
  qn = 1 ./ (1 + sum (dY .^ 2, 2));

  Fattr = zeros (N, ydims);
  w = pv .* qn;
  for m = 1:ydims
    Fattr(:,m) = accumarray (ii, w .* dY(:,m), [N, 1]);
  endfor

  [Frep, Z] = __bhtsne__ (Y, theta);

  grad = 4 * (Fattr - Frep / Z);

  if (nargout > 1)
    q = max (qn / Z, realmin);
    loss = sum (pv .* log (max (pv, realmin))) - sum (pv .* log (q));
  endif

endfunction

## ---------------------------------------------------------------------------
## Gradient descent with adaptive learning rate (Jacobi) and momentum.
function [Y, loss] = tsne_embedding (Y, P, exaggeration, learnrate, opts, ...
                                     N, ydims, algorithm, theta)

  adp = ones (N, ydims);
  minRate = 0.01;
  momentums = [0.5, 0.8];
  momentumChange = 250;
  exaggerationStop = 100;
  kk = 0.15;
  phi = 0.85;

  exact = strcmp (algorithm, "exact");

  P = exaggeration * P;
  Ychange = zeros (N, ydims);
  Q = [];
  lossbh = [];
  for iter = 1:opts.MaxIter
    if (iter == exaggerationStop)
      P = P / exaggeration;
      exaggeration = 1;
    endif

    if (exact)
      [grad, Q] = tsne_gradient (P, Y, N);
    else
      [grad, lossbh] = tsne_gradient_bh (P, Y, N, ydims, theta);
    endif

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
    if (exact)
      loss = P(:)' * log (P(:)) - P(:)' * log (Q(:));
    else
      loss = lossbh;
    endif
  endif

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

%!test
%! ## The Barnes-Hut tree reproduces the exact pairwise repulsion when nothing
%! ## is collapsed, which is what pins the summation without an oracle.
%! Y = reshape (mod ((1:60)*7, 13), 30, 2) - 6;
%! [F, Z] = __bhtsne__ (Y, 0);
%! n = rows (Y);
%! Fe = zeros (n, 2);
%! Ze = 0;
%! for i = 1:n
%!   d = Y(i,:) - Y;
%!   q = 1 ./ (1 + sum (d .^ 2, 2));
%!   q(i) = 0;
%!   Ze += sum (q);
%!   Fe(i,:) = sum ((q .^ 2) .* d, 1);
%! endfor
%! assert_equal (F, Fe, 1e-12);
%! assert_equal (Z, Ze, 1e-10);

%!test
%! ## A larger Theta collapses cells, so it approximates rather than reproduces.
%! Y = reshape (mod ((1:60)*7, 13), 30, 2) - 6;
%! [F0, Z0] = __bhtsne__ (Y, 0);
%! [F5, Z5] = __bhtsne__ (Y, 0.5);
%! assert_equal (norm (F5 - F0, "fro") / norm (F0, "fro") < 0.1, true);
%! assert_equal (abs (Z5 - Z0) / Z0 < 0.1, true);
%! assert_equal (isequal (F5, F0), false);

%!test
%! ## The tree handles coincident points, which an embedding starts out full
%! ## of and which a subdivision to one point per cell never separates.
%! Y = zeros (12, 2);
%! [F, Z] = __bhtsne__ (Y, 0.5);
%! assert_equal (F, zeros (12, 2), 1e-12);
%! assert_equal (Z, 132, 1e-10);

%!test
%! ## 'barneshut' returns an embedding of the right shape and a finite loss.
%! X = [reshape(mod((1:100)*7, 13), 20, 5) - 6; ...
%!      reshape(mod((1:100)*7, 13), 20, 5) + 30];
%! [Y, loss] = tsne (X, "Algorithm", "barneshut", "Perplexity", 3, ...
%!                   "NumDimensions", 2, "Options", struct ("MaxIter", 200));
%! assert_equal (size (Y), [40, 2]);
%! assert_equal (all (isfinite (Y(:))), true);
%! assert_equal (isfinite (loss) && loss > 0, true);

%!test
%! ## It optimizes: the divergence falls as the iterations run.
%! X = [reshape(mod((1:100)*7, 13), 20, 5) - 6; ...
%!      reshape(mod((1:100)*7, 13), 20, 5) + 30];
%! args = {"Algorithm", "barneshut", "Perplexity", 3, ...
%!         "InitialY", reshape(mod((1:80)*3, 7), 40, 2) - 3};
%! [~, l1] = tsne (X, args{:}, "Options", struct ("MaxIter", 120));
%! [~, l2] = tsne (X, args{:}, "Options", struct ("MaxIter", 400));
%! assert_equal (l2 < l1, true);

%!test
%! ## It tracks the exact algorithm it approximates, from the same start.
%! X = [reshape(mod((1:100)*7, 13), 20, 5) - 6; ...
%!      reshape(mod((1:100)*7, 13), 20, 5) + 30];
%! o = struct ("MaxIter", 400);
%! args = {"Perplexity", 3, "InitialY", reshape(mod((1:80)*3, 7), 40, 2) - 3, ...
%!         "Options", o};
%! [~, lb] = tsne (X, args{:}, "Algorithm", "barneshut");
%! [~, le] = tsne (X, args{:}, "Algorithm", "exact");
%! assert_equal (abs (lb - le) / le < 0.3, true);

%!test
%! ## Theta reaches the fit: two values give two embeddings.
%! X = reshape (mod ((1:100)*7, 13), 20, 5) - 6;
%! o = struct ("MaxIter", 50);
%! args = {"Algorithm", "barneshut", "Perplexity", 4, "NumDimensions", 2, ...
%!         "InitialY", reshape(mod((1:40)*3, 7), 20, 2) - 3, "Options", o};
%! Ya = tsne (X, args{:}, "Theta", 0);
%! Yb = tsne (X, args{:}, "Theta", 0.8);
%! assert_equal (isequal (Ya, Yb), false);

%!test
%! ## 'barneshut' takes the same NumDimensions the tree is built for.
%! X = reshape (mod ((1:100)*7, 13), 20, 5) - 6;
%! for d = 1:3
%!   Y = tsne (X, "Algorithm", "barneshut", "NumDimensions", d, ...
%!             "Perplexity", 4, "Options", struct ("MaxIter", 20));
%!   assert_equal (size (Y), [20, d]);
%! endfor

## Test input validation
%!error<Invalid call to tsne> tsne ()
%!error<tsne: X must be a numeric matrix.> tsne ({1, 2})
%!error<tsne: X must be real and finite.> tsne ([1, Inf; 2, 3])
%!error<tsne: 'Algorithm' must be 'exact' or 'barneshut'.> ...
%! tsne (ones (5, 3), "Algorithm", "bogus")
%!error<tsne: 'Theta' must be a non-negative scalar.> ...
%! tsne (ones (5, 3), "Algorithm", "barneshut", "Theta", -1)
%!error<tsne: 'Theta' must be a non-negative scalar.> ...
%! tsne (ones (5, 3), "Algorithm", "barneshut", "Theta", [1, 2])
%!error<tsne: 'NumDimensions' must not be greater than 3 for the 'barneshut' Algorithm.> ...
%! tsne (ones (8, 5), "Algorithm", "barneshut", "NumDimensions", 4, "Perplexity", 2)
%!error<tsne: 'Perplexity' must be a positive scalar smaller than N.> ...
%! tsne (ones (5, 3), "Perplexity", 5)
%!error<tsne: 'Exaggeration' must be a scalar not less than 1.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "Exaggeration", 0)
%!error<tsne: 'LearnRate' must be a positive scalar.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "LearnRate", -1)
%!error<tsne: unknown parameter name 'bogus'.> tsne (ones (5, 3), "bogus", 1)
%!error<tsne: 'NumDimensions' must be a positive integer.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "NumDimensions", 0)
%!error<tsne: 'NumDimensions' must be a positive integer.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "NumDimensions", -1)
%!error<tsne: 'NumDimensions' must be a positive integer.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "NumDimensions", 2.5)
%!error<tsne: 'NumDimensions' must be a positive integer.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "NumDimensions", [2, 3])
%!error<tsne: 'NumDimensions' must be a positive integer.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "NumDimensions", "a")
%!error<tsne: 'NumDimensions' must be a positive integer.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "NumDimensions", true)
%!error<tsne: 'NumDimensions' must be a positive integer.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "NumDimensions", Inf)
%!error<tsne: 'NumDimensions' must not be greater than the number of columns of X.> ...
%! tsne (ones (5, 3), "Perplexity", 2, "NumDimensions", 4)

## An empty NumDimensions falls back to the default, and an integer-typed one
## is accepted, both as in MATLAB.
%!test
%! X = [1 2 3; 2 4 6; 3 5 9; 4 8 11; 5 9 14; 6 11 17];
%! Y = tsne (X, "Perplexity", 2, "NumDimensions", []);
%! assert_equal (columns (Y), 2);
%! Y = tsne (X, "Perplexity", 2, "NumDimensions", int32 (2));
%! assert_equal (columns (Y), 2);
%! Y = tsne (X, "Perplexity", 2, "NumDimensions", 3);
%! assert_equal (columns (Y), 3);
