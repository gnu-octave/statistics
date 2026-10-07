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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{knnstat} =} knntest (@var{X}, @var{Y})
## @deftypefnx {statistics} {@var{knnstat} =} knntest (@var{X}, @var{Y}, @var{Name}, @var{Value})
## @deftypefnx {statistics} {[@var{knnstat}, @var{p}] =} knntest (@dots{})
## @deftypefnx {statistics} {[@var{knnstat}, @var{p}, @var{h}] =} knntest (@dots{})
##
## Two-sample multivariate test based on nearest neighbours.
##
## @code{@var{knnstat} = knntest (@var{X}, @var{Y})} returns the nearest
## neighbour statistic of the two samples @var{X} and @var{Y}, whose rows are
## observations and whose columns are variables.  The two samples are pooled,
## each observation's @math{k} nearest neighbours are found among the others,
## and @var{knnstat} is the share of those neighbours that belong to the same
## sample as the observation.  It lies in @math{[0, 1]}: near 1 the samples
## are well separated, and near 0.5 two samples of equal size are alike.
##
## @var{X} and @var{Y} are both numeric matrices with the same number of
## columns, or both tables.  Two tables are read over the variables they share,
## which must be all the variables of one of them.  An observation holding a
## missing value in a variable used is left out.
##
## @code{[@var{knnstat}, @var{p}, @var{h}] = knntest (@dots{})} also returns
## the p-value of the right-tailed test of the null hypothesis that @var{X} and
## @var{Y} come from the same distribution, and @var{h}, which is 1 where the
## null hypothesis is rejected at the significance level @qcode{'Alpha'} and 0
## otherwise.  The p-value takes @math{(knnstat - mu) / sigma} as standard
## normal, with @math{m_x} and @math{m_y} the sizes of the samples, @math{m}
## their sum, and @math{q = m_x m_y / m^2}:
## @tex
## $$\mu = {m_x (m_x - 1) + m_y (m_y - 1) \over m (m - 1)}, \qquad
## \sigma^2 = {q + 4 q^2 \over m k}.$$
## @end tex
## @ifnottex
## @math{mu = (m_x (m_x - 1) + m_y (m_y - 1)) / (m (m - 1))} and
## @math{sigma^2 = (q + 4 q^2) / (m k)}.
## @end ifnottex
##
## @code{@var{knnstat} = knntest (@dots{}, @var{Name}, @var{Value})} takes the
## following options.
##
## @multitable @columnfractions 0.26 0.02 0.72
## @headitem @var{Name} @tab @tab @var{Value}
## @item @qcode{'Alpha'} @tab @tab The significance level, a scalar between 0
## and 1, 0.05 by default.
## @item @qcode{'NumNeighbors'} @tab @tab The number of nearest neighbours
## @math{k}, a positive integer, 10 by default.  It is taken as one less than
## the number of observations where it exceeds that.
## @item @qcode{'Distance'} @tab @tab The distance between observations; see
## below.
## @item @qcode{'VariableNames'} @tab @tab The variables to use, among those
## the two tables share, as a character vector, a string array or a cell array
## of character vectors.  All the shared variables by default.
## @item @qcode{'CategoricalVariables'} @tab @tab The variables holding
## levels: @qcode{'all'}, their indices, a logical vector over the variables,
## or, for tables, their names.  A table variable holding logical values, an
## unordered @code{categorical} array, a @code{string} array or a cell array of
## character vectors holds levels whether named here or not.  A matrix holds
## none unless named.
## @end multitable
##
## @qcode{'Distance'} is one of the following, @qcode{'seuclidean'} by default
## where every variable is continuous, @qcode{'hamming'} where every variable
## holds levels, and @qcode{'goodall3'} where the two are mixed.
##
## @multitable @columnfractions 0.22 0.02 0.76
## @headitem @var{Distance} @tab @tab Description
## @item @qcode{'euclidean'} @tab @tab Euclidean distance.
## @item @qcode{'seuclidean'} @tab @tab Euclidean distance with each variable
## divided by its standard deviation over @code{[@var{X}; @var{Y}]}.
## @item @qcode{'cityblock'} @tab @tab City block distance.
## @item @qcode{'cosine'} @tab @tab One minus the cosine of the angle between
## two observations.
## @item @qcode{'correlation'} @tab @tab One minus the correlation between two
## observations.
## @item @qcode{'fasteuclidean'}, @qcode{'fastseuclidean'} @tab @tab The same
## distances as @qcode{'euclidean'} and @qcode{'seuclidean'}, computed exactly.
## @item @qcode{'hamming'} @tab @tab The share of variables that differ.
## @item @qcode{'goodall3'} @tab @tab The Goodall 3 dissimilarity of
## @code{nomdist}, with level frequencies counted over
## @code{[@var{X}; @var{Y}]} and the observation whose neighbours are sought
## counted once more, as @code{nomdist2} counts a query with
## @qcode{'CountQuery'}.  A continuous variable adds its absolute difference
## divided by its range over @code{[@var{X}; @var{Y}]}, and the sum over
## every variable is divided by their number.
## @end multitable
##
## The first five treat a variable holding levels as numbers, its levels coded
## @code{1}, @code{2}, @dots{}; @qcode{'hamming'} treats every distinct value
## of a continuous variable as a level.
##
## A tie between two neighbours at the same distance is broken by the order of
## the pooled sample, the rows of @var{X} coming before those of @var{Y}, and
## distances equal to within rounding count as tied, so the result does not
## depend on the platform.  MATLAB lets rounding decide such ties, so on data
## where many distances are equal, such as values on a coarse grid,
## @var{knnstat} can differ slightly from MATLAB's.
##
## MATLAB defines its @qcode{'goodall3'} for a continuous variable in a way it
## does not document, and which none of the constructions tried here
## reproduces, so on mixed data @var{knnstat} differs from MATLAB's.  On
## levels alone the two agree.
##
## Reference: M. Williams (2010).  How good are your fits? Unbinned
## multivariate goodness-of-fit tests in high energy physics.  Journal of
## Instrumentation, 5(09), P09004.
##
## @seealso{nomdist, nomdist2, kstest2}
## @end deftypefn

function [knnstat, p, h] = knntest (X, Y, varargin)

  ## Input validation
  if (nargin < 2)
    error ("knntest: too few input arguments.");
  endif
  optNames = {'Alpha', 'NumNeighbors', 'Distance', 'VariableNames', ...
              'CategoricalVariables'};
  ## An empty Distance resolves to the default for the variables' kinds
  dfValues = {0.05, 10, [], [], []};
  [alpha, k, dist, vnames, catvars, args] = ...
                        parsePairedArguments (optNames, dfValues, varargin(:));
  if (! isempty (args))
    error ("knntest: invalid optional paired argument.");
  endif
  if (! (isnumeric (alpha) && isreal (alpha) && isscalar (alpha)
         && alpha > 0 && alpha < 1))
    error ("knntest: 'Alpha' must be a scalar between 0 and 1.");
  endif
  if (! (isnumeric (k) && isreal (k) && isscalar (k) && isfinite (k)
         && k >= 1 && k == fix (k)))
    error ("knntest: 'NumNeighbors' must be a positive integer.");
  endif

  ## The pooled sample coded one column per variable, rows holding a missing
  ## value left out, and which variables hold levels
  [Z, iscat, mx, my, errmsg] = __twosample__ (X, Y, vnames, catvars);
  if (! isempty (errmsg))
    error ("knntest: %s", errmsg);
  endif
  m = mx + my;

  ## The distance, whose default depends on the kinds of variable
  if (isempty (dist))
    if (all (iscat))
      dist = 'hamming';
    elseif (any (iscat))
      dist = 'goodall3';
    else
      dist = 'seuclidean';
    endif
  endif
  dists = {'cityblock', 'correlation', 'cosine', 'euclidean', ...
           'fasteuclidean', 'seuclidean', 'fastseuclidean', 'hamming', ...
           'goodall3'};
  if (isa (dist, 'string') && isscalar (dist))
    dist = char (dist);
  endif
  if (! (ischar (dist) && isrow (dist) && any (strcmpi (dist, dists))))
    error ("knntest: 'Distance' must be one of %s.", ...
           strjoin (strcat ("'", dists, "'"), ', '));
  endif
  dist = lower (dist);

  g = [zeros(mx, 1); ones(my, 1)];

  ## The same-sample share among each observation's k nearest neighbours
  k = min (k, m - 1);
  I = zeros (k, m);
  step = max (1, min (m, floor (2e6 / m)));
  for b1 = 1:step:m
    b = b1:min (b1 + step - 1, m);
    D = knnDistances (Z, b, dist, iscat);
    D(sub2ind (size (D), b, 1:numel (b))) = Inf;
    ## Distances equal to within rounding tie, and sort keeps the lower row
    s = max (abs (D(isfinite (D))));
    if (! isempty (s) && s > 0)
      D = round (D / s * 1e12);
    endif
    [~, o] = sort (D, 1);
    I(:,b) = reshape (g(o(1:k,:)), k, numel (b)) == g(b)';
  endfor
  knnstat = mean (I(:));

  ## Right tail of the normal approximation
  mu = (mx * (mx - 1) + my * (my - 1)) / (m * (m - 1));
  q = mx * my / m ^ 2;
  sigma = sqrt ((q + 4 * q ^ 2) / (m * k));
  p = 0.5 * erfc ((knnstat - mu) / sigma / sqrt (2));
  h = double (p <= alpha);

endfunction

## Distances from every observation to the observations B, one column each.
function D = knnDistances (Z, b, dist, iscat)
  switch (dist)
    case {'euclidean', 'fasteuclidean'}
      D = pdist2 (Z, Z(b,:), 'euclidean');
    case {'seuclidean', 'fastseuclidean'}
      D = pdist2 (Z, Z(b,:), 'seuclidean', std (Z, [], 1));
    case {'cityblock', 'cosine', 'correlation'}
      D = pdist2 (Z, Z(b,:), dist);
    case 'hamming'
      D = pdist2 (Z, Z(b,:), 'hamming');
    case 'goodall3'
      K = columns (Z);
      D = zeros (rows (Z), numel (b));
      if (any (iscat))
        C = Z(:,iscat);
        D = columns (C) * nomdist2 (C, C(b,:), 'goodall3', 'CountQuery', true);
      endif
      if (! all (iscat))
        R = range (Z(:,! iscat), 1);
        R(R == 0) = 1;
        W = Z(:,! iscat) ./ R;
        D += pdist2 (W, W(b,:), 'cityblock');
      endif
      D /= K;
  endswitch
endfunction

%!shared XS, YS, XA, YA
%! XS = [8, 4; 13, 1; 14, 2; 3, 3];
%! YS = [3, 19; 1, 7; 1, 17; 4, 10; 2, 17];
%! XA = [2, 2, 3; 3, 1, 1; 3, 1, 3; 3, 3, 1; 3, 3, 2; 1, 1, 1];
%! YA = [3, 2, 2; 1, 2, 2; 3, 2, 2; 2, 3, 3; 2, 1, 1; 3, 1, 1];

## Expected values from MATLAB R2026a unless stated otherwise
%!test
%! ## A tie goes to the lower pooled index
%! [s, p, h] = knntest ([0; 20], [1; 2], 'NumNeighbors', 1, ...
%!                      'Distance', 'euclidean');
%! assert_equal (s, 0.25);
%! assert_equal (p, 0.593168142116604, -1e-13);
%! assert_equal (h, 0);
%!test
%! [s, p] = knntest ([1; 2], [0; 20], 'NumNeighbors', 1, ...
%!                   'Distance', 'euclidean');
%! assert_equal (s, 0.5);
%! assert_equal (p, 0.318675944116969, -1e-13);
%!test
%! ## An observation is never its own neighbour, even beside a duplicate
%! [s, p] = knntest ([0; 5], [0; 9], 'NumNeighbors', 1, ...
%!                   'Distance', 'euclidean');
%! assert_equal (s, 0);
%! assert_equal (p, 0.827110706924420, -1e-13);
%!test
%! ## seuclidean, the default, scales by the standard deviation over [X; Y]
%! [s, p, h] = knntest (XS, YS, 'NumNeighbors', 1);
%! assert_equal (s, 7 / 9, -1e-14);
%! assert_equal (p, 0.076726923031103, -1e-13);
%! assert_equal (h, 0);
%!assert_equal (knntest (XS, YS, 'NumNeighbors', 1, ...
%!                       'Distance', 'SEuclidean'), ...
%!              knntest (XS, YS, 'NumNeighbors', 1))
%!assert_equal (knntest (XS, YS, 'NumNeighbors', 1, 'Distance', ...
%!                       'fastseuclidean'), knntest (XS, YS, 'NumNeighbors', 1))
%!assert_equal (knntest (XS, YS, 'NumNeighbors', 1, ...
%!                       'Distance', 'fasteuclidean'), ...
%!              knntest (XS, YS, 'NumNeighbors', 1, 'Distance', 'euclidean'))
%!assert_equal (knntest ([XS; NaN, 2], YS, 'NumNeighbors', 1), 7 / 9, -1e-14)
%!test
%! [~, ~, h] = knntest (XS, YS, 'NumNeighbors', 1, 'Alpha', 0.5);
%! assert_equal (h, 1);
%!test
%! ## Goodall 3 counts the query into the frequencies over [X; Y]
%! [s, p] = knntest (XA, YA, 'NumNeighbors', 1, 'CategoricalVariables', ...
%!                   'all', 'Distance', 'goodall3');
%! assert_equal (s, 7 / 12, -1e-14);
%! assert_equal (p, 0.264043416942331, -1e-13);
%!test
%! ## hamming is the default where every variable holds levels
%! [s, p] = knntest (XA, YA, 'NumNeighbors', 1, 'CategoricalVariables', 'all');
%! assert_equal (s, 7 / 12, -1e-14);
%! assert_equal (p, 0.264043416942331, -1e-13);
%!test
%! ## A level of Y that X lacks keeps its row
%! [s, p] = knntest (XA, [YA; 4, 1, 1], 'NumNeighbors', 1, ...
%!                   'CategoricalVariables', 'all', 'Distance', 'goodall3');
%! assert_equal (s, 7 / 13, -1e-14);
%! assert_equal (p, 0.346797481676874, -1e-13);
%!assert_equal (knntest (XA, [YA; 4, 1, 1], 'NumNeighbors', 1, ...
%!                      'CategoricalVariables', 'all'), 7 / 13, -1e-14)
%!test
%! X = [4, 1, 4, 4, 4, 1; 2, 1, 2, 2, 2, 3; 1, 3, 1, 2, 2, 3; 2, 2, 1, 4, 4, 2;
%!      2, 3, 3, 1, 3, 4; 2, 4, 1, 2, 2, 3; 3, 1, 2, 3, 4, 1];
%! Y = [1, 4, 4, 2, 4, 1; 1, 3, 1, 1, 4, 2; 1, 3, 3, 1, 3, 2; 1, 3, 2, 4, 2, 3;
%!      2, 3, 4, 1, 2, 3; 1, 4, 2, 3, 4, 2; 2, 3, 2, 3, 2, 1];
%! [s, p] = knntest (X, Y, 'NumNeighbors', 1, 'CategoricalVariables', 'all', ...
%!                   'Distance', 'goodall3');
%! assert_equal (s, 5 / 14, -1e-14);
%! assert_equal (p, 0.709666127521303, -1e-13);
%!test
%! X = [1, 2, 1, 2, 1, 4; 3, 1, 4, 1, 4, 4; 4, 1, 3, 2, 2, 3; 3, 3, 4, 2, 4, 3;
%!      4, 3, 2, 4, 3, 4; 2, 3, 4, 3, 1, 1; 3, 4, 1, 2, 3, 2];
%! Y = [4, 2, 4, 2, 1, 2; 4, 1, 1, 3, 4, 1; 4, 3, 1, 3, 4, 4; 4, 4, 2, 4, 2, 1;
%!      2, 2, 2, 3, 3, 1; 1, 3, 4, 4, 4, 3; 2, 2, 4, 3, 3, 3];
%! s = knntest (X, Y, 'NumNeighbors', 1, 'CategoricalVariables', 'all', ...
%!              'Distance', 'goodall3');
%! assert_equal (s, 5 / 14, -1e-14);
%!test
%! X = [3, 1, 2, 2; 1, 1, 3, 2; 2, 1, 3, 1; 3, 3, 1, 3; 3, 2, 2, 1;
%!      2, 3, 2, 1; 1, 1, 2, 2];
%! Y = [3, 2, 3, 3; 2, 3, 1, 3; 2, 3, 3, 3; 2, 2, 3, 2; 3, 1, 1, 3;
%!      3, 2, 1, 3; 3, 1, 1, 2];
%! [s, p] = knntest (X, Y, 'NumNeighbors', 1, 'CategoricalVariables', 'all', ...
%!                   'Distance', 'goodall3');
%! assert_equal (s, 5 / 7, -1e-14);
%! assert_equal (p, 0.090543972339391, -1e-13);
%!test
%! load fisheriris
%! [s, p, h] = knntest (meas(1:25,:), meas(26:50,:));
%! assert_equal (s, 0.438, -1e-14);
%! assert_equal (p, 0.949281930319421, -1e-13);
%! assert_equal (h, 0);
%!test
%! ## Ties decided by rounding in MATLAB, which gives 0.909
%! load fisheriris
%! [s, p, h] = knntest (meas(51:100,:), meas(101:150,:));
%! assert_equal (s, 0.908, -1e-14);
%! assert_equal (p, 1.72915613588074e-76, -1e-10);
%! assert_equal (h, 1);
%!test
%! ## Ties decided by rounding in MATLAB, which gives 0.92
%! load fisheriris
%! s = knntest (meas(51:100,:), meas(101:150,:), 'NumNeighbors', 5, ...
%!              'Distance', 'cityblock');
%! assert_equal (s, 0.916, -1e-14);
%!test
%! ## NumNeighbors above the observations but one is taken as that
%! [s, p] = knntest ([0; 20], [1; 2]);
%! assert_equal (s, 1 / 3, -1e-14);
%! assert_equal (p, 0.5, -1e-14);
%!test
%! ## Mixed data: ours; MATLAB's undocumented rule gives 0.583333
%! X = [3, 51; 3, 53; 2, 44; 2, 56; 1, 46; 3, 45];
%! Y = [1, 69; 3, 8; 2, 23; 2, 1; 1, 5; 1, 86];
%! [s, p] = knntest (X, Y, 'NumNeighbors', 1, 'CategoricalVariables', 1);
%! assert_equal (s, 2 / 3, -1e-14);
%! assert_equal (p, 0.14936110409448, -1e-12);
%!test
%! X = [3, 51; 3, 53; 2, 44; 2, 56; 1, 46; 3, 45];
%! Y = [1, 69; 3, 8; 2, 23; 2, 1; 1, 5; 1, 86];
%! s = knntest (X, Y, 'NumNeighbors', 1, 'CategoricalVariables', [true, false]);
%! assert_equal (s, knntest (X, Y, 'NumNeighbors', 1, ...
%!                          'CategoricalVariables', 1));
%!test
%! ## A table's categorical variable holds levels without being named
%! X = [1, 89; 2, 58; 1, 99; 2, 52; 1, 84; 3, 49];
%! Y = [1, 72; 2, 58; 3, 100; 3, 97; 1, 68; 2, 70];
%! TX = table (categorical (X(:,1)), X(:,2));
%! TY = table (categorical (Y(:,1)), Y(:,2));
%! assert_equal (knntest (TX, TY, 'NumNeighbors', 1), ...
%!               knntest (X, Y, 'NumNeighbors', 1, 'CategoricalVariables', 1));
%!test
%! C = {'a', 'b', 'c'};
%! TX = table (C(XA(:,1))', XA(:,2) == 1, string (C(XA(:,3)))', ...
%!             'VariableNames', {'A', 'B', 'C'});
%! TY = table (C(YA(:,1))', YA(:,2) == 1, string (C(YA(:,3)))', ...
%!             'VariableNames', {'A', 'B', 'C'});
%! assert_equal (knntest (TX, TY, 'NumNeighbors', 1), ...
%!               knntest ([XA(:,1), XA(:,2) == 1, XA(:,3)], ...
%!                        [YA(:,1), YA(:,2) == 1, YA(:,3)], ...
%!                        'NumNeighbors', 1, 'CategoricalVariables', 'all'));
%!test
%! ## Tables are read over the variables they share
%! TX = table (XS(:,1), XS(:,2), (1:4)', 'VariableNames', {'A', 'B', 'Z'});
%! TY = table (YS(:,2), YS(:,1), 'VariableNames', {'B', 'A'});
%! assert_equal (knntest (TX, TY, 'NumNeighbors', 1), 7 / 9, -1e-14);
%!test
%! TX = table (XS(:,1), XS(:,2), 'VariableNames', {'A', 'B'});
%! TY = table (YS(:,1), YS(:,2), 'VariableNames', {'A', 'B'});
%! assert_equal (knntest (TX, TY, 'NumNeighbors', 1, 'VariableNames', 'A'), ...
%!               knntest (XS(:,1), YS(:,1), 'NumNeighbors', 1));

%!error<knntest: too few input arguments.> knntest (1)
%!error<knntest: invalid optional paired argument.> knntest (XS, YS, 'Tail', 1)
%!error<knntest: 'Alpha' must be a scalar between 0 and 1.> ...
%! knntest (XS, YS, 'Alpha', 1)
%!error<knntest: 'NumNeighbors' must be a positive integer.> ...
%! knntest (XS, YS, 'NumNeighbors', 1.5)
%!error<knntest: X and Y must both be matrices or both be tables.> ...
%! knntest (table (XS), YS)
%!error<knntest: 'VariableNames' applies only where X and Y are tables.> ...
%! knntest (XS, YS, 'VariableNames', 'A')
%!error<knntest: X and Y must be real numeric matrices or tables.> ...
%! knntest ({1}, YS)
%!error<knntest: X and Y must have the same number of columns.> ...
%! knntest (XS, YS(:,1))
%!error<knntest: the variable names of one of X and Y must all be names of the other.> ...
%! knntest (table ([1; 2], 'VariableNames', {'A'}), ...
%!          table ([1; 2], 'VariableNames', {'B'}))
%!error<knntest: 'VariableNames' must name variables X and Y share.> ...
%! knntest (table ([1; 2], 'VariableNames', {'A'}), ...
%!          table ([1; 2], 'VariableNames', {'A'}), 'VariableNames', 'B')
%!error<knntest: variable 'A' must be one column.> ...
%! knntest (table ([1, 2; 3, 4], 'VariableNames', {'A'}), ...
%!          table ([1, 2; 3, 4], 'VariableNames', {'A'}))
%!error<knntest: variable 'A' must hold numbers, logical values, categories or text.> ...
%! knntest (table ({1; 2}, 'VariableNames', {'A'}), ...
%!          table ({1; 2}, 'VariableNames', {'A'}))
%!error<knntest: variable 'A' must hold the same kind of values in X and Y.> ...
%! knntest (table ([1; 2], 'VariableNames', {'A'}), ...
%!          table ({'a'; 'b'}, 'VariableNames', {'A'}))
%!error<knntest: 'CategoricalVariables' can name variables only where X and Y are tables.> ...
%! knntest (XS, YS, 'CategoricalVariables', 'A')
%!error<knntest: 'CategoricalVariables' must name variables in use.> ...
%! knntest (table ([1; 2], 'VariableNames', {'A'}), ...
%!          table ([1; 2], 'VariableNames', {'A'}), 'CategoricalVariables', 'B')
%!error<knntest: 'CategoricalVariables' must hold one logical value for each of the 2 variables.> ...
%! knntest (XS, YS, 'CategoricalVariables', true)
%!error<knntest: 'CategoricalVariables' must be 'all', indices from 1 to 2, a logical vector or variable names.> ...
%! knntest (XS, YS, 'CategoricalVariables', 3)
%!error<knntest: X and Y must each hold an observation with no missing value.> ...
%! knntest ([NaN, 1], YS)
%!error<knntest: 'Distance' must be one of 'cityblock', 'correlation', 'cosine', 'euclidean', 'fasteuclidean', 'seuclidean', 'fastseuclidean', 'hamming', 'goodall3'.> ...
%! knntest (XS, YS, 'Distance', 'minkowski')
