## Copyright (C) 2014-2019 Piotr Dollar <pdollar@gmail.com>
## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation; either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see
## <http:##www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{d} =} pdist2 (@var{x}, @var{y})
## @deftypefnx {statistics} {@var{d} =} pdist2 (@var{x}, @var{y}, @var{metric})
##
## Compute pairwise distance between two sets of vectors.
##
## Let @var{X} be an MxP matrix representing m points in P-dimensional space
## and @var{Y} be an NxP matrix representing another set of points in the same
## space.  This function computes the M-by-N distance matrix @var{D} where
## @code{@var{D}(i,j)} is the distance between @code{@var{X}(i,:)} and
## @code{@var{Y}(j,:)}.
##
## The optional argument @var{metric} can be used to select different
## distances:
##
## @table @asis
## @item @qcode{"euclidean"} (default)
##
## @item @qcode{"sqeuclidean"}
## Compute the squared euclidean distance, i.e., the euclidean distance
## before computing square root.  This is ideal when the interest is on the
## order of the euclidean distances rather than the actual distance value
## because it performs significantly faster while preserving the order.
##
## @item @qcode{"chisq'"}
## The chi-squared distance between two vectors is defined as:
## @code{d(x, y) = sum ((xi-yi)^2 / (xi+yi)) / 2}.
## The chi-squared distance is useful when comparing histograms.
##
## @item @qcode{"cosine"}
## Distance is defined as the cosine of the angle between two vectors.
##
## @item @qcode{"emd"}
## Earth Mover's Distance (EMD) between positive vectors (histograms).
## Note for 1D, with all histograms having equal weight, there is a simple
## closed form for the calculation of the EMD.  The EMD between histograms
## @var{x} and @var{y} is given by @code{sum (abs (cdf (x) - cdf (y)))},
## where @code{cdf} is the cumulative distribution function (computed
## simply by @code{cumsum}).
##
## @item @qcode{"L1"}
## The L1 distance between two vectors is defined as:  @code{sum (abs (x-y))}
##
## @end table
##
## @seealso{pdist}
## @end deftypefn

function [D, I] = pdist2 (X, Y, varargin)

  ## Check input data
  if (nargin < 2)
	  error ("pdist2: too few input arguments.");
  endif
  if (size (X, 2) != size (Y, 2))
	  error ("pdist2: X and Y must have equal number of columns.");
  endif
  if (ndims (X) != 2 || ndims (Y) != 2)
    error ("pdist2: X and Y must be 2 dimensional matrices.");
  endif

  ## Add default values
  Distance = "euclidean";   # Distance metric
  DistParameter = [];       # Distance parameter
  SortOrder = [];           # Flag for sorting distances to find

  ## Parse additional Distance metric and Distance parameter (if available)
  DMs = {"euclidean", "squaredeuclidean", "seuclidean", ...
         "mahalanobis", "cityblock", "minkowski", "chebychev", ...
         "cosine", "correlation", "hamming", "jaccard", "spearman"};
  if (numel (varargin) > 0)
    if (any (strcmpi (DMs, varargin{1})))
      Distance = tolower (varargin{1});
      varargin(1) = [];
      if (numel (varargin) > 0)
        if (isnumeric (varargin{1}))
          DistParameter = varargin{1};
          varargin(1) = [];
        endif
      endif
    endif
  endif

  ## Parse additional parameters in Name/Value pairs
  parcount = 0;
  while (numel (varargin) > 0)
    if (numel (varargin) < 2)
      error ("pdist2: missing value in optional name/value paired arguments.");
    endif
    switch (tolower (varargin{1}))
      case "smallest"
        SortOrder = "ascend";
        K = varargin{2};
        parcount += 1;
      case "largest"
        SortOrder = "descend";
        K = varargin{2};
        parcount += 1;
      otherwise
        error ("pdist2: invalid NAME in optional pairs of arguments.");
    endswitch
    varargin(1:2) = [];
  endwhile

  ## Check additional arguments
  if (parcount > 1)
    error ("pdist2: you can only use either Smallest or Largest.");
  endif
  if (isempty (SortOrder) && nargout > 1)
    error (strcat (["pdist2: Smallest or Largest must be specified"], ...
                   [" to compute second output."]));
  endif

  ## Calculate selected distance
  [ix, iy] = meshgrid (1:size (X, 1), 1:size (Y, 1));
  switch (Distance)
  case "euclidean"
    D = sqrt (sum ((X(ix(:),:) - Y(iy(:),:)) .^ 2, 2));

  case "squaredeuclidean"
    D = sum ((X(ix(:),:) - Y(iy(:),:)) .^ 2, 2);

  case "seuclidean"
    if (isempty (DistParameter))
      DistParameter = std (X, [], 1);
    else
      if (numel (DistParameter) != columns (X))
        error (strcat (["pdist2: DistParameter for standardized euclidean"], ...
                       [" must be a vector of equal length to the number"], ...
                       [" of columns in X."]));
      endif
      if (any (DistParameter < 0))
        error (strcat (["pdist2: DistParameter for standardized euclidean"], ...
                       [" must be a nonnegative vector."]));
      endif
    endif
    DistParameter(DistParameter == 0) = 1;  # fix constant variable
    D = sqrt (sum (((X(ix(:),:) - Y(iy(:),:)) ./ DistParameter) .^ 2, 2));

  case "mahalanobis"
    if (isempty (DistParameter))
      DistParameter = cov (X(! any (isnan (X), 2),:));
    else
      if (columns (DistParameter) != columns (X))
        error (strcat (["pdist2: DistParameter for mahalanobis distance"], ...
                       [" must be a covariance matrix with the same"], ...
                       [" number of columns as X."]));
      endif
      try
        chol (DistParameter);
      catch ME
        error (strcat (["pdist2: covariance matrix for mahalanobis"],...
                       [" distance must be symmetric and positive definite."]));
      end_try_catch
    endif
    dxy = X(ix(:),:) - Y(iy(:),:);
    D   = sqrt (sum ((dxy  * inv (DistParameter)) .* dxy, 2));

  case "cityblock"
    D = sum (abs (X(ix(:),:) - Y(iy(:),:)), 2);

  case "minkowski"
    if (isempty (DistParameter))
      DistParameter = 2;
    else
      if (! (isnumeric (DistParameter) && isscalar (DistParameter)
                                       && DistParameter > 0))
        error (strcat (["pdist2: DistParameter for minkowski distance"],...
                       [" must be a positive scalar."]));
      endif
    endif
    D = sum (abs (X(ix(:),:) - Y(iy(:),:)) .^ DistParameter, 2) .^ ...
            (1 / DistParameter);

  case "chebychev"
    D = max (abs (X(ix(:),:) - Y(iy(:),:)), [], 2);

  case "cosine"
    sx = sum (X .^ 2, 2) .^ (-1 / 2);
    sy = sum (Y .^ 2, 2) .^ (-1 / 2);
    D  = 1 - sum (X(ix(:),:) .* Y(iy(:),:), 2) .* sx(ix(:)) .* sy(iy(:));

  case "correlation"
    mX = mean (X(ix(:),:), 2);
    mY = mean (Y(iy(:),:), 2);
    xy = sum ((X(ix(:),:) - mX) .* (Y(iy(:),:) - mY), 2);
    xx = sqrt (sum ((X(ix(:),:) - mX) .* (X(ix(:),:) - mX), 2));
    yy = sqrt (sum ((Y(iy(:),:) - mY) .* (Y(iy(:),:) - mY), 2));
    D = 1 - (xy ./ (xx .* yy));

  case "hamming"
    D = mean (abs (X(ix(:),:) != Y(iy(:),:)), 2);

  case "jaccard"
    xy0 = (X(ix(:),:) != 0 | Y(iy(:),:) != 0);
    D = sum ((X(ix(:),:) != Y(iy(:),:)) & xy0, 2) ./ sum (xy0, 2);

  case "spearman"
    for i = 1:size (X, 1)
      rX(i,:) = tiedrank (X(i,:));
    endfor
    for i = 1:size (Y, 1)
      rY(i,:) = tiedrank (Y(i,:));
    endfor
    rM = (size (X, 2) + 1) / 2;
    xy = sum ((rX(ix(:),:) - rM) .* (rY(iy(:),:) - rM), 2);
    xx = sqrt (sum ((rX(ix(:),:) - rM) .* (rX(ix(:),:) - rM), 2));
    yy = sqrt (sum ((rY(iy(:),:) - rM) .* (rY(iy(:),:) - rM), 2));
    D = 1 - (xy ./ (xx .* yy));

  endswitch

  D = reshape (D, size (Y, 1), size (X, 1))';

  if (nargout > 1)
    [D, I] = sort (D', 2, SortOrder);
    K = min (size (D, 2), K);   # fix max K to avoid out of bound error
    D = D(:,1:K)';
    I = I(:,1:K)';
  endif

endfunction


## Test output
%!shared x, y
%! x = [1, 1, 1; 2, 2, 2; 3, 3, 3];
%! y = [0, 0, 0; 1, 2, 3; 0, 2, 4; 4, 7, 1];
%!test
%! d = sqrt([3, 5, 11, 45; 12, 2, 8, 30; 27, 5, 11, 21]);
%! assert (pdist2 (x, y), d);
%!test
%! d = [5.1962, 2.2361, 3.3166, 6.7082; ...
%!      3.4641, 2.2361, 3.3166, 5.4772];
%! i = [3, 1, 1, 1; 2, 3, 3, 2];
%! [D, I] = pdist2 (x, y, "euclidean", "largest", 2);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! d = [1.7321, 1.4142, 2.8284, 4.5826; ...
%!      3.4641, 2.2361, 3.3166, 5.4772];
%! i = [1, 2, 2, 3;2, 1, 1, 2];
%! [D, I] = pdist2 (x, y, "euclidean", "smallest", 2);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 6.1644, 5.3852; 1.4142, 6.9282, 8.7750; ...
%!      3.7417, 7.0711, 9.9499; 6.1644, 10.4881, 10.3441];
%! i = [2, 4, 4; 3, 2, 2; 1, 3, 3; 4, 1, 1];
%! [D, I] = pdist2 (y, yy, "euclidean", "smallest", 4);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 38, 29; 2, 48, 77; 14, 50, 99; 38, 110, 107];
%! i = [2, 4, 4; 3, 2, 2; 1, 3, 3; 4, 1, 1];
%! [D, I] = pdist2 (y, yy, "squaredeuclidean", "smallest", 4);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! d = [2.1213, 4.2426, 6.3640; 1.2247, 2.4495, 4.4159; ...
%!      3.2404, 4.8990, 6.8191; 2.7386, 4.2426, 6.1237];
%! assert (pdist2 (y, x, "mahalanobis"), d, 1e-4);
%!test
%! d = [1.2247, 2.4495, 4.4159; 2.1213, 4.2426, 6.1237];
%! i = [2, 2, 2; 1, 4, 4];
%! [D, I] = pdist2 (y, x, "mahalanobis", "smallest", 2);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! d = [3.2404, 4.8990, 6.8191; 2.7386, 4.2426, 6.3640];
%! i = [3, 3, 3; 4, 1, 1];
%! [D, I] = pdist2 (y, x, "mahalanobis", "largest", 2);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 8.4853, 18.0416; 2.4495, 10.0995, 19.4808; ...
%!      2.4495, 10.6771, 19.7104; 2.4495, 10.6771, 20.4573];
%! i = [2, 2, 2; 1, 4, 4; 4, 1, 1; 3, 3, 3];
%! [D, I] = pdist2 (y, yy, "mahalanobis", "smallest", 4);
%! assert ({D, I}, {d, i}, 1e-4);
%!test
%! d = [3, 3, 5, 9; 6, 2, 4, 8; 9, 3, 5, 7];
%! assert (pdist2 (x, y, "cityblock"), d);
%!test
%! d = [1, 2, 3, 6; 2, 1, 2, 5; 3, 2, 3, 4];
%! assert (pdist2 (x, y, "chebychev"), d);
%!test
%! d = repmat ([NaN, 0.0742, 0.2254, 0.1472], [3, 1]);
%! assert (pdist2 (x, y, "cosine"), d, 1e-4);
%!test
%! yy = [1 2 3;5 6 7;9 5 1];
%! d = [0, 0, 0.5; 0, 0, 2; 1.5, 1.5, 2; NaN, NaN, NaN];
%! i = [2, 2, 4; 3, 3, 2; 4, 4, 3; 1, 1, 1];
%! [D, I] = pdist2 (y, yy, "correlation", "smallest", 4);
%! assert ({D, I}, {d, i}, eps);
%! [D, I] = pdist2 (y, yy, "spearman", "smallest", 4);
%! assert ({D, I}, {d, i}, eps);
%!test
%! d = [1, 2/3, 1, 1; 1, 2/3, 1, 1; 1, 2/3, 2/3, 2/3];
%! i = [1, 1, 1, 2; 2, 2, 3, 3; 3, 3, 2, 1];
%! [D, I] = pdist2 (x, y, "hamming", "largest", 4);
%! assert ({D, I}, {d, i}, eps);
%! [D, I] = pdist2 (x, y, "jaccard", "largest", 4);
%! assert ({D, I}, {d, i}, eps);

## Test input validation
%!error<pdist2: too few input arguments.> pdist2 (1)
%!error<pdist2: X and Y must have equal number of columns.> ...
%! pdist2 (ones (4, 5), ones (4))
%!error<pdist2: X and Y must be 2 dimensional matrices.> ...
%! pdist2 (ones (4, 2, 3), ones (3, 2))
%!error<pdist2: missing value in optional name/value paired arguments.> ...
%! pdist2(ones (3), ones (3), "euclidean", "Largest")
%!error<pdist2: missing value in optional name/value paired arguments.> ...
%! pdist2(ones (3), ones (3), "minkowski", 3, "Largest")
%!error<pdist2: invalid NAME in optional pairs of arguments.> ...
%! pdist2(ones (3), ones (3), "minkowski", 3, "large", 4)
%!error<pdist2: you can only use either Smallest or Largest.> ...
%! pdist2(ones (3), ones (3), "minkowski", 3, "Largest", 4, "smallest", 5)
%!error<pdist2: Smallest or Largest must be specified to compute second output.> ...
%! [d, i] = pdist2(ones (3), ones (3), "minkowski", 3)
%!error<pdist2: DistParameter for standardized euclidean must be a vector of> ...
%! pdist2(ones (3), ones (3), "seuclidean", 3)
%!error<pdist2: DistParameter for standardized euclidean must be a nonnegative> ...
%! pdist2(ones (3), ones (3), "seuclidean", [1, -1, 3])
%!error<pdist2: DistParameter for mahalanobis distance must be a covariance> ...
%! pdist2(ones (3), eye (3), "mahalanobis", eye(2))
%!error<pdist2: covariance matrix for mahalanobis distance must be symmetric> ...
%! pdist2(ones (3), eye (3), "mahalanobis", ones(3))
%!error<pdist2: DistParameter for minkowski distance must be a positive scalar.> ...
%! pdist2(ones (3), eye (3), "minkowski", 0)
%!error<pdist2: DistParameter for minkowski distance must be a positive scalar.> ...
%! pdist2(ones (3), eye (3), "minkowski", -5)
%!error<pdist2: DistParameter for minkowski distance must be a positive scalar.> ...
%! pdist2(ones (3), eye (3), "minkowski", [1, 2])
