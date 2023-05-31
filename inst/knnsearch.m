## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{idx} =} knnsearch (@var{x}, @var{y})
## @deftypefnx {statistics} {@var{idx} =} knnsearch (@var{x}, @var{y}, @var{name}, @var{value})
## @deftypefnx {statistics} {[@var{idx}, @var{D}] =} knnsearch (@dots{})
##
## Find k-nearest neighbors using input data
##
##
## @seealso{rangesearch}
## @end deftypefn

function [idx, dist] = knnsearch (X, Y, varargin)

  ## Check input data
  if (nargin < 2)
	  error ("knnsearch: too few input arguments.");
  endif

  if (size (X, 2) != size (Y, 2))
	  error ("knnsearch: number of columns in X and Y must match.");
  endif

  ## Add default values
  K = 1;                    # Number of nearest neighbors
  P = 2;                    # Exponent for Minkowski distance
  S = [];                   # Scale for the standardized Euclidean distance
  C = [];                   # Covariance matrix for Mahalanobis distance
  BS = 50;                  # Maximum number of points per leaf node for Kd-tree
  SI = true;                # Sort returned indices according to distance
  Distance = "euclidean";   # Distance metric to be used
  NSMethod = "exhaustive";  # Nearest neighbor search method

  ## Parse additional parameters in Name/Value pairs
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))
      case "k"
        K = varargin{2};
      case "p"
        P = varargin{2};
      case "scale"
        S = varargin{2};
      case "cov"
        C = varargin{2};
      case "bucketsize"
        BS = varargin{2};
      case "sortindices"
        SI = varargin{2};
      case "distance"
        Distance = varargin{2};
      case "nsmethod"
        NSMethod = varargin{2};
      otherwise
        error ("knnsearch: invalid NAME in optional pairs of arguments.");
    endswitch
    varargin(1:2) = [];
  endwhile

  [ix, iy] = meshgrid (1:size (X, 1), 1:size (Y, 1));
  if strcmpi (Distance, "euclidean")
    D = sqrt (sum ((X(ix(:),:) - Y(iy(:),:)) .^ 2, 2));

  elseif strcmpi (Distance, "seuclidean")
    if (isempty (S))
      S = std (X, [], 1);
    endif
    mu = mean (X, 1);
    sx = (X - mu) ./ S;
    sy = (Y - mu) ./ S;
    D = sqrt (sum ((sx(ix(:),:) - sy(iy(:),:)) .^ 2, 2));

  elseif (strcmpi (Distance, "mahalanobis"))
    if isempty(C)
      C = cov (X(! any (isnan (X), 2),:));
    endif
    dxy = X(ix(:),:) - Y(iy(:),:);
    D   = sqrt (sum ((dxy  * inv (C)) .* dxy, 2));

  elseif (strcmpi (Distance, "minkowski"))
    D = sum (abs (X(ix(:),:) - Y(iy(:),:)) .^ P, 2) .^ (1 / P);

  elseif (strcmpi (Distance, "cityblock") || strcmpi (Distance, "manhattan"))
    D = sum (abs (X(ix(:),:) - Y(iy(:),:)), 2);

  elseif (strcmpi (Distance, "chebychev"))
    D = max (abs (X(ix(:),:) - Y(iy(:),:)), [], 2);

  elseif (strcmpi (Distance, "cosine"))
    sx = sum (X .^ 2, 2) .^ (-1 / 2);
    sy = sum (Y .^ 2, 2) .^ (-1 / 2);
    D  = 1 - sum (X(ix(:),:) .* Y(iy(:),:), 2) .* sx(ix(:)) .* sy(iy(:));

  elseif (strcmp (Distance, "correlation"))
    mX = mean (X(ix(:),:), 2);
    mY = mean (Y(iy(:),:), 2);
    xy = sum ((X(ix(:),:) - mX) .* (Y(iy(:),:) - mY), 2);
    xx = sqrt (sum ((X(ix(:),:) - mX) .* (X(ix(:),:) - mX), 2));
    yy = sqrt (sum ((Y(iy(:),:) - mY) .* (Y(iy(:),:) - mY), 2));
    D = 1 - (xy ./ (xx .* yy));

  elseif (strcmpi (Distance, "spearman"))
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

  elseif (strcmpi (Distance, "hamming"))
    D = mean (abs (X(ix(:),:) != Y(iy(:),:)), 2);

  elseif (strcmpi (Distance, "jaccard"))
    xy0 = (X(ix(:),:) != 0 | Y(iy(:),:) != 0);
    D = sum ((X(ix(:),:) != Y(iy(:),:)) & xy0, 2) ./ sum (xy0, 2);
  endif

  D = reshape (D, size (Y, 1), size (X, 1));
  if (K == 1)
    [dist, idx] = min (D, [], 2);
  else
    [dist, idx] = sort (D, 2);
    dist = dist(:,1:K);
    idx = idx(:,1:K);
  endif

endfunction

## Test output
%!shared x, y
%! x = [1, 2, 3, 4; 2, 3, 4, 5; 3, 4, 5, 6];
%! y = [1, 2, 2, 3; 2, 3, 3, 4];
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "euclidean");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "euclidean", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
%! assert (D, [sqrt(2), 3.162277660168380; sqrt(2), sqrt(2)], 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "seuclidean");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "seuclidean", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
%! assert (D, [sqrt(2), 3.162277660168380; sqrt(2), sqrt(2)], 1e-14);
%!test
%! xx = [1, 2; 1, 3; 2, 4; 3, 6];
%! yy = [2, 4; 2, 6];
%! [idx, D] = knnsearch (xx, yy, "Distance", "mahalanobis");
%! assert (idx, [3; 2]);
%! assert (D, [0; 3.162277660168377], 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "minkowski");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * sqrt (2));
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "minkowski", "p", 3);
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 1.259921049894873, 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "cityblock");
%! assert (idx, [1; 1]);
%! assert (D, [2; 2]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "manhattan");
%! assert (idx, [1; 1]);
%! assert (D, [2; 2]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "chebychev");
%! assert (idx, [1; 1]);
%! assert (D, [1; 1]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "cosine");
%! assert (idx, [2; 3]);
%! assert (D, [0.005674536395645; 0.002911214328620], 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "correlation");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 0.051316701949486, 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "spearman");
%! assert (idx, [1; 1]);
%! assert (D, ones (2, 1) * 0.051316701949486, 1e-14);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "hamming");
%! assert (idx, [1; 1]);
%! assert (D, [0.5; 0.5]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "jaccard");
%! assert (idx, [1; 1]);
%! assert (D, [0.5; 0.5]);
%!test
%! [idx, D] = knnsearch (x, y, "Distance", "jaccard", "k", 2);
%! assert (idx, [1, 2; 1, 2]);
%! assert (D, [0.5, 1; 0.5, 0.5]);

## Test input validation
%!error<knnsearch: too few input arguments.> knnsearch (1)
%!error<knnsearch: number of columns in X and Y must match.> ...
%! knnsearch (ones (4, 5), ones (4))
%!error<knnsearch: invalid NAME in optional pairs of arguments.> ...
%! knnsearch (ones (4, 2), ones (3, 2), "Distance", "euclidean", "some", "some")
