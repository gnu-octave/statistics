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
	  error ("knnsearch: number of rows in X and Y must match.");
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
    IS  = S(:) .^ (-2);
    dxy = X(ix(:),:) - Y(iy(:),:);
    D   = sqrt ((dxy .^ 2) .* IS);

  elseif (strcmpi (Distance, "mahalanobis"))
    if isempty(C)
      C = cov (X(! any (isnan (X), 2),:));
    endif
    dxy = X(ix(:),:) - Y(iy(:),:);
    D   = sqrt ( sum ((dxy  *inv (C)) .* dxy, 2));

  elseif (strcmpi (Distance, "minkowski"))
    D = sum (abs (X(ix(:),:) - Y(iy(:),:)) .^ P, 2) .^ (1 / P);

  elseif (strcmpi (Distance, "cityblock") || strcmpi (Distance, "manhattan'"))
    D = sum (abs (X(ix(:),:) - Y(iy(:),:)), 2);

  elseif (strcmpi (Distance, "cosine"))
    sx = sum (X .^ 2, 2) .^ (-1 / 2);
    sy = sum (Y .^ 2, 2) .^ (-1 / 2);
    D  = 1 - sum (X(ix(:),:) .* Y(iy(:),:), 2) .* sx(ix(:)) .* sy(iy(:));

  elseif (strcmp (Distance, "correlation"))
    D = 1 - corrcoef (Y', X');

  elseif (strcmpi (Distance, "spearman"))
    D = 1 - corrcoef (Y', X', "Rank");

  elseif (strcmpi (Distance, "hamming"))
    D = mean (abs (X(ix(:),:) != Y(iy(:),:)), 2);

  else
    error ("knnsearch: distance metric '%s' not supported yet", Distance);
  endif

  D = reshape (D, size (Y,1), size (X,1));
  if (K == 1)
    [dist,idx] = min (D, [], 2);
  else
    [dist,idx] = sort (D, 2);
    dist = dist(:,1:K);
    idx = idx(:,1:K);
  endif

endfunction
