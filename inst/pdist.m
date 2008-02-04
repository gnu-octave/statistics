## Copyright (C) 2006, 2008  Bill Denney  <bill@denney.ws>
##
## This software is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This software is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this software; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {@var{y} =} pdist (@var{x})
## @deftypefnx {Function File} {@var{y} =} pdist (@var{x},
## @var{distfun})
## @deftypefnx {Function File} {@var{y} =} pdist (@var{x},
## @var{distfun}, @var{distfunarg}, @dots{})
## Return the distance between any two rows in @var{x}.
##
## @var{x} is the matrix (n x m) to determine the distance between.  If
## no @var{distfun} is given, then the 'euclidean' distance is assumed.
## @var{distfun} may be any of these or a function handle to a user
## defined function that takes two arguments distfun (@var{u}, @var{V})
## where @var{u} is a the row (1 x m) that is having its distance taken
## relative to @var{V} (a p x m matrix).
##
## The output vector, @var{y}, is (n - 1) * (n / 2) long where the
## distances are in the order [(1, 2); (1, 3); @dots{}; (2, 3); @dots{};
## (n-1, n)].
##
## Any additional arguments after the @var{distfun} are passed as
## distfun (@var{u}, @var{V}, @var{distfunarg1}, @var{distfunarg2}
## @dots{}).
##
## Pre-defined distance functions are:
##
## @table @samp
## @item "euclidean"
## Euclidean distance (default)
##
## @item "seuclidean"
## Standardized Euclidean distance. Each coordinate in the sum of
## squares is inverse weighted by the sample variance of that
## coordinate.
##
## @item "mahalanobis"
## Mahalanobis distance
##
## @item "cityblock"
## City Block metric (aka manhattan distance)
##
## @item "minkowski"
## Minkowski metric (with a default parameter 2)
##
## @item "cosine"
## One minus the cosine of the included angle between points (treated as
## vectors)
##
## @item "correlation"
## One minus the sample correlation between points (treated as
## sequences of values).
##
## @item "spearman"
## One minus the sample Spearman's rank correlation between
## observations, treated as sequences of values
##
## @item "hamming"
## Hamming distance, the percentage of coordinates that differ
##
## @item "jaccard"
## One minus the Jaccard coefficient, the percentage of nonzero
## coordinates that differ
##
## @item "chebychev"
## Chebychev distance (maximum coordinate difference)
## @end table
## @seealso{cluster,squareform}
## @end deftypefn

## Author: Bill Denney <denney@...>

function y = pdist (x, distfun, varargin)

  if (nargin < 1)
    print_usage ();
  elseif (nargin > 1) && ...
        ! (ischar (distfun) || ...
           strcmp (class(distfun), "function_handle"))
    error ("pdist: the distance function must be either a string or a \
	function handle.");
  endif

  if (nargin < 2)
    distfun = "euclidean";
  endif

  if (isempty (x))
    error ("pdist: x cannot be empty");
  elseif (length (size (x)) > 2)
    error ("pdist: x must be 1 or 2 dimensional");
  endif

  sx1 = size (x, 1);
  y = [];
  ## compute the distance
  for i = 1:sx1
    tmpd = feval (distfun, x(i,:), x(i+1:sx1,:), varargin{:});
    y = [y;tmpd(:)];
  endfor

endfunction

## the different standardized distance functions

function d = euclidean(u, v)
  d = sqrt (sum ((repmat (u, size (v,1), 1) - v).^2, 2));
endfunction

function d = seuclidean(u, v)
  ## FIXME
  error("Not implemented")
endfunction

function d = mahalanobis(u, v, p)
  repu = repmat (u, size (v,1), 1);
  d = (repu - v)' * inv (cov (repu, v)) * (repu - v);
  d = d.^(0.5);
endfunction

function d = cityblock(u, v)
  d = sum (abs (repmat (u, size(v,1), 1) - v), 2);
endfunction

function d = minkowski
  if (nargin < 3)
    p = 2;
  endif

  d = (sum (abs (repmat (u, size(v,1), 1) - v).^p, 2)).^(1/p);
endfunction

function d = cosine(u, v)
  repu = repmat (u, size (v,1), 1);
  d = dot (repu, v, 2) ./ (dot(repu, repu).*dot(v, v));
endfunction

function d = correlation(u, v)
  repu = repmat (u, size (v,1), 1);
  d = cor(repu, v);
endfunction

function d = spearman(u, v)
  repu = repmat (u, size (v,1), 1);
  d = spearman (repu, v);
endfunction

function d = hamming(u, v)
  ## Hamming distance, the percentage of coordinates that differ
  sv2 = size(v, 2);
  for i = 1:sv2
    v(:,i) = (v(:,i) == u(i));
  endfor
  d = sum (v,2)./sv2;
endfunction

function d = jaccard(u, v)
  ## Jaccard distance, one minus the percentage of non-zero coordinates
  ## that differ
  sv2 = size(v, 2);
  for i = 1:sv2
    v(:,i) = (v(:,i) == u(i)) && (u(i) || v(:,i));
  endfor
  d = 1 - sum (v,2)./sv2;
endfunction

function d = chebychev(u, v)
  repu = repmat (u, size (v,1), 1);
  d = max (abs (repu - v), [], 2);
endfunction
