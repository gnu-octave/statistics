## Copyright (C) 2014 Piotr Dollar <pdollar@gmail.com>
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
## @deftypefn  {Function File} {} pdist2 (@var{x}, @var{y})
## @deftypefnx {Function File} {} pdist2 (@var{x}, @var{y}, @var{metric})
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

## Taken from Piotr's Computer Vision Matlab Toolbox Version 2.52, with
## author permission to distribute under GPLv3

function D = pdist2 (X, Y, metric = "euclidean")

  if (nargin < 2 || nargin > 3)
    print_usage ();
  elseif (columns (X) != columns (Y))
    error ("pdist2: X and Y must have equal number of columns");
  elseif (ndims (X) != 2 || ndims (Y) != 2)
    error ("pdist2: X and Y must be 2 dimensional matrices");
  endif

  switch (tolower (metric))
    case "sqeuclidean",  D = distEucSq (X, Y);
    case "euclidean",     D = sqrt (distEucSq (X, Y));
    case "l1",            D = distL1 (X, Y);
    case "cosine",        D = distCosine (X, Y);
    case "emd",           D = distEmd (X, Y);
    case "chisq",         D = distChiSq (X, Y);
    otherwise
      error ("pdist2: unknown distance METRIC %s", metric);
  endswitch
  D = max (0, D);

endfunction

## TODO we could check the value of p and n first, and choose one
## or the other loop accordingly.
## L1 COMPUTATION WITH LOOP OVER p, FAST FOR SMALL p.
## function D = distL1( X, Y )
## m = size(X,1); n = size(Y,1); p = size(X,2);
## mOnes = ones(1,m); nOnes = ones(1,n); D = zeros(m,n);
## for i=1:p
##   yi = Y(:,i);  yi = yi( :, mOnes );
##   xi = X(:,i);  xi = xi( :, nOnes );
##   D = D + abs( xi-yi' );
## end

function D = distL1 (X, Y)
  m = rows (X);
  n = rows (Y);
  mOnes = ones (1, m);
  D = zeros (m, n);
  for i = 1:n
    yi = Y(i,:);
    yi = yi(mOnes,:);
    D(:,i) = sum (abs (X-yi), 2);
  endfor
endfunction

function D = distCosine (X, Y)
  p = columns (X);
  X = X ./ repmat (sqrt (sumsq (X, 2)), [1 p]);
  Y = Y ./ repmat (sqrt (sumsq (Y, 2)), [1 p]);
  D = 1 - X*Y';
endfunction

function D = distEmd (X, Y)
  Xcdf = cumsum (X,2);
  Ycdf = cumsum (Y,2);
  m = rows (X);
  n = rows (Y);
  mOnes = ones (1, m);
  D = zeros (m, n);
  for i=1:n
    ycdf = Ycdf(i,:);
    ycdfRep = ycdf(mOnes,:);
    D(:,i) = sum (abs (Xcdf - ycdfRep), 2);
  endfor
endfunction

function D = distChiSq (X, Y)
  ## note: supposedly it's possible to implement this without a loop!
  m = rows (X);
  n = rows (Y);
  mOnes = ones (1, m);
  D = zeros (m, n);
  for i = 1:n
    yi = Y(i, :);
    yiRep = yi(mOnes, :);
    s = yiRep + X;
    d = yiRep - X;
    D(:,i) = sum (d.^2 ./ (s+eps), 2);
  endfor
  D = D/2;
endfunction

function dists = distEucSq (x, y)
  xx = sumsq (x, 2);
  yy = sumsq (y, 2)';
  dists = bsxfun (@plus, xx, yy) - 2 * x * (y');
endfunction

## euclidean distance as loop for testing purposes
%!function dist = euclidean_distance (x, y)
%!  [m, p] = size (X);
%!  [n, p] = size (Y);
%!  D = zeros (m, n);
%!  for i = 1:n
%!    d = X - repmat (Y(i,:), [m 1]);
%!    D(:,i) = sumsq (d, 2);
%!  endfor
%!endfunction

%!test
%! x = [1 1 1; 2 2 2; 3 3 3];
%! y = [0 0 0; 1 2 3; 0 2 4; 4 7 1];
%! d = sqrt([  3    5   11   45
%!            12    2    8   30
%!            27    5   11   21]);
%! assert (pdist2 (x, y), d)

