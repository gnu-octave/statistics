## Copyright (C) 2014 JD Walsh <walsh@math.gatech.edu>
##
## This program is free software; you can redistribute it and/or modify 
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} @var{Y} = cmdscale (@var{D})
## @deftypefnx{Function File} [@var{Y}, @var{e} ] = cmdscale (@var{D})
## Classical multidimensional scaling of a matrix.
##
## Takes an @var{n} by @var{n} distance (or difference, similarity, or
## dissimilarity) matrix @var{D}. Returns @var{Y}, a matrix of @var{n} points
## with coordinates in @var{p} dimensional space which approximate those
## distances (or differences, similarities, or dissimilarities). Also returns
## the eigenvalues @var{e} of
## @code{@var{B} = -1/2 * @var{J} * (@var{D}.^2) * @var{J}}, where
## @code{J = eye(@var{n}) - ones(@var{n},@var{n})/@var{n}}. @var{p}, the number
## of columns of @var{Y}, is equal to the number of positive real eigenvalues of
## @var{B}.
##
## @var{D} can be a full or sparse matrix or a vector of length
## @code{@var{n}*(@var{n}-1)/2} containing the upper triangular elements (like
## the output of the @code{pdist} function). It must be symmetric with
## non-negative entries whose values are further restricted by the type of
## matrix being represented:
##
## * If @var{D} is either a distance, dissimilarity, or difference matrix, then
## it must have zero entries along the main diagonal. In this case the points
## @var{Y} equal or approximate the distances given by @var{D}.
##
## * If @var{D} is a similarity matrix, the elements must all be less than or
## equal to one, with ones along the the main diagonal. In this case the points
## @var{Y} equal or approximate the distances given by
## @code{@var{D} = sqrt(ones(@var{n},@var{n})-@var{D})}.
##
## @var{D} is a Euclidean matrix if and only if @var{B} is positive
## semi-definite. When this is the case, then @var{Y} is an exact representation
## of the distances given in @var{D}. If @var{D} is non-Euclidean, @var{Y} only
## approximates the distance given in @var{D}. The approximation used by
## @code{cmdscale} minimizes the statistical loss function known as
## @var{strain}.
##
## The returned @var{Y} is an @var{n} by @var{p} matrix showing possible
## coordinates of the points in @var{p} dimensional space
## (@code{@var{p} < @var{n}}). The columns are correspond to the positive
## eigenvalues of @var{B} in descending order. A translation, rotation, or
## reflection of the coordinates given by @var{Y} will satisfy the same distance
## matrix up to the limits of machine precision.
##
## For any @code{@var{k} <= @var{p}}, if the largest @var{k} positive
## eigenvalues of @var{B} are significantly greater in absolute magnitude than
## its other eigenvalues, the first @var{k} columns of @var{Y} provide a
## @var{k}-dimensional reduction of @var{Y} which approximates the distances
## given by @var{D}. The optional return @var{e} can be used to consider various
## values of @var{k}, or to evaluate the accuracy of specific dimension
## reductions (e.g., @code{@var{k} = 2}).
##
## Reference: Ingwer Borg and Patrick J.F. Groenen (2005), Modern
## Multidimensional Scaling, Second Edition, Springer, ISBN: 978-0-387-25150-9
## (Print) 978-0-387-28981-6 (Online)
##
## @seealso{pdist}
## @end deftypefn

## Author: JD Walsh <walsh@math.gatech.edu>
## Created: 2014-10-31
## Description: Classical multidimensional scaling
## Keywords: multidimensional-scaling mds distance clustering

## TO DO: include missing functions `mdscale' and `procrustes' in @seealso

function [Y, e] = cmdscale (D)

  % Check for matrix input
  if ((nargin ~= 1) || ...
      (~any(strcmp ({'matrix' 'scalar' 'range'}, typeinfo(D)))))
    usage ('cmdscale: input must be vector or matrix; see help');
  endif

  % If vector, convert to matrix; otherwise, check for square symmetric input
  if (isvector (D))
    D = squareform (D);
  elseif ((~issquare (D)) || (norm (D - D', 1) > 0))
    usage ('cmdscale: matrix input must be square symmetric; see help');
  endif

  n = size (D,1);
  % Check for valid format (see help above); If similarity matrix, convert
  if (any (any (D < 0)))
    usage ('cmdscale: entries must be nonnegative; see help');
  elseif (trace (D) ~= 0)
      if ((~all (diag (D) == 1)) || (~all (D <= 1)))
        usage ('cmdscale: input must be distance vector or matrix; see help');
      endif
      D = sqrt (ones (n,n) - D);
  endif

  % Build centering matrix, perform double centering, extract eigenpairs
  J = eye (n) - ones (n,n) / n;
  B = -1 / 2 * J * (D .^ 2) * J;
  [Q, e] = eig (B);
  e = diag (e);
  etmp = e;
  e = sort(e, 'descend');

  % Remove complex eigenpairs (only possible due to machine approximation)
  if (iscomplex (etmp))
    for i = 1 : size (etmp,1)
      cmp(i) = (isreal (etmp(i)));
    endfor
    etmp = etmp(cmp);
    Q = Q(:,cmp);
  endif

  % Order eigenpairs
  [etmp, ord] = sort (etmp, 'descend');
  Q = Q(:,ord);

  % Remove negative eigenpairs
  cmp = (etmp > 0);
  etmp = etmp(cmp);
  Q = Q(:,cmp);

  % Test for n-dimensional results
  if (size(etmp,1) == n)
    etmp = etmp(1:n-1);
    Q = Q(:, 1:n-1);
  endif

  % Build output matrix Y
  Y = Q * diag (sqrt (etmp));

endfunction

%!shared m, n, X, D
%! m = randi(100) + 1; n = randi(100) + 1; X = rand(m, n); D = pdist(X);
%!assert(norm(pdist(cmdscale(D))), norm(D), sqrt(eps))
%!assert(norm(pdist(cmdscale(squareform(D)))), norm(D), sqrt(eps))

