## Copyright (C) 2014 JD Walsh <walsh@math.gatech.edu>
## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn  {statistics} {@var{Y} =} cmdscale (@var{D})
## @deftypefnx {statistics} {[@var{Y}, @var{e}] =} cmdscale (@var{D})
## @deftypefnx {statistics} {@var{Y} =} cmdscale (@var{D}, @var{p})
## @deftypefnx {statistics} {[@var{Y}, @var{e}] =} cmdscale (@var{D}, @var{p})
##
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
## The optional argument @var{p} is a positive integer between 1 and @var{n}
## that specifies the maximum dimensionality of the desired embedding @var{Y}.
## If specified, @var{Y} will have at most @var{p} columns, and the returned
## eigenvalues @var{e} will be a vector of exactly length @var{p}. Specifying
## @var{p} can be useful for reducing dimensions for visualization (e.g., @code{@var{p} = 2}).
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
## equal to one, with ones along the main diagonal. In this case the points
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
## (@code{@var{p} < @var{n}}). The columns correspond to the positive
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

function [Y, e] = cmdscale (D, p)

  ## Check for matrix input and valid number of arguments
  if (nargin < 1 || nargin > 2 || ! isnumeric (D) || ! ismatrix (D))
    error ("cmdscale: input must be a numeric vector or matrix.");
  endif

  ## If vector, convert to matrix; otherwise, check for square symmetric input
  if (isvector (D))
    D = squareform (D);
  elseif (! issquare (D) || norm (D - D', 1) > 0)
    error ("cmdscale: matrix input must be square symmetric.");
  endif

  n = size (D, 1);

  if (nargin > 1 && ! isempty (p))
    if (! isscalar (p) || ! isreal (p) || p != round (p) || p < 1 || p > n)
      error ("cmdscale: p must be an integer between 1 and %d.", n);
    endif
  endif

  ## Check for valid format and if similarity matrix, convert
  if (any (any (D < 0)))
    error ("cmdscale: entries must be nonnegative.");
  elseif (trace (D) != 0)
    if (! all (diag (D) == 1) || ! all (D(:) <= 1))
      error ("cmdscale: input must be a distance vector or matrix.");
    endif
    D = sqrt (ones (n,n) - D);
  endif

  ## Build centering matrix, perform double centering
  J = eye (n) - ones (n, n) / n;
  B = -1 / 2 * J * (D .^ 2) * J;

  B = (B + B') / 2;

  ## extract and sort eigenpairs
  [Q, e] = eig (B);
  etmp = diag (e);
  [etmp, ord] = sort (etmp, 'descend');
  Q = Q(:, ord);

  ## determine total strictly positive eigenvalues
  tol = n * max (abs (etmp)) * eps;
  n_pos = sum (etmp > tol);

  if (n_pos == n)
    n_pos = n - 1;
  endif

  if (nargin > 1 && ! isempty (p))
    ## if p is given, eigenvalue output length is exactly p
    e = etmp(1:p);
    k = min (p, n_pos);
  else
    e = etmp;
    k = n_pos;
  endif

  ## build output matrix Y by safely slicing the top k elements
  etmp = etmp(1:k);
  Q = Q(:, 1:k);
  Y = Q * diag (sqrt (etmp));

  [~, maxind] = max (abs (Y), [], 1);
  d = size (Y, 2);
  idx = maxind + (0 : n : (d - 1) * n);
  colsign = sign (Y(idx));
  colsign(colsign == 0) = 1; 
  Y = Y .* colsign;

endfunction

%!test
%! m = randi (100) + 1;
%! n = randi (100) + 1;
%! X = rand (m, n);
%! D = pdist (X);
%! assert (norm (pdist (cmdscale (D))), norm (D), sqrt (eps));
%! assert (norm (pdist (cmdscale (squareform (D)))), norm (D), sqrt (eps));
%!test
%! ## test output
%! X = [
%!   0.8147, 0.1576, 0.6557, 0.7060;
%!   0.9058, 0.9706, 0.0357, 0.0318;
%!   0.1270, 0.9572, 0.8491, 0.2769;
%!   0.9134, 0.4854, 0.9340, 0.0462;
%!   0.6324, 0.8003, 0.6787, 0.0971
%! ];
%! D = pdist (X);
%! p = 2;
%! [Y, e] = cmdscale (D, p);
%! expected_Y = [
%!    0.635444598081665, -0.209808014423477;
%!   -0.558655450609184, -0.457908993032377;
%!   -0.158680352453745,  0.622280326562354;
%!    0.222509398493731, -0.047804408953240;
%!   -0.140618193512467,  0.093241089846740
%! ];
%! expected_e = [0.810349112746116; 0.651912015993974];
%! assert (Y, expected_Y, 1e-14);
%! assert (e, expected_e, 1e-14);
%!test
%! ## basic dimentionality reduction
%! D = [0 2 3; 2 0 4; 3 4 0];
%! [Y, e] = cmdscale (D, 2);
%! assert (size (Y, 2), 2);
%! assert (length (e), 2);
%!test
%! ## oversized dimension
%! X = [0 0; 1 0; 0 1; 1 1];
%! D = pdist (X);
%! [Y, e] = cmdscale (D, 3);
%! assert (size (Y, 2), 2);
%! assert (length (e), 3);
%!test
%! ## non euclidean distance.
%! X = [1 2; 3 4; 5 6; 7 8; 9 10];
%! D = pdist (X, "cityblock");
%! [Y, e] = cmdscale (D, 2);
%! assert (size (Y, 2), 1);
%! assert (length (e), 2);
%!test
%! ## compatability with p
%! X = rand (10, 4);
%! D = pdist (X);
%! [Y, e] = cmdscale (D, 3);
%! assert (size (Y, 2), 3);
%! assert (length (e), 3);
%! assert (size (Y, 1), 10);
%!test
%! ## sign convention.
%! rng (0, "twister");
%! X = rand (10, 3);
%! D = pdist (X);
%! Y = cmdscale (D);
%! [~, maxind] = max (abs (Y), [], 1);
%! d = size (Y, 2);
%! n = size (Y, 1);
%! idx = maxind + (0 : n : (d - 1) * n);
%! assert (all (Y(idx) >= 0));
%!test
%! ## testing with p = n and without p
%! rng (1, "twister");
%! X = rand (10, 4);
%! D = pdist (X);
%! n_points = size (X, 1);
%! [Y1, e1] = cmdscale (D);
%! [Y2, e2] = cmdscale (D, n_points);
%! assert (size (Y1, 2), size (Y2, 2));
%!error <cmdscale: input must be a numeric vector or matrix.> cmdscale ({'not', 'a', 'matrix'})
%!error <matrix input must be square symmetric> cmdscale (rand (3, 4))
%!error <entries must be nonnegative> cmdscale (-ones (3))
%!error <p must be an integer> cmdscale (eye (3), 0)
%!error <p must be an integer> cmdscale (eye (3), 4)
%!error <p must be an integer> cmdscale (eye (3), 1.5)
%!error <p must be an integer> cmdscale (eye (3), [1, 2])
%!error <p must be an integer> cmdscale (eye (3), 2 + 1i)

