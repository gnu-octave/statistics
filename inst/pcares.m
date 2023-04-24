## Copyright (C) 2013-2019 Fernando Damian Nieuwveldt <fdnieuwveldt@gmail.com>
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 3
## of the License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{residuals} =} pcares (@var{x}, @var{ndim})
## @deftypefnx {statistics} {[@var{residuals}, @var{reconstructed}] =} pcares (@var{x}, @var{ndim})
##
## Calculate residuals from principal component analysis.
##
## @code{@var{residuals} = pcares (@var{x}, @var{ndim})} returns the residuals
## obtained by retaining @var{ndim} principal components of the @math{NxD}
## matrix @var{x}. Rows of @var{x} correspond to observations, columns of
## @var{x} correspond to variables.  @var{ndim} is a scalar and must be less
## than or equal to @math{D}.  @var{residuals} is a matrix of the same size as
## @var{x}.  Use the data matrix, not the covariance matrix, with this function.
##
## @code{[@var{residuals}, @var{reconstructed}] = pcares (@var{x}, @var{ndim})}
## returns the reconstructed observations, i.e. the approximation to @var{x}
## obtained by retaining its first @var{ndim} principal components.
##
## @code{pcares} does not normalize the columns of @var{x}.  Use
## @qcode{pcares (zscore (@var{x}), @var{ndim})} in order to perform the
## principal components analysis based on standardized variables, i.e. based on
## correlations.  Use @code{pcacov} in order to perform principal components
## analysis directly on a covariance or correlation matrix without constructing
## residuals.
##
## @seealso{factoran, pcacov, pca}
## @end deftypefn

function [residuals, reconstructed] = pcares (x, ndim)

  ## Check input arguments
  if (nargin < 2)
    error ("pcares: too few input arguments.");
  endif
  if (ndim > size (x, 2))
    error ("pcares: NDIM must be less than or equal to the column of X.");
  endif

  ## Mean center data
  Xcentered = bsxfun (@minus, x, mean (x));

  ## Apply svd to get the principal component coefficients
  [U, S, V] = svd (Xcentered);

  ## Use only the first ndim PCA components
  v = V(:,1:ndim);

  ## Calculate the residuals
  residuals = Xcentered - Xcentered * (v * v');

  ## Compute extra output arguments
  if (nargout > 1)
    ## Reconstructed data using ndim PCA components
    reconstructed = x - residuals;
  endif

endfunction

%!demo
%! x = [ 7    26     6    60;
%!       1    29    15    52;
%!      11    56     8    20;
%!      11    31     8    47;
%!       7    52     6    33;
%!      11    55     9    22;
%!       3    71    17     6;
%!       1    31    22    44;
%!       2    54    18    22;
%!      21    47     4    26;
%!       1    40    23    34;
%!      11    66     9    12;
%!      10    68     8    12];
%!
%! ## As we increase the number of principal components, the norm
%! ## of the residuals matrix will decrease
%! r1 = pcares (x,1);
%! n1 = norm (r1)
%! r2 = pcares (x,2);
%! n2 = norm (r2)
%! r3 = pcares (x,3);
%! n3 = norm (r3)
%! r4 = pcares (x,4);
%! n4 = norm (r4)

## Test output
%!test
%! load hald
%! r1 = pcares (ingredients,1);
%! r2 = pcares (ingredients,2);
%! r3 = pcares (ingredients,3);
%! assert (r1(1,:), [2.0350,  2.8304, -6.8378, 3.0879], 1e-4);
%! assert (r2(1,:), [-2.4037, 2.6930, -1.6482, 2.3425], 1e-4);
%! assert (r3(1,:), [ 0.2008, 0.1957,  0.2045, 0.1921], 1e-4);

## Test input validation
%!error<pcares: too few input arguments.> pcares (ones (20, 3))
%!error<pcares: NDIM must be less than or equal to the column of X.> ...
%! pcares (ones (30, 2), 3)

