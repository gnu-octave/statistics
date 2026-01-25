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
## @deftypefn  {statistics} {@var{coeff} =} pcacov (@var{K})
## @deftypefnx {statistics} {[@var{coeff}, @var{latent}] =} pcacov (@var{K})
## @deftypefnx {statistics} {[@var{coeff}, @var{latent}, @var{explained}] =} pcacov (@var{K})
##
## Perform principal component analysis on covariance matrix
##
## @code{@var{coeff} = pcacov (@var{K})} performs principal component analysis
## on the square covariance matrix @var{K} and returns the principal component
## coefficients, also known as loadings.  The columns are in order of decreasing
## component variance.
##
## @code{[@var{coeff}, @var{latent}] = pcacov (@var{K})} also returns a vector
## with the principal component variances, i.e. the eigenvalues of @var{K}.
## @var{latent} has a length of @qcode{size (@var{coeff}, 1)}.
##
## @code{[@var{coeff}, @var{latent}, @var{explained}] = pcacov (@var{K})} also
## returns a vector with the percentage of the total variance explained by each
## principal component.  @var{explained} has the same size as @var{latent}.
## The entries in @var{explained} range from 0 (none of the variance is
## explained) to 100 (all of the variance is explained).
##
## @code{pcacov} does not standardize @var{K} to have unit variances.  In order
## to perform principal component analysis on standardized variables, use the
## correlation matrix @qcode{@var{R} = @var{K} ./ (@var{SD} * @var{SD}')}, where
## @qcode{@var{SD} = sqrt (diag (@var{K}))}, in place of @var{K}.  To perform
## principal component analysis directly on the data matrix, use @code{pca}.
##
## @subheading References
## @enumerate
## @item
## Jolliffe, I. T., Principal Component Analysis, 2nd Edition, Springer, 2002
## @end enumerate
##
## @seealso{barttest, factoran, pcares, pca}
## @end deftypefn

function [coeff, latent, explained] = pcacov (K)

  ## Check X being a square matrix
  if (ndims (K) != 2 || size (K, 1) != size (K, 2))
    error ("pcacov: K must be a square matrix.");
  endif

  ## Check X being a symmetric matrix
  if (! issymmetric (K))
    error ("pcacov: K must be a symmetric matrix.");
  endif

  [U, S, V] = svd (K);
  s_vals = diag (S);

  col_dot_prod = sum (real (conj (U) .* V), 1);
  tol = size (K, 1) * eps (max (s_vals));
  is_negative = col_dot_prod < -0.9;
  is_significant = s_vals' > tol;

  ## Check for positive semi-definiteness
  if (any (is_negative & is_significant))
     error ("pcacov: K must be a positive semi-definite matrix.");
  endif

  ## Force a sign convention on the coefficients so that
  ## the largest element in each column has a positive sign
  [row, col] = size (U);
  [~, m_ind] = max (abs (U), [], 1);
  csign = sign (U (m_ind + (0:row:(col - 1) * row)));
  coeff = bsxfun (@times, U, csign);

  ## Compute extra output arguments
  if (nargout > 1)
    latent = diag (S);
  endif
  if (nargout > 2)
    explained = 100 * latent ./ sum (latent);
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
%!      10    68     8    12
%!     ];
%! Kxx = cov (x);
%! [coeff, latent, explained] = pcacov (Kxx)

## Test output
%!test
%! load hald
%! Kxx = cov (ingredients);
%! [coeff,latent,explained] = pcacov(Kxx);
%! c_out = [-0.0678, -0.6460,  0.5673, 0.5062; ...
%!          -0.6785, -0.0200, -0.5440, 0.4933; ...
%!           0.0290,  0.7553,  0.4036, 0.5156; ...
%!           0.7309, -0.1085, -0.4684, 0.4844];
%! l_out = [517.7969; 67.4964; 12.4054; 0.2372];
%! e_out = [ 86.5974; 11.2882;  2.0747; 0.0397];
%! assert (coeff, c_out, 1e-4);
%! assert (latent, l_out, 1e-4);
%! assert (explained, e_out, 1e-4);

## Test input validation
%!error <pcacov: K must be a square matrix.> pcacov (ones (2, 3))
%!error <pcacov: K must be a square matrix.> pcacov (ones (3, 3, 3))
%!error <pcacov: K must be a symmetric matrix.> pcacov ([1, 2; 0, 1])
%!error <pcacov: K must be a positive semi-definite matrix.> pcacov ([10, 0; 0, -1])
