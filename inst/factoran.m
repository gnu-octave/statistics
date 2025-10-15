## Copyright (C) 2025 Swayam Shah <swayamshah66@gmail.com>
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
## @deftypefn  {statistics} {@var{loadings} =} factoran (@var{X}, @var{nfac})
## @deftypefnx {statistics} {[@var{loadings}, @var{specvar}] =} factoran (@var{X}, @var{nfac})
## @deftypefnx {statistics} {[@var{loadings}, @var{specvar}, @var{fscores}] =} factoran (@var{X}, @var{nfac})
##
## Perform principal axis factor analysis on data matrix.
##
## @code{@var{loadings} = factoran (@var{X}, @var{nfac})} performs principal
## axis factoring to extract @var{nfac} factors from the @math{N x P} data
## matrix @var{X}, where rows correspond to observations and columns to
## variables.  The output @var{loadings} is a @math{P x @var{nfac}} matrix
## whose columns contain the loadings on each factor, in decreasing order of
## importance.
##
## @code{[@var{loadings}, @var{specvar}] = factoran (@dots{})} also returns a
## @math{P x 1} vector @var{specvar} containing the specific variances (unique
## variances) for each variable.
##
## @code{[@var{loadings}, @var{specvar}, @var{fscores}] = factoran (@dots{})}
## also returns the @math{N x @var{nfac}} matrix @var{fscores} of estimated
## factor scores, computed using the regression method.
##
## The analysis is performed on the correlation matrix of the standardized
## @var{X}.  Initial communalities are set to 1.  Iterations continue until the
## maximum change in communality is less than 1e-4 or 50 iterations are reached.
## The sign of each loading vector is chosen so that the element with largest
## absolute value is positive.
##
## @subheading References
## @enumerate
## @item
## Harman, H. H., Modern Factor Analysis, 3rd Edition, University of Chicago
## Press, 1976.
## @end enumerate
##
## @seealso{barttest, pca, pcacov, pcares}
## @end deftypefn

function [loadings, specvar, fscores] = factoran (X, nfac)

  if (nargin < 2)
    print_usage ();
  endif

  [nobs, nvars] = size (X);

  if (! isnumeric (X) || ndims (X) != 2)
    error ("factoran: X must be a numeric matrix.");
  endif

  if (nfac < 1 || nfac >= nvars)
    error ("factoran: NFAC must be an integer between 1 and %d.", nvars - 1);
  endif

  ## Standardize the data
  Z = zscore (X);

  ## Correlation matrix
  R = corrcoef (Z);

  ## Initial communalities (squared diagonals, but since 1, set to 1)
  h2 = ones (nvars, 1);

  ## Iteration parameters
  maxit = 50;
  tol = 1e-4;

  converged = false;
  for it = 1 : maxit
    Rstar = R - diag (1 - h2);
    [V, D] = eig (Rstar);
    [~, idx] = sort (diag (D), "descend");
    V = V(:, idx);
    ev = diag (D(idx, idx));
    Ltmp = V(:, 1:nfac) * diag (sqrt (ev(1:nfac)));
    h2_new = sum (Ltmp .^ 2, 2);
    if (max (abs (h2_new - h2)) < tol)
      converged = true;
      break;
    endif
    h2 = h2_new;
  endfor

  if (! converged)
    warning ("factoran: did not converge in %d iterations.", maxit);
  endif

  loadings = Ltmp;

  ## Force sign convention: largest absolute value in each column positive
  [~, m_ind] = max (abs (loadings), [], 1);
  for j = 1 : nfac
    if (loadings (m_ind(j), j) < 0)
      loadings (:, j) = - loadings (:, j);
    endif
  endfor

  specvar = 1 - h2;

  ## Factor scores using regression method
  if (nargout > 2)
    if (any (specvar <= 0))
      warning ("factoran: non-positive specific variance, scores may be unreliable.");
    endif
    invPsi = diag (1 ./ specvar);
    tmp = loadings' * (invPsi * loadings);
    if (rcond (tmp) < eps)
      warning ("factoran: ill-conditioned matrix for factor scores.");
      fscores = NaN (nobs, nfac);
    else
      score_coeff = (tmp \ (loadings' * invPsi));  ## k x p
      B = score_coeff';  ## p x k
      fscores = Z * B;
    endif
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
%! [loadings, specvar, fscores] = factoran (x, 2);

## Test output
%!test
%! x = [1, 2; 2, 1; 3, 3];
%! [loadings, specvar, fscores] = factoran (x, 1);
%! l_out = [0.7071; 0.7071];
%! s_out = [0.5000; 0.5000];
%! f_out = [-0.7071; -0.7071; 1.4142];
%! assert (loadings, l_out, 1e-4);
%! assert (specvar, s_out, 1e-4);
%! assert (fscores, f_out, 1e-4);

## Test input validation
%!error factoran ()
%!error factoran (ones (5,3), 0)
%!error factoran (ones (5,3), 3)
%!error<factoran: X must be a numeric matrix.> factoran ({1,2}, 1)
%!error<factoran: X must be a numeric matrix.> factoran (ones (2,2,2), 1)
%!error<variables with zero variance> x=ones(3,2); x(:,2)=0; factoran (x,1)
