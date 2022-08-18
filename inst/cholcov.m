## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{T} = cholcov (@var{sigma})
## @deftypefnx {Function File} [@var{T}, @var{p} = cholcov (@var{sigma})
## @deftypefnx {Function File} [@dots{}] = cholcov (@var{sigma}, @var{flag})
##
## Cholesky-like decomposition for covariance matrix.
##
## @code{@var{T} = cholcov (@var{sigma})} computes matrix @var{T} such that
## @var{sigma} = @var{T}'  @var{T}.  @var{sigma} must be square, symmetric, and
## positive semi-definite.
##
## If @var{sigma} is positive definite, then @var{T} is the square, upper
## triangular Cholesky factor.  If @var{sigma} is not positive definite, @var{T}
## is computed with an eigenvalue decomposition of @var{sigma}, but in this case
## @var{T} is not necessarily triangular or square.  Any eigenvectors whose
## corresponding eigenvalue is close to zero (within a tolerance) are omitted.
## If any remaining eigenvalues are negative, @var{T} is empty.
##
## The tolerance is calculated as @code{10 * eps (max (abs (diag (sigma))))}.
##
## @code{[@var{T}, @var{p} = cholcov (@var{sigma})} returns in @var{p} the
## number of negative eigenvalues of @var{sigma}.  If @var{p} > 0, then @var{T}
## is empty, whereas if @var{p} = 0, @var{sigma}) is positive semi-definite.
##
## If @var{sigma} is not square and symmetric, P is NaN and T is empty.
##
## @code{[@var{T}, @var{p} = cholcov (@var{sigma}, 0)} returns @var{p} = 0 if
## @var{sigma} is positive definite, in which case @var{T} is the Cholesky
## factor.  If @var{sigma} is not positive definite, @var{p} is a positive
## integer and @var{T} is empty.
##
## @code{[@dots{}] = cholcov (@var{sigma}, 1)} is equivalent to
## @code{ [@dots{}] = cholcov (@var{sigma})}.
##
## @seealso{chov}
## @end deftypefn

function [T, p] = cholcov (sigma, flag)
  ## Check number of input arguments
  narginchk (1,2)
  ## Add default flag if not givens
  if (nargin < 2)
    flag = 1;
  endif
  ## Check if sigma is a sparse matrix
  is_sparse = issparse (sigma);
  ## Check if sigma is single or double class
  is_type = "double";
  if (isa (sigma, "single"))
    is_type = "single";
  endif
  ## Test for sigma being square and symmetric
  [col, row] = size (sigma);
  ## Add tolerance
  Tol = 10 * eps (max (abs (diag (sigma))));
  if ((row == col) && all (all (abs (sigma - sigma') < col * Tol)))
    ## Check if positive definite
    [T, p] = chol (sigma);
    if (p > 0)
      ## Check flag for factoring using eigenvalue decomposition
      if (flag)
        [V, LAMBDA] = eig (full ((sigma + sigma') / 2));
        [~, EIGMAX] = max (abs (V), [], 1);
        neg_idx = (V(EIGMAX + (0:row:(col-1)*row)) < 0);
        V(:,neg_idx) = -V(:,neg_idx);
        LAMBDA = diag(LAMBDA);
        Tol = eps (max (LAMBDA)) * length (LAMBDA);
        t = (abs (LAMBDA) > Tol);
        LAMBDA = LAMBDA(t);
        p = sum (LAMBDA < 0);
        ## Check for negative eigenvalues
        if (p == 0)
          T = diag (sqrt (LAMBDA)) * V(:,t)';
        else
          T = zeros (0, is_type);
        endif
      else
        T = zeros (0, is_type);
      endif
    endif
  else
    T = zeros (0, is_type);
    p = NaN (is_type);
  endif
  if (is_sparse)
    T = sparse(T);
  endif
endfunction

%!demo
%! C1 = [2, 1, 1, 2; 1, 2, 1, 2; 1, 1, 2, 2; 2, 2, 2, 3]
%! T = cholcov (C1)
%! C2 = T'*T

%!test
%! C1 = [2, 1, 1, 2; 1, 2, 1, 2; 1, 1, 2, 2; 2, 2, 2, 3];
%! T = cholcov (C1);
%! assert (C1, T'*T, 1e-15 * ones (size (C1)));

