## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software: you can redistribute it and/or
## modify it under the terms of the GNU General Public License as
## published by the Free Software Foundation, either version 3 of the
## License, or (at your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} @var{r} = mvnrnd (@var{mu}, @var{sigma})
## @deftypefnx {statistics} @var{r} = mvnrnd (@var{mu}, @var{sigma}, @var{n})
## @deftypefnx {statistics} @var{r} = mvnrnd (@var{mu}, @var{sigma}, @var{n}, @var{T})
## @deftypefnx {statistics} [@var{r}, @var{T}] = mvnrnd (@dots{})
##
## Random vectors from the multivariate normal distribution.
##
## @code{@var{r} = mvnrnd (@var{mu}, @var{sigma})} returns an N-by-D matrix
## @var{r} of random vectors chosen from the multivariate normal distribution
## with mean vector @var{mu} and covariance matrix @var{sigma}.  @var{mu} is an
## N-by-D matrix, and @code{mvnrnd} generates each N of @var{r} using the
## corresponding N of @var{mu}.  @var{sigma} is a D-by-D symmetric positive
## semi-definite matrix, or a D-by-D-by-N array.  If @var{sigma} is an array,
## @code{mvnrnd} generates each N of @var{r} using the corresponding page of
## @var{sigma}, i.e., @code{mvnrnd} computes @var{r(i,:)} using @var{mu(i,:)}
## and @var{sigma(:,:,i)}.  If the covariance matrix is diagonal, containing
## variances along the diagonal and zero covariances off the diagonal,
## @var{sigma} may also be specified as a 1-by-D matrix or a 1-by-D-by-N array,
## containing just the diagonal.  If @var{mu} is a 1-by-D vector, @code{mvnrnd}
## replicates it to match the trailing dimension of SIGMA.
##
## @code{@var{r} = mvnrnd (@var{mu}, @var{sigma}, @var{n})} returns a N-by-D
## matrix R of random vectors chosen from the multivariate normal distribution
## with 1-by-D mean vector @var{mu}, and D-by-D covariance matrix @var{sigma}.
##
## @code{@var{r} = mvnrnd (@var{mu}, @var{sigma}, @var{n}, @var{T})} supplies
## the Cholesky factor @var{T} of @var{sigma}, so that @var{sigma(:,:,J)} ==
## @var{T(:,:,J)}'*@var{T(:,:,J)} if @var{sigma} is a 3D array or @var{sigma} ==
## @var{T}'*@var{T} if @var{sigma} is a matrix.  No error checking is done on
## @var{T}.
##
## @code{[@var{r}, @var{T}] = mvnrnd (@dots{})} returns the Cholesky factor
## @var{T}, so it can be re-used to make later calls more efficient, although
## there are greater efficiency gains when SIGMA can be specified as a diagonal
## instead.
##
## @seealso{mvncdf, mvnpdf}
## @end deftypefn

function [r, T] = mvnrnd (mu, sigma, N, T)

  ## Check input arguments
  if (nargin < 2 || isempty (mu) || isempty (sigma))
    error ("mvnrnd: too few input arguments.");
  elseif (ndims (mu) > 2)
    error ("mvnrnd: wrong size of MU.");
  elseif (ndims (sigma) > 3)
    error ("mvnrnd: wrong size of SIGMA.");
  endif

  ## Get data type
  if (isa (mu, "single") || isa (sigma, "single"))
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Check whether sigma is passed as a diagonal or a matrix
  sd = size (sigma);
  if (sd(1) == 1 && sd(2) > 1)
    sd(1) = sd(2);
    is_diag = true;
  else
    is_diag = false;
  endif

  ## Get size of mean vector
  [rm, cm] = size (mu);
  ## Make sure MU is a row vector
  if (cm == 1 && rm == sd(1))
    mu = mu';
    [rm, cm] = size (mu);
  endif

  ## Check for valid N input argument
  if (nargin < 3 || isempty (N))
    N_empty = true;
  else
    N_empty = false;
    ## If MU is a row vector, rep it out to match N
    if (rm == 1)
      rm = N;
      mu = repmat (mu, rm, 1);
    elseif (rm != N)
      error ("mvnrnd: size mismatch of N and MU.");
    endif
  endif

  ## For single covariance matrix
  if (ndims (sigma) == 2)
    ## Check sigma for correct size
    if (sd(1) != sd(2))
      error ("mvnpdf: bad covariance matrix.");
    elseif (! sd(1) == cm)
      error ("mvnpdf: covariance matrix mismatch.");
    endif
    ## Check for Cholesky factor T
    if (nargin > 3)
      r = randn (rm, size (T, 1), is_class) * T + mu;
    elseif (is_diag)
      ## Check sigma for invalid values
      if (any (sigma <= 0))
        error ("mvnpdf: SIGMA diagonal contains negative or zero values.");
      endif
      t = sqrt (sigma);
      if (nargout > 1)
        T = diag (t);
      endif
      r = bsxfun (@times, randn (rm, cm, is_class), t) + mu;
    else
      ## Compute a Cholesky factorization
      [T, err] = cholcov (sigma);
      if (err != 0)
        error ("mvnrnd: covariance matrix is not positive definite.");
      endif
      r = randn (rm, size (T, 1), is_class) * T + mu;
    endif
  endif

  ## For multiple covariance matrices
  if (ndims (sigma) == 3)
    ## If MU is a row vector, rep it out to match sigma
    if (rm == 1 && N_empty)
      rm = sd(3);
      mu = repmat (mu, rm, 1);
    endif
    ## Check sigma for correct size
    if (sd(1) != sd(2))
      error ("mvnpdf: bad multiple covariance matrix.");
    elseif (sd(1) != cm)
      error ("mvnpdf: multiple covariance matrix mismatch.");
    elseif (sd(3) != rm)
      error ("mvnpdf: multiple covariance pages mismatch.");
    endif
    ## Check for Cholesky factor T
    if (nargin < 4)   # T not present
      if (nargout > 1)
        T = zeros (sd, is_class);
      endif
      if (is_diag)
        sigma = reshape(sigma,sd(2),sd(3))';
        ## Check sigma for invalid values
        if (any (sigma(:) <= 0))
          error ("mvnpdf: SIGMA diagonals contain negative or zero values.");
        endif
        R = sqrt(sigma);
        r = bsxfun (@times, randn (rm, cm, is_class), R) + mu;
        if (nargout > 1)
          for i = 1:rm
            T(:,:,i) = diag (R(i,:));
          endfor
        endif
      else
        r = zeros (rm, cm, is_class);
          for i = 1:rm
            [R, err] = cholcov (sigma(:,:,i));
            if (err != 0)
              error (strcat (["mvnrnd: multiple covariance matrix"], ...
                             [" is not positive definite."]));
            endif
            Rrows = size (R,1);
            r(i,:) = randn (1, Rrows, is_class) * R + mu(i,:);
            if (nargout > 1)
              T(1:Rrows,:,i) = R;
            endif
          endfor
      endif
    else              # T present
      r = zeros (rm, cm, is_class);
      for i = 1:rm
        r(i,:) = randn (1, cm, is_class) * T(:,:,i) + mu(i,:);
      endfor
    endif
  endif

endfunction

## Test input validation
%!error<mvnrnd: too few input arguments.> mvnrnd ()
%!error<mvnrnd: too few input arguments.> mvnrnd ([2, 3, 4])
%!error<mvnrnd: wrong size of MU.> mvnrnd (ones (2, 2, 2), ones (1, 2, 3, 4))
%!error<mvnrnd: wrong size of SIGMA.> mvnrnd (ones (1, 3), ones (1, 2, 3, 4))

## Output validation tests
%!assert (size (mvnrnd ([2, 3, 4], [2, 2, 2])), [1, 3])
%!assert (size (mvnrnd ([2, 3, 4], [2, 2, 2], 10)), [10, 3])
