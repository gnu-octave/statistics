## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{r} =} lognrnd (@var{mu}, @var{sigma})
## @deftypefnx {statistics} {@var{r} =} lognrnd (@var{mu}, @var{sigma}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} lognrnd (@var{mu}, @var{sigma}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} lognrnd (@var{mu}, @var{sigma}, [@var{sz}])
##
## Random arrays from the lognormal distribution.
##
## @code{@var{r} = laplace_rnd (@var{mu}, @var{sigma})} returns an array of
## random numbers chosen from the lognormal distribution with parameters @var{mu}
## and @var{sigma}.  The size of @var{r} is the common size of @var{mu} and
## @var{sigma}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.  Both parameters must be reals and @var{sigma} > 0.  For
## @var{sigma} <= 0, NaN is returned.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{logncdf, logninv, lognpdf, lognstat}
## @end deftypefn

function r = lognrnd (mu, sigma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for common size of P, MU, and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("lognrnd: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (mu) || iscomplex (sigma))
    error ("lognrnd: MU and SIGMA must not be complex.");
  endif

  ## Check for SIZE vector or DIMENSION input arguments
  if (nargin == 2)
    sz = size (mu);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["lognrnd: dimension vector must be a row vector of"], ...
                     [" non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("lognrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("lognrnd: MU and SIGMA must be scalar or of size SZ.");
  endif

  ## Check for appropriate class
  if (isa (mu, "single") || isa (sigma, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from lognormal distribution
  if (isscalar (mu) && isscalar (sigma))
    if ((sigma > 0) && (sigma < Inf))
      r = exp (mu + sigma * randn (sz, cls));
    else
      r = NaN (sz, cls);
    endif
  else
    r = exp (mu + sigma .* randn (sz, cls));

    k = (sigma < 0) | (sigma == Inf);
    r(k) = NaN;
  endif

endfunction


%!assert (size (lognrnd (1,2)), [1, 1])
%!assert (size (lognrnd (ones (2,1), 2)), [2, 1])
%!assert (size (lognrnd (ones (2,2), 2)), [2, 2])
%!assert (size (lognrnd (1, 2*ones (2,1))), [2, 1])
%!assert (size (lognrnd (1, 2*ones (2,2))), [2, 2])
%!assert (size (lognrnd (1, 2, 3)), [3, 3])
%!assert (size (lognrnd (1, 2, [4 1])), [4, 1])
%!assert (size (lognrnd (1, 2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (lognrnd (1, 2)), "double")
%!assert (class (lognrnd (single (1), 2)), "single")
%!assert (class (lognrnd (single ([1 1]), 2)), "single")
%!assert (class (lognrnd (1, single (2))), "single")
%!assert (class (lognrnd (1, single ([2 2]))), "single")

## Test input validation
%!error lognrnd ()
%!error lognrnd (1)
%!error lognrnd (ones (3), ones (2))
%!error lognrnd (ones (2), ones (3))
%!error lognrnd (i, 2)
%!error lognrnd (2, i)
%!error lognrnd (1,2, -1)
%!error lognrnd (1,2, ones (2))
%!error lognrnd (1, 2, [2 -1 2])
%!error lognrnd (1,2, 1, ones (2))
%!error lognrnd (1,2, 1, -1)
%!error lognrnd (ones (2,2), 2, 3)
%!error lognrnd (ones (2,2), 2, [3, 2])
%!error lognrnd (ones (2,2), 2, 2, 3)
