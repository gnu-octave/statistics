## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} @var{r} = laplace_rnd (@var{mu}, @var{beta})
## @deftypefnx {statistics} @var{r} = laplace_rnd (@var{mu}, @var{beta}, @var{rows})
## @deftypefnx {statistics} @var{r} = laplace_rnd (@var{mu}, @var{beta}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} @var{r} = laplace_rnd (@var{mu}, @var{beta}, [@var{sz}])
##
## Random arrays from the Laplace distribution.
##
## @code{@var{r} = laplace_rnd (@var{mu}, @var{beta})} returns an array of
## random numbers chosen from the Laplace distribution with parameters @var{mu}
## and @var{beta}.  The size of @var{r} is the common size of @var{mu} and
## @var{beta}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.  Both parameters must be reals and @var{beta} > 0.  For
## @var{beta} <= 0, NaN is returned.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{laplace_cdf, laplace_inv, laplace_rnd}
## @end deftypefn

function r = laplace_rnd (mu, beta, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for common size of MU, and BETA
  if (! isscalar (mu) || ! isscalar (beta))
    [retval, mu, beta] = common_size (mu, beta);
    if (retval > 0)
      error ("laplace_rnd: MU and BETA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (mu) || iscomplex (beta))
    error ("laplace_rnd: MU and BETA must not be complex.");
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
      error (strcat (["laplace_rnd: dimension vector must be row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("laplace_rnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("laplace_rnd: MU and BETA must be scalar or of size SZ.");
  endif

  ## Check for appropriate class
  if (isa (mu, "single") || isa (beta, "single"))
    is_type = "single";
  else
    is_type = "double";
  endif

  ## Generate random sample from Laplace distribution
  tmp = rand (sz, is_type);
  r = ((tmp < 1/2) .* log (2 * tmp) - ...
       (tmp > 1/2) .* log (2 * (1 - tmp))) .* beta + mu;

  ## Force output to NaN for invalid parameter BETA <= 0
  k = (beta <= 0);
  r(k) = NaN;

endfunction

## Test results
%!assert (size (laplace_rnd (1, 1, 1)), [1, 1])
%!assert (size (laplace_rnd (1, 1, 2)), [2, 2])
%!assert (size (laplace_rnd (1, 1, [2, 1])), [2, 1])
%!assert (size (laplace_rnd (1, zeros (2, 2))), [2, 2])
%!assert (size (laplace_rnd (1, ones (2, 1))), [2, 1])
%!assert (size (laplace_rnd (1, ones (2, 2))), [2, 2])
%!assert (size (laplace_rnd (ones (2, 1), 1)), [2, 1])
%!assert (size (laplace_rnd (ones (2, 2), 1)), [2, 2])
%!assert (size (laplace_rnd (1, 1, 3)), [3, 3])
%!assert (size (laplace_rnd (1, 1, [4 1])), [4, 1])
%!assert (size (laplace_rnd (1, 1, 4, 1)), [4, 1])
%!test
%! r =  laplace_rnd (1, [1, 0, -1]);
%! assert (r([2:3]), [NaN, NaN])

## Test class of input preserved
%!assert (class (laplace_rnd (1, 0)), "double")
%!assert (class (laplace_rnd (1, single (0))), "single")
%!assert (class (laplace_rnd (1, single ([0 0]))), "single")
%!assert (class (laplace_rnd (1, single (1))), "single")
%!assert (class (laplace_rnd (1, single ([1 1]))), "single")
%!assert (class (laplace_rnd (single (1), 1)), "single")
%!assert (class (laplace_rnd (single ([1 1]), 1)), "single")

## Test input validation
%!error laplace_rnd ()
%!error laplace_rnd (1)
%!error<laplace_rnd: MU and BETA must be of common size or scalars.> ...
%! laplace_rnd (ones (3), ones (2))
%!error<laplace_rnd: MU and BETA must be of common size or scalars.> ...
%! laplace_rnd (ones (2), ones (3))
%!error<laplace_rnd: MU and BETA must not be complex.> laplace_rnd (i, 2)
%!error<laplace_rnd: MU and BETA must not be complex.> laplace_rnd (1, i)
%!error<laplace_rnd: dimension vector must be row vector of non-negative> ...
%! laplace_rnd (0, 1, [3, -1])
%!error<laplace_rnd: dimension vector must be row vector of non-negative> ...
%! laplace_rnd (0, 1, -1)
%!error<laplace_rnd: dimensions must be non-negative integers.> ...
%! laplace_rnd (0, 1, 3, -1)
%!error<laplace_rnd: MU and BETA must be scalar or of size SZ.> ...
%! laplace_rnd (2, ones (2), 3)
%!error<laplace_rnd: MU and BETA must be scalar or of size SZ.> ...
%! laplace_rnd (2, ones (2), [3, 2])
%!error<laplace_rnd: MU and BETA must be scalar or of size SZ.> ...
%! laplace_rnd (2, ones (2), 3, 2)
