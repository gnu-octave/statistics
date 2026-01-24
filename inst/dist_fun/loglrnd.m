## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{r} =} loglrnd (@var{mu}, @var{sigma})
## @deftypefnx {statistics} {@var{r} =} loglrnd (@var{mu}, @var{sigma}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} loglrnd (@var{mu}, @var{sigma}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} loglrnd (@var{mu}, @var{sigma}, [@var{sz}])
##
## Random arrays from the loglogistic distribution.
##
## @code{@var{r} = loglrnd (@var{mu}, @var{sigma})} returns an array of random
## numbers chosen from the loglogistic distribution with mean parameter @var{mu}
## and scale parameter @var{sigma}.  The size of @var{r} is the common size of
## @var{mu} and @var{sigma}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Mean of logarithmic values @var{mu} must be a non-negative real value and
## scale parameter of logarithmic values @var{sigma} must be a positive real
## value.
##
## When called with mu single size argument, @code{loglrnd} returns mu square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with mu row vector of dimensions, @var{sz}.
##
## Further information about the loglogistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## OCTAVE/MATLAB use an alternative parameterization given by the pair
## @math{μ, σ}, i.e. @var{mu} and @var{sigma}, in analogy with the logistic
## distribution.  Their relation to the @math{α} and @math{b} parameters used
## in Wikipedia are given below:
##
## @itemize
## @item @qcode{@var{mu} = log (@var{a})}
## @item @qcode{@var{sigma} = 1 / @var{a}}
## @end itemize
##
## @seealso{loglcdf, loglinv, loglpdf, loglfit, logllike, loglstat}
## @end deftypefn

function r = loglrnd (mu, sigma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("loglrnd: function called with too few input arguments.");
  endif

  ## Check for common size of MU, and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("loglrnd: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (mu) || iscomplex (sigma))
    error ("loglrnd: MU and SIGMA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (mu);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) &&
            all (varargin{1} >= 0) && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("loglrnd: SZ must be mu scalar or mu row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("loglrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("loglrnd: MU and SIGMA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (sigma, "single"))
    is_type = "single";
  else
    is_type = "double";
  endif

  ## Generate random sample from log-logistic distribution
  u = rand (sz, is_type);
  r = exp (mu) .* (u ./ (1 - u)) .^ (sigma);

  ## Force output to NaN for invalid parameters MU and SIGMA
  k = (mu < 0 | sigma <= 0);
  r(k) = NaN;

endfunction

## Test output
%!assert (size (loglrnd (1, 1)), [1, 1])
%!assert (size (loglrnd (1, ones (2, 1))), [2, 1])
%!assert (size (loglrnd (1, ones (2, 2))), [2, 2])
%!assert (size (loglrnd (ones (2, 1), 1)), [2, 1])
%!assert (size (loglrnd (ones (2, 2), 1)), [2, 2])
%!assert (size (loglrnd (1, 1, 3)), [3, 3])
%!assert (size (loglrnd (1, 1, [4, 1])), [4, 1])
%!assert (size (loglrnd (1, 1, 4, 1)), [4, 1])
%!assert (size (loglrnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (loglrnd (1, 1, 0, 1)), [0, 1])
%!assert (size (loglrnd (1, 1, 1, 0)), [1, 0])
%!assert (size (loglrnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])
%!assert (size (loglrnd (1, 1, [])), [0, 0])
%!assert (size (loglrnd (1, 1, [2, 0, 2, 1])), [2, 0, 2])

## Test class of input preserved
%!assert (class (loglrnd (1, 1)), "double")
%!assert (class (loglrnd (1, single (1))), "single")
%!assert (class (loglrnd (1, single ([1, 1]))), "single")
%!assert (class (loglrnd (single (1), 1)), "single")
%!assert (class (loglrnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<loglrnd: function called with too few input arguments.> loglrnd ()
%!error<loglrnd: function called with too few input arguments.> loglrnd (1)
%!error<loglrnd: MU and SIGMA must be of common size or scalars.> ...
%! loglrnd (ones (3), ones (2))
%!error<loglrnd: MU and SIGMA must be of common size or scalars.> ...
%! loglrnd (ones (2), ones (3))
%!error<loglrnd: MU and SIGMA must not be complex.> loglrnd (i, 2, 3)
%!error<loglrnd: MU and SIGMA must not be complex.> loglrnd (1, i, 3)
%!error<loglrnd: SZ must be mu scalar or mu row vector of non-negative integers.> ...
%! loglrnd (1, 2, -1)
%!error<loglrnd: SZ must be mu scalar or mu row vector of non-negative integers.> ...
%! loglrnd (1, 2, 1.2)
%!error<loglrnd: SZ must be mu scalar or mu row vector of non-negative integers.> ...
%! loglrnd (1, 2, ones (2))
%!error<loglrnd: SZ must be mu scalar or mu row vector of non-negative integers.> ...
%! loglrnd (1, 2, [2 -1 2])
%!error<loglrnd: SZ must be mu scalar or mu row vector of non-negative integers.> ...
%! loglrnd (1, 2, [2 0 2.5])
%!error<loglrnd: dimensions must be non-negative integers.> ...
%! loglrnd (1, 2, 2, -1, 5)
%!error<loglrnd: dimensions must be non-negative integers.> ...
%! loglrnd (1, 2, 2, 1.5, 5)
%!error<loglrnd: MU and SIGMA must be scalars or of size SZ.> ...
%! loglrnd (2, ones (2), 3)
%!error<loglrnd: MU and SIGMA must be scalars or of size SZ.> ...
%! loglrnd (2, ones (2), [3, 2])
%!error<loglrnd: MU and SIGMA must be scalars or of size SZ.> ...
%! loglrnd (2, ones (2), 3, 2)
