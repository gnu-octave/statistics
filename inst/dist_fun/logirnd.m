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
## @deftypefn  {statistics} {@var{r} =} logirnd (@var{mu}, @var{scale})
## @deftypefnx {statistics} {@var{r} =} logirnd (@var{mu}, @var{scale}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} logirnd (@var{mu}, @var{scale}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} logirnd (@var{mu}, @var{scale}, [@var{sz}])
##
## Random arrays from the logistic distribution.
##
## @code{@var{r} = logirnd (@var{mu}, @var{scale})} returns an array of
## random numbers chosen from the logistic distribution with parameters @var{mu}
## and @var{scale}.  The size of @var{r} is the common size of @var{mu} and
## @var{scale}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.  Both parameters must be reals and @var{scale} > 0.  For
## @var{scale} <= 0, NaN is returned.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{logistic_cdf, logistic_inv, logistic_pdf}
## @end deftypefn

function r = logirnd (mu, scale, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for common size of MU, and SCALE
  if (! isscalar (mu) || ! isscalar (scale))
    [retval, mu, scale] = common_size (mu, scale);
    if (retval > 0)
      error ("logirnd: MU and SCALE must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SCALE being reals
  if (iscomplex (mu) || iscomplex (scale))
    error ("logirnd: MU and SCALE must not be complex.");
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
      error (strcat (["logirnd: dimension vector must be row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("logirnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("logirnd: MU and SCALE must be scalar or of size SZ.");
  endif

  ## Check for appropriate class
  if (isa (mu, "single") || isa (scale, "single"))
    is_type = "single";
  else
    is_type = "double";
  endif

  ## Generate random sample from Laplace distribution

  r = - log (1 ./ rand (sz, is_type) - 1) .* scale + mu;

  ## Force output to NaN for invalid parameter SCALE <= 0
  k = (scale <= 0);
  r(k) = NaN;

endfunction

## Test results
%!assert (size (logirnd (1, 1, 1)), [1, 1])
%!assert (size (logirnd (1, 1, 2)), [2, 2])
%!assert (size (logirnd (1, 1, [2, 1])), [2, 1])
%!assert (size (logirnd (1, zeros (2, 2))), [2, 2])
%!assert (size (logirnd (1, ones (2, 1))), [2, 1])
%!assert (size (logirnd (1, ones (2, 2))), [2, 2])
%!assert (size (logirnd (ones (2, 1), 1)), [2, 1])
%!assert (size (logirnd (ones (2, 2), 1)), [2, 2])
%!assert (size (logirnd (1, 1, 3)), [3, 3])
%!assert (size (logirnd (1, 1, [4 1])), [4, 1])
%!assert (size (logirnd (1, 1, 4, 1)), [4, 1])
%!test
%! r =  logirnd (1, [1, 0, -1]);
%! assert (r([2:3]), [NaN, NaN])

## Test class of input preserved
%!assert (class (logirnd (1, 0)), "double")
%!assert (class (logirnd (1, single (0))), "single")
%!assert (class (logirnd (1, single ([0 0]))), "single")
%!assert (class (logirnd (1, single (1))), "single")
%!assert (class (logirnd (1, single ([1 1]))), "single")
%!assert (class (logirnd (single (1), 1)), "single")
%!assert (class (logirnd (single ([1 1]), 1)), "single")

## Test input validation
%!error logirnd ()
%!error logirnd (1)
%!error<logirnd: MU and SCALE must be of common size or scalars.> ...
%! logirnd (ones (3), ones (2))
%!error<logirnd: MU and SCALE must be of common size or scalars.> ...
%! logirnd (ones (2), ones (3))
%!error<logirnd: MU and SCALE must not be complex.> logirnd (i, 2)
%!error<logirnd: MU and SCALE must not be complex.> logirnd (1, i)
%!error<logirnd: dimension vector must be row vector of non-negative> ...
%! logirnd (0, 1, [3, -1])
%!error<logirnd: dimension vector must be row vector of non-negative> ...
%! logirnd (0, 1, -1)
%!error<logirnd: dimensions must be non-negative integers.> ...
%! logirnd (0, 1, 3, -1)
%!error<logirnd: MU and SCALE must be scalar or of size SZ.> ...
%! logirnd (2, ones (2), 3)
%!error<logirnd: MU and SCALE must be scalar or of size SZ.> ...
%! logirnd (2, ones (2), [3, 2])
%!error<logirnd: MU and SCALE must be scalar or of size SZ.> ...
%! logirnd (2, ones (2), 3, 2)
