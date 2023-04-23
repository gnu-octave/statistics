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
## @deftypefn  {statistics} {@var{r} =} loglrnd (@var{alpha}, @var{beta})
## @deftypefnx {statistics} {@var{r} =} loglrnd (@var{alpha}, @var{beta}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} loglrnd (@var{alpha}, @var{beta}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} loglrnd (@var{alpha}, @var{beta}, [@var{sz}])
##
## Random arrays from the log-logistic distribution.
##
## @code{@var{r} = loglrnd (@var{alpha}, @var{beta})} returns an array of random
## numbers chosen from the log-logistic distribution with scale parameter
## @var{alpha} and shape parameter @var{beta}.  The size of @var{r} is the
## common size of @var{alpha} and @var{beta}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.  Both parameters must
## be positive reals, otherwise @qcode{NaN} is returned.
##
## When called with a single size argument, it returns a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## Further information about the log-logistic distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-logistic_distribution}
##
## MATLAB compatibility: MATLAB uses an alternative parameterization given by
## the pair @math{μ, s}, i.e. @var{mu} and @var{scale}, in analogy with the
## logistic distribution.  Their relation to the @var{alpha} and @var{beta}
## parameters is given below:
##
## @itemize
## @item @qcode{@var{alpha} = exp (@var{mu})}
## @item @qcode{@var{beta} = 1 / @var{scale}}
## @end itemize
##
## @seealso{loglcdf, loglinv, loglpdf, loglfit, logllike, loglstat}
## @end deftypefn

function r = loglrnd (alpha, beta, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for common size of ALPHA, and BETA
  if (! isscalar (alpha) || ! isscalar (beta))
    [retval, alpha, beta] = common_size (alpha, beta);
    if (retval > 0)
      error ("loglrnd: ALPHA and BETA must be of common size or scalars.");
    endif
  endif

  ## Check for X, ALPHA, and BETA being reals
  if (iscomplex (alpha) || iscomplex (beta))
    error ("loglrnd: ALPHA and BETA must not be complex.");
  endif

  ## Check for SIZE vector or DIMENSION input arguments
  if (nargin == 2)
    sz = size (alpha);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["loglrnd: dimension vector must be a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("loglrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (alpha) && ! isequal (size (alpha), sz))
    error ("loglrnd: ALPHA and BETA must be scalars or of size SZ.");
  endif

  ## Check for appropriate class
  if (isa (alpha, "single") || isa (beta, "single"))
    is_type = "single";
  else
    is_type = "double";
  endif

  ## Generate random sample from Laplace distribution
  p = rand (sz, is_type);
  r = alpha .* (p ./ (1 - p)) .^ (1 ./ beta);

  ## Force output to NaN for invalid parameters ALPHA and BETA
  k = (alpha <= 0 | beta <= 0);
  r(k) = NaN;

endfunction

## Test results
%!assert (size (loglrnd (1, 1, 1)), [1, 1])
%!assert (size (loglrnd (1, 1, 2)), [2, 2])
%!assert (size (loglrnd (1, 1, [2, 1])), [2, 1])
%!assert (size (loglrnd (1, zeros (2, 2))), [2, 2])
%!assert (size (loglrnd (1, ones (2, 1))), [2, 1])
%!assert (size (loglrnd (1, ones (2, 2))), [2, 2])
%!assert (size (loglrnd (ones (2, 1), 1)), [2, 1])
%!assert (size (loglrnd (ones (2, 2), 1)), [2, 2])
%!assert (size (loglrnd (1, 1, 3)), [3, 3])
%!assert (size (loglrnd (1, 1, [4 1])), [4, 1])
%!assert (size (loglrnd (1, 1, 4, 1)), [4, 1])
%!test
%! r =  loglrnd (1, [1, 0, -1]);
%! assert (r([2:3]), [NaN, NaN])

## Test class of input preserved
%!assert (class (loglrnd (1, 0)), "double")
%!assert (class (loglrnd (1, single (0))), "single")
%!assert (class (loglrnd (1, single ([0 0]))), "single")
%!assert (class (loglrnd (1, single (1))), "single")
%!assert (class (loglrnd (1, single ([1 1]))), "single")
%!assert (class (loglrnd (single (1), 1)), "single")
%!assert (class (loglrnd (single ([1 1]), 1)), "single")

## Test input validation
%!error loglrnd ()
%!error loglrnd (1)
%!error<loglrnd: ALPHA and BETA must be of common size or scalars.> ...
%! loglrnd (ones (3), ones (2))
%!error<loglrnd: ALPHA and BETA must be of common size or scalars.> ...
%! loglrnd (ones (2), ones (3))
%!error<loglrnd: ALPHA and BETA must not be complex.> loglrnd (i, 2)
%!error<loglrnd: ALPHA and BETA must not be complex.> loglrnd (1, i)
%!error<loglrnd: dimension vector must be a row vector of non-negative> ...
%! loglrnd (0, 1, [3, -1])
%!error<loglrnd: dimension vector must be a row vector of non-negative> ...
%! loglrnd (0, 1, -1)
%!error<loglrnd: dimensions must be non-negative integers.> ...
%! loglrnd (0, 1, 3, -1)
%!error<loglrnd: ALPHA and BETA must be scalars or of size SZ.> ...
%! loglrnd (2, ones (2), 3)
%!error<loglrnd: ALPHA and BETA must be scalars or of size SZ.> ...
%! loglrnd (2, ones (2), [3, 2])
%!error<loglrnd: ALPHA and BETA must be scalars or of size SZ.> ...
%! loglrnd (2, ones (2), 3, 2)
