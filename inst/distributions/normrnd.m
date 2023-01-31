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
## @deftypefn  {statistics} @var{r} = normrnd (@var{mu}, @var{sigma})
## @deftypefnx {statistics} @var{r} = normrnd (@var{mu}, @var{sigma}, @var{rows})
## @deftypefnx {statistics} @var{r} = normrnd (@var{mu}, @var{sigma}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} @var{r} = normrnd (@var{mu}, @var{sigma}, [@var{sz}])
##
## Random arrays from the normal distribution.
##
## @code{@var{r} = normrnd (@var{mu}, @var{sigma})} returns an array of random
## numbers chosen from the normal distribution with parameters @var{mu} and
## @var{sigma}.  The size of @var{r} is the common size of @var{mu} and
## @var{sigma}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.  Both parameters must be finite real numbers and
## @var{sigma} > 0, otherwise NaN is returned.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{norminv, norminv, normpdf, normfit, normlike, normstat}
## @end deftypefn

function r = normrnd (mu, sigma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    print_usage ();
  endif

  ## Check for common MU and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [retval, mu, sigma] = common_size (mu, sigma);
    if (retval > 0)
      error ("normrnd: MU and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for MU and SIGMA being reals
  if (iscomplex (mu) || iscomplex (sigma))
    error ("normrnd: MU and SIGMA must not be complex.");
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
      error (strcat (["normrnd: dimension vector must be row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("normrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("normrnd: MU and SIGMA must be scalar or of size SZ.");
  endif

  ## Check for appropriate class
  if (isa (mu, "single") || isa (sigma, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Laplace distribution
  if (isscalar (mu) && isscalar (sigma))
    if (isfinite (mu) && (sigma >= 0) && (sigma < Inf))
      r = mu + sigma * randn (sz, cls);
    else
      r = NaN (sz, cls);
    endif
  else
    r = mu + sigma .* randn (sz, cls);
    k = ! isfinite (mu) | !(sigma >= 0) | !(sigma < Inf);
    r(k) = NaN;
  endif

endfunction


%!assert (size (normrnd (1,2)), [1, 1])
%!assert (size (normrnd (ones (2,1), 2)), [2, 1])
%!assert (size (normrnd (ones (2,2), 2)), [2, 2])
%!assert (size (normrnd (1, 2*ones (2,1))), [2, 1])
%!assert (size (normrnd (1, 2*ones (2,2))), [2, 2])
%!assert (size (normrnd (1, 2, 3)), [3, 3])
%!assert (size (normrnd (1, 2, [4 1])), [4, 1])
%!assert (size (normrnd (1, 2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (normrnd (1, 2)), "double")
%!assert (class (normrnd (single (1), 2)), "single")
%!assert (class (normrnd (single ([1 1]), 2)), "single")
%!assert (class (normrnd (1, single (2))), "single")
%!assert (class (normrnd (1, single ([2 2]))), "single")

## Test input validation
%!error normrnd ()
%!error normrnd (1)
%!error normrnd (ones (3), ones (2))
%!error normrnd (ones (2), ones (3))
%!error normrnd (i, 2)
%!error normrnd (2, i)
%!error normrnd (1,2, -1)
%!error normrnd (1,2, ones (2))
%!error normrnd (1, 2, [2 -1 2])
%!error normrnd (1,2, 1, ones (2))
%!error normrnd (1,2, 1, -1)
%!error normrnd (ones (2,2), 2, 3)
%!error normrnd (ones (2,2), 2, [3, 2])
%!error normrnd (ones (2,2), 2, 2, 3)
