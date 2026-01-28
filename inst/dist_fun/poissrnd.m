## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{r} =} poissrnd (@var{lambda})
## @deftypefnx {statistics} {@var{r} =} poissrnd (@var{lambda}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} poissrnd (@var{lambda}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} poissrnd (@var{lambda}, [@var{sz}])
##
## Random arrays from the Poisson distribution.
##
## @code{@var{r} = normrnd (@var{lambda})} returns an array of random numbers
## chosen from the Poisson distribution with rate parameter @var{lambda}.  The
## size of @var{r} is the common size of @var{lambda}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.  @var{lambda} must
## be a finite real number and greater or equal to 0, otherwise @qcode{NaN} is
## returned.
##
## When called with a single size argument, @code{poissrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the Poisson distribution can be found at
## @url{https://en.wikipedia.org/wiki/Poisson_distribution}
##
## @seealso{poisscdf, poissinv, poisspdf, poissfit, poisslike, poisstat}
## @end deftypefn

function r = poissrnd (lambda, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1)
    error ("poissrnd: function called with too few input arguments.");
  endif

  ## Check for LAMBDA being real
  if (iscomplex (lambda))
    error ("poissrnd: LAMBDA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 1)
    sz = size (lambda);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0)
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif (isempty (varargin{1}))
      r = [];
      return;
    else
      error (strcat ("poissrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 2)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("poissrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  ## Use 'size (ones (sz))' to ignore any trailing singleton dimensions in SZ
  if (! isscalar (lambda) && ! isequal (size (lambda), size (ones (sz))))
    error ("poissrnd: LAMBDA must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (lambda, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Poisson distribution
  if (isscalar (lambda))
    if (lambda >= 0 && lambda < Inf)
      r = randp (lambda, sz, cls);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);
    k = (lambda >= 0) & (lambda < Inf);
    r(k) = randp (lambda(k), cls);
  endif

endfunction

## Test output
%!assert (size (poissrnd (2)), [1, 1])
%!assert (size (poissrnd (ones (2, 1))), [2, 1])
%!assert (size (poissrnd (ones (2, 2))), [2, 2])
%!assert (size (poissrnd (1, 3)), [3, 3])
%!assert (size (poissrnd (1, [4, 1])), [4, 1])
%!assert (size (poissrnd (1, 4, 1)), [4, 1])
%!assert (size (poissrnd (1, 4, 1)), [4, 1])
%!assert (size (poissrnd (1, 4, 1, 5)), [4, 1, 5])
%!assert (size (poissrnd (1, 0, 1)), [0, 1])
%!assert (size (poissrnd (1, 1, 0)), [1, 0])
%!assert (size (poissrnd (1, 1, 2, 0, 5)), [1, 2, 0, 5])
%!assert (size (poissrnd (1, [])), [0, 0])
%!assert (size (poissrnd (1, [2, 0, 2, 1])), [2, 0, 2])
%!assert (poissrnd (0, 1, 1), 0)
%!assert (poissrnd ([0, 0, 0], [1, 3]), [0 0 0])

## Test class of input preserved
%!assert (class (poissrnd (2)), "double")
%!assert (class (poissrnd (single (2))), "single")
%!assert (class (poissrnd (single ([2 2]))), "single")

## Test input validation
%!error<poissrnd: function called with too few input arguments.> poissrnd ()
%!error<poissrnd: LAMBDA must not be complex.> poissrnd (i)
%!error<poissrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! poissrnd (1, -1)
%!error<poissrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! poissrnd (1, 1.2)
%!error<poissrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! poissrnd (1, ones (2))
%!error<poissrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! poissrnd (1, [2 -1 2])
%!error<poissrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! poissrnd (1, [2 0 2.5])
%!error<poissrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! poissrnd (ones (2), ones (2))
%!error<poissrnd: dimensions must be non-negative integers.> ...
%! poissrnd (1, 2, -1, 5)
%!error<poissrnd: dimensions must be non-negative integers.> ...
%! poissrnd (1, 2, 1.5, 5)
%!error<poissrnd: LAMBDA must be scalar or of size SZ.> poissrnd (ones (2,2), 3)
%!error<poissrnd: LAMBDA must be scalar or of size SZ.> poissrnd (ones (2,2), [3, 2])
%!error<poissrnd: LAMBDA must be scalar or of size SZ.> poissrnd (ones (2,2), 2, 3)
