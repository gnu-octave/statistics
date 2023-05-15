## Copyright (C) 2016 Dag Lyberg
## Copyright (C) 1995-2015 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
##
## This file is part of the statistics package for GNU Octave.
##
## Octave is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or (at
## your option) any later version.
##
## Octave is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} nakarnd (@var{mu}, @var{omega})
## @deftypefnx {statistics} {@var{r} =} nakarnd (@var{mu}, @var{omega}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} nakarnd (@var{mu}, @var{omega}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} nakarnd (@var{mu}, @var{omega}, [@var{sz}])
##
## Random arrays from the Nakagami distribution.
##
## @code{@var{r} = nakarnd (@var{mu}, @var{omega})} returns an array of random
## numbers chosen from the Nakagami distribution with shape parameter @var{mu}
## and spread parameter @var{omega}.  The size of @var{r} is the common size of
## @var{mu} and @var{omega}.  A scalar input functions as a constant matrix of
## the same size as the other inputs.
##
## Both parameters must be positive reals and @qcode{@var{mu} >= 0.5}.  For
## @qcode{@var{mu} < 0.5} or @qcode{@var{omega} <= 0}, @qcode{NaN} is returned.
##
## When called with a single size argument, @code{nakarnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the Nakagami distribution can be found at
## @url{https://en.wikipedia.org/wiki/Nakagami_distribution}
##
## @seealso{nakacdf, nakainv, nakapdf}
## @end deftypefn

function r = nakarnd (mu, omega, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("nakarnd: function called with too few input arguments.");
  endif

  ## Check for common size of MU and OMEGA
  if (! isscalar (mu) || ! isscalar (omega))
    [retval, mu, omega] = common_size (mu, omega);
    if (retval > 0)
      error ("nakarnd: MU and OMEGA must be of common size or scalars.");
    endif
  endif

  ## Check for MU and OMEGA being reals
  if (iscomplex (mu) || iscomplex (omega))
    error ("nakarnd: MU and OMEGA must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (mu);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["nakarnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("nakarnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (mu), sz))
    error ("nakarnd: MU and OMEGA must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (omega, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from Nakagami distribution
  if (isscalar (mu) && isscalar (omega))
    if ((0.5 <= mu) && (mu < Inf) && (0 < omega) && (omega < Inf))
      m_gamma = mu;
      w_gamma = omega / mu;
      r = gamrnd (m_gamma, w_gamma, sz);
      r = sqrt (r);
    else
      r = NaN (sz, cls);
    endif
  else
    r = NaN (sz, cls);
    k = (0.5 <= mu) & (mu < Inf) & (0 < omega) & (omega < Inf);
    m_gamma = mu;
    w_gamma = omega ./ mu;
    r(k) = gamrnd (m_gamma(k), w_gamma(k));
    r(k) = sqrt (r(k));
  endif

endfunction

## Test output
%!assert (size (nakarnd (1, 1)), [1 1])
%!assert (size (nakarnd (1, ones (2,1))), [2, 1])
%!assert (size (nakarnd (1, ones (2,2))), [2, 2])
%!assert (size (nakarnd (ones (2,1), 1)), [2, 1])
%!assert (size (nakarnd (ones (2,2), 1)), [2, 2])
%!assert (size (nakarnd (1, 1, 3)), [3, 3])
%!assert (size (nakarnd (1, 1, [4, 1])), [4, 1])
%!assert (size (nakarnd (1, 1, 4, 1)), [4, 1])
%!assert (size (nakarnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (nakarnd (1, 1, 0, 1)), [0, 1])
%!assert (size (nakarnd (1, 1, 1, 0)), [1, 0])
%!assert (size (nakarnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (nakarnd (1, 1)), "double")
%!assert (class (nakarnd (1, single (1))), "single")
%!assert (class (nakarnd (1, single ([1, 1]))), "single")
%!assert (class (nakarnd (single (1), 1)), "single")
%!assert (class (nakarnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<nakarnd: function called with too few input arguments.> nakarnd ()
%!error<nakarnd: function called with too few input arguments.> nakarnd (1)
%!error<nakarnd: MU and OMEGA must be of common size or scalars.> ...
%! nakarnd (ones (3), ones (2))
%!error<nakarnd: MU and OMEGA must be of common size or scalars.> ...
%! nakarnd (ones (2), ones (3))
%!error<nakarnd: MU and OMEGA must not be complex.> nakarnd (i, 2, 3)
%!error<nakarnd: MU and OMEGA must not be complex.> nakarnd (1, i, 3)
%!error<nakarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nakarnd (1, 2, -1)
%!error<nakarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nakarnd (1, 2, 1.2)
%!error<nakarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nakarnd (1, 2, ones (2))
%!error<nakarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nakarnd (1, 2, [2 -1 2])
%!error<nakarnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nakarnd (1, 2, [2 0 2.5])
%!error<nakarnd: dimensions must be non-negative integers.> ...
%! nakarnd (1, 2, 2, -1, 5)
%!error<nakarnd: dimensions must be non-negative integers.> ...
%! nakarnd (1, 2, 2, 1.5, 5)
%!error<nakarnd: MU and OMEGA must be scalars or of size SZ.> ...
%! nakarnd (2, ones (2), 3)
%!error<nakarnd: MU and OMEGA must be scalars or of size SZ.> ...
%! nakarnd (2, ones (2), [3, 2])
%!error<nakarnd: MU and OMEGA must be scalars or of size SZ.> ...
%! nakarnd (2, ones (2), 3, 2)
