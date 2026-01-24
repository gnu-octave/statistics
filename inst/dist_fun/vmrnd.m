## Copyright (C) 2009 Soren Hauberg <soren@hauberg.org>
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
## @deftypefn  {statistics} {@var{r} =} vmrnd (@var{mu}, @var{k})
## @deftypefnx {statistics} {@var{r} =} vmrnd (@var{mu}, @var{k}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} vmrnd (@var{mu}, @var{k}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} vmrnd (@var{mu}, @var{k}, [@var{sz}])
##
## Random arrays from the von Mises distribution.
##
## @code{@var{r} = vmrnd (@var{mu}, @var{k})} returns an array of random angles
## chosen from a von Mises distribution with location parameter @var{mu} and
## concentration parameter @var{k} on the interval [-pi, pi].  The size of
## @var{r} is the common size of @var{mu} and @var{k}.  A scalar input functions
## as a constant matrix of the same size as the other inputs.  Both parameters
## must be finite real numbers and @var{k} > 0, otherwise NaN is returned.
##
## When called with a single size argument, @code{vmrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the von Mises distribution can be found at
## @url{https://en.wikipedia.org/wiki/Von_Mises_distribution}
##
## @seealso{vmcdf, vminv, vmpdf}
## @end deftypefn

function r = vmrnd (mu, k, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("vmrnd: function called with too few input arguments.");
  endif

  ## Check for common size of MU and Κ
  if (! isscalar (mu) || ! isscalar (k))
    [retval, mu, k] = common_size (mu, k);
    if (retval > 0)
      error ("vmrnd: MU and K must be of common size or scalars.");
    endif
  endif

  ## Check for MU and Κ being reals
  if (iscomplex (mu) || iscomplex (k))
    error ("vmrnd: MU and K must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (mu);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif ((isrow (varargin{1}) || isempty (varargin{1})) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat ("vmrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("vmrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (mu) && ! isequal (size (k), sz))
    error ("vmrnd: MU and K must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (k, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Handle zero size dimensions
  if (any (sz == 0))
    r = nan (sz, cls);
    return
  endif

  ## Simulate!
  if (all (k < 1e-6))
    ## k is small: sample uniformly on circle
    r = mu + (2 * pi * rand (sz) - pi);

  else
    a = 1 + sqrt (1 + 4 .* k .^ 2);
    b = (a - sqrt (2 .* a)) ./ (2 .* k);
    r_tmp = (1 + b .^ 2) ./ (2 .* b);

    N = prod (sz);
    if (isscalar (k))
      r_tmp = repmat (r_tmp, 1, N);
      k_tmp = repmat (k, 1, N);
      mu_rs = repmat (mu, 1, N);
    else
      r_tmp = reshape (r_tmp, 1, N);
      k_tmp = reshape (k, 1, N);
      mu_rs = reshape (mu, 1, N);
    endif
    notdone = true (N, 1);
    while (any (notdone))
      u (:, notdone) = (rand (3, N))(:,notdone);

      z (notdone) = (cos (pi .* u (1, :)))(notdone);
      f (notdone) = ((1 + r_tmp(notdone) .* z(notdone)) ./ (r_tmp(notdone) + z(notdone)));
      c (notdone) = (k_tmp(notdone) .* (r_tmp(notdone) - f(notdone)));

      notdone = (u (2, :) >= c .* (2 - c)) & (log (c) - log (u (2, :)) + 1 - c < 0);
      #N = sum (notdone);
    endwhile

    r = mu_rs + sign (u (3, :) - 0.5) .* acos (f);
    r = reshape (r, sz);
  endif

  ## Cast to appropriate class
  r = cast (r, cls);

endfunction

## Test output
%!assert (size (vmrnd (1, 1)), [1 1])
%!assert (size (vmrnd (1, ones (2,1))), [2, 1])
%!assert (size (vmrnd (1, ones (2,2))), [2, 2])
%!assert (size (vmrnd (ones (2,1), 1)), [2, 1])
%!assert (size (vmrnd (ones (2,2), 1)), [2, 2])
%!assert (size (vmrnd (1, 1, 3)), [3, 3])
%!assert (size (vmrnd (1, 1, [4, 1])), [4, 1])
%!assert (size (vmrnd (1, 1, 4, 1)), [4, 1])
%!assert (size (vmrnd (1, 1, 4, 1, 5)), [4, 1, 5])
%!assert (size (vmrnd (1, 1, 0, 1)), [0, 1])
%!assert (size (vmrnd (1, 1, 1, 0)), [1, 0])
%!assert (size (vmrnd (1, 1, 1, 2, 0, 5)), [1, 2, 0, 5])

## Test class of input preserved
%!assert (class (vmrnd (1, 1)), "double")
%!assert (class (vmrnd (1, single (1))), "single")
%!assert (class (vmrnd (1, single ([1, 1]))), "single")
%!assert (class (vmrnd (single (1), 1)), "single")
%!assert (class (vmrnd (single ([1, 1]), 1)), "single")

## Test input validation
%!error<vmrnd: function called with too few input arguments.> vmrnd ()
%!error<vmrnd: function called with too few input arguments.> vmrnd (1)
%!error<vmrnd: MU and K must be of common size or scalars.> ...
%! vmrnd (ones (3), ones (2))
%!error<vmrnd: MU and K must be of common size or scalars.> ...
%! vmrnd (ones (2), ones (3))
%!error<vmrnd: MU and K must not be complex.> vmrnd (i, 2, 3)
%!error<vmrnd: MU and K must not be complex.> vmrnd (1, i, 3)
%!error<vmrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! vmrnd (1, 2, -1)
%!error<vmrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! vmrnd (1, 2, 1.2)
%!error<vmrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! vmrnd (1, 2, ones (2))
%!error<vmrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! vmrnd (1, 2, [2 -1 2])
%!error<vmrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! vmrnd (1, 2, [2 0 2.5])
%!error<vmrnd: dimensions must be non-negative integers.> ...
%! vmrnd (1, 2, 2, -1, 5)
%!error<vmrnd: dimensions must be non-negative integers.> ...
%! vmrnd (1, 2, 2, 1.5, 5)
%!error<vmrnd: MU and K must be scalars or of size SZ.> ...
%! vmrnd (2, ones (2), 3)
%!error<vmrnd: MU and K must be scalars or of size SZ.> ...
%! vmrnd (2, ones (2), [3, 2])
%!error<vmrnd: MU and K must be scalars or of size SZ.> ...
%! vmrnd (2, ones (2), 3, 2)
