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
## chosen from a von Mises distribution with mean direction parameter @var{mu}
## and concentration parameter @var{k} on the interval [-pi, pi].  The size of
## @var{r} is the common size of the input arguments. A scalar input functions
## as a constant matrix of the same size as the other inputs.
##
## When called with a single size argument, return a square matrix with the
## dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @seealso{vmcdf, vmpdf}
## @end deftypefn

function r = vmrnd (mu, k, varargin)

  if (nargin < 2)
    print_usage ();
  endif

  if (! isscalar (mu) || ! isscalar (k))
    [retval, mu, k] = common_size (mu, k);
    if (retval > 0)
      error ("vmrnd: MU and K must be of common size or scalars.");
    endif
  endif

  if (iscomplex (mu) || iscomplex (k))
    error ("unifrnd: MU and K must not be complex.");
  endif

  if (nargin == 2)
    sz = size (mu);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["vmrnd: dimension vector must be row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 3)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("vmrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  if (! isscalar (mu) && ! isequal (size (k), sz))
    error ("vmrnd: MU and K must be scalar or of size SZ.");
  endif

  ## Simulate!
  if (k < 1e-6)
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
endfunction


%!assert (size (vmrnd (1, 2)), [1, 1])
%!assert (size (vmrnd (1 * ones (2, 1), 2)), [2, 1])
%!assert (size (vmrnd (1 * ones (2, 2), 2)), [2, 2])
%!assert (size (vmrnd (1, 2 * ones (2, 1))), [2, 1])
%!assert (size (vmrnd (1, 2 * ones (2, 2))), [2, 2])
%!assert (size (vmrnd (1, 2, 3)), [3, 3])
%!assert (size (vmrnd (1, 2, [4, 1])), [4, 1])
%!assert (size (vmrnd (1, 2, 4, 1)), [4, 1])

## Test class of input preserved
%!assert (class (vmrnd (1, 2)), "double")
%!assert (class (vmrnd (single (1), 2)), "single")
%!assert (class (vmrnd (single ([1, 1]), 2)), "single")
%!assert (class (vmrnd (1, single (2))), "single")
%!assert (class (vmrnd (1, single ([2, 2]))), "single")

## Test input validation
%!error vmrnd ()
%!error vmrnd (1)
%!error vmrnd (ones (3), 2 * ones (2))
%!error vmrnd (ones (2), 2 * ones (3))
%!error vmrnd (i, 2)
%!error vmrnd (1, i)
%!error vmrnd (1, 2, -1)
%!error vmrnd (1, 2, ones (2))
%!error vmrnd (1, 2, [2 -1 2])
%!error vmrnd (1 * ones (2),2, 3)
%!error vmrnd (1 * ones (2),2, [3, 2])
%!error vmrnd (1 * ones (2),2, 3, 2)

