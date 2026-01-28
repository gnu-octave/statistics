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
## @deftypefn  {statistics} {@var{rnd} =} nbinrnd (@var{r}, @var{ps})
## @deftypefnx {statistics} {@var{rnd} =} nbinrnd (@var{r}, @var{ps}, @var{rows})
## @deftypefnx {statistics} {@var{rnd} =} nbinrnd (@var{r}, @var{ps}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{rnd} =} nbinrnd (@var{r}, @var{ps}, [@var{sz}])
##
## Random arrays from the negative binomial distribution.
##
## @code{@var{rnd} = nbinrnd (@var{r}, @var{ps})} returns an array of random
## numbers chosen from the negative binomial distribution with parameters
## @var{r} and @var{ps}, where @var{r} is the number of successes until the
## experiment is stopped and @var{ps} is the probability of success in each
## experiment, given the number of failures in @var{x}.  The size of @var{rnd}
## is the common size of @var{r} and @var{ps}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## When @var{r} is an integer, the negative binomial distribution is also known
## as the Pascal distribution and it models the number of failures in @var{x}
## before a specified number of successes is reached in a series of independent,
## identical trials.  Its parameters are the probability of success in a single
## trial, @var{ps}, and the number of successes, @var{r}.  A special case of the
## negative binomial distribution, when @qcode{@var{r} = 1}, is the geometric
## distribution, which models the number of failures before the first success.
##
## @var{r} can also have non-integer positive values, in which form the negative
## binomial distribution, also known as the Polya distribution, has no
## interpretation in terms of repeated trials, but, like the Poisson
## distribution, it is useful in modeling count data.  The negative binomial
## distribution is more general than the Poisson distribution because it has a
## variance that is greater than its mean, making it suitable for count data
## that do not meet the assumptions of the Poisson distribution.  In the limit,
## as @var{r} increases to infinity, the negative binomial distribution
## approaches the Poisson distribution.
##
## Further information about the negative binomial distribution can be found at
## @url{https://en.wikipedia.org/wiki/Negative_binomial_distribution}
##
## @seealso{nbininv, nbininv, nbinpdf, nbinfit, nbinlike, nbinstat}
## @end deftypefn

function rnd = nbinrnd (r, ps, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("nbinrnd: function called with too few input arguments.");
  endif

  ## Check for common size R and PS
  if (! isscalar (r) || ! isscalar (ps))
    [retval, r, ps] = common_size (r, ps);
    if (retval > 0)
      error ("nbinrnd: R and PS must be of common size or scalars.");
    endif
  endif

  ## Check for R and PS being reals
  if (iscomplex (r) || iscomplex (ps))
    error ("nbinrnd: R and PS must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (r);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} >= 0
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0)
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif (isempty (varargin{1}))
      rnd = [];
      return;
    else
      error (strcat ("nbinrnd: SZ must be a scalar or a row vector", ...
                     " of non-negative integers."));
    endif
  elseif (nargin > 3)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("nbinrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  ## Use 'size (ones (sz))' to ignore any trailing singleton dimensions in SZ
  if (! isscalar (r) && ! isequal (size (r), size (ones (sz))))
    error ("nbinrnd: R and PS must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (r, "single") || isa (ps, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Generate random sample from negative binomial distribution
  if (isscalar (r) && isscalar (ps))
    if ((r > 0) && (r < Inf) && (ps > 0) && (ps <= 1))
      rnd = randp ((1 - ps) ./ ps .* randg (r, sz, cls), cls);
    elseif ((r > 0) && (r < Inf) && (ps == 0))
      rnd = zeros (sz, cls);
    else
      rnd = NaN (sz, cls);
    endif
  else
    rnd = NaN (sz, cls);

    k = (r > 0) & (r < Inf) & (ps == 0);
    rnd(k) = 0;

    k = (r > 0) & (r < Inf) & (ps > 0) & (ps <= 1);
    rnd(k) = randp ((1 - ps(k)) ./ ps(k) .* randg (r(k), cls));
  endif

endfunction

## Test output
%!assert (size (nbinrnd (1, 0.5)), [1, 1])
%!assert (size (nbinrnd (1, 0.5 * ones (2, 1))), [2, 1])
%!assert (size (nbinrnd (1, 0.5 * ones (2, 2))), [2, 2])
%!assert (size (nbinrnd (ones (2, 1), 0.5)), [2, 1])
%!assert (size (nbinrnd (ones (2, 2), 0.5)), [2, 2])
%!assert (size (nbinrnd (1, 0.5, 3)), [3, 3])
%!assert (size (nbinrnd (1, 0.5, [4, 1])), [4, 1])
%!assert (size (nbinrnd (1, 0.5, 4, 1)), [4, 1])
%!assert (size (nbinrnd (1, 0.5, 4, 1, 5)), [4, 1, 5])
%!assert (size (nbinrnd (1, 0.5, 0, 1)), [0, 1])
%!assert (size (nbinrnd (1, 0.5, 1, 0)), [1, 0])
%!assert (size (nbinrnd (1, 0.5, 1, 2, 0, 5)), [1, 2, 0, 5])
%!assert (size (nbinrnd (1, 0.5, [])), [0, 0])
%!assert (size (nbinrnd (1, 0.5, [2, 0, 2, 1])), [2, 0, 2])

## Test class of input preserved
%!assert (class (nbinrnd (1, 0.5)), "double")
%!assert (class (nbinrnd (1, single (0.5))), "single")
%!assert (class (nbinrnd (1, single ([0.5, 0.5]))), "single")
%!assert (class (nbinrnd (single (1), 0.5)), "single")
%!assert (class (nbinrnd (single ([1, 1]), 0.5)), "single")

## Test input validation
%!error<nbinrnd: function called with too few input arguments.> nbinrnd ()
%!error<nbinrnd: function called with too few input arguments.> nbinrnd (1)
%!error<nbinrnd: R and PS must be of common size or scalars.> ...
%! nbinrnd (ones (3), ones (2))
%!error<nbinrnd: R and PS must be of common size or scalars.> ...
%! nbinrnd (ones (2), ones (3))
%!error<nbinrnd: R and PS must not be complex.> nbinrnd (i, 2, 3)
%!error<nbinrnd: R and PS must not be complex.> nbinrnd (1, i, 3)
%!error<nbinrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nbinrnd (1, 2, -1)
%!error<nbinrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nbinrnd (1, 2, 1.2)
%!error<nbinrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nbinrnd (1, 2, ones (2))
%!error<nbinrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nbinrnd (1, 2, [2 -1 2])
%!error<nbinrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! nbinrnd (1, 2, [2 0 2.5])
%!error<nbinrnd: dimensions must be non-negative integers.> ...
%! nbinrnd (1, 2, 2, -1, 5)
%!error<nbinrnd: dimensions must be non-negative integers.> ...
%! nbinrnd (1, 2, 2, 1.5, 5)
%!error<nbinrnd: R and PS must be scalars or of size SZ.> ...
%! nbinrnd (2, ones (2), 3)
%!error<nbinrnd: R and PS must be scalars or of size SZ.> ...
%! nbinrnd (2, ones (2), [3, 2])
%!error<nbinrnd: R and PS must be scalars or of size SZ.> ...
%! nbinrnd (2, ones (2), 3, 2)
