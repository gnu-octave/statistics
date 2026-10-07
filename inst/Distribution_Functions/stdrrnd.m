## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{r} =} stdrrnd (@var{k}, @var{df})
## @deftypefnx {statistics} {@var{r} =} stdrrnd (@var{k}, @var{df}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} stdrrnd (@var{k}, @var{df}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} stdrrnd (@var{k}, @var{df}, [@var{sz}])
##
## Random arrays from the studentized range distribution.
##
## @code{@var{r} = stdrrnd (@var{k}, @var{df})} returns an array of random
## numbers chosen from the studentized range distribution for @var{k} groups
## and @var{df} degrees of freedom: the range of @var{k} standard normal
## variables divided by an independent @math{sqrt (chi2 / df)}.  The size of
## @var{r} is the common size of @var{k} and @var{df}.  A scalar input functions
## as a constant matrix of the same size as the other input.  @var{k} must be an
## integer of at least 2 and @var{df} positive, @code{Inf} included; otherwise
## @code{NaN} is returned.
##
## When called with a single size argument, @code{stdrrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## MATLAB has no public counterpart of this function.
##
## Further information about the studentized range distribution can be found
## at @url{https://en.wikipedia.org/wiki/Studentized_range_distribution}
##
## @seealso{stdrcdf, stdrinv, stdrpdf, stdrstat}
## @end deftypefn

function r = stdrrnd (k, df, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("stdrrnd: function called with too few input arguments.");
  endif

  ## Check for common size of K and DF
  if (! isscalar (k) || ! isscalar (df))
    [retval, k, df] = common_size (k, df);
    if (retval > 0)
      error ("stdrrnd: K and DF must be of common size or scalars.");
    endif
  endif

  ## Check for K and DF being reals
  if (iscomplex (k) || iscomplex (df))
    error ("stdrrnd: K and DF must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 2)
    sz = size (k);
  elseif (nargin == 3)
    if (isscalar (varargin{1}) && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif (isempty (varargin{1}))
      r = [];
      return;
    else
      error (strcat ("stdrrnd: SZ must be a scalar or a row vector", ...
                     " of integers."));
    endif
  elseif (nargin > 3)
    notint = cellfun (@(x) (! isscalar (x) || x != fix (x)), varargin);
    if (any (notint))
      error ("stdrrnd: dimensions must be integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Negative dimensions are treated as zero, as in core Octave and MATLAB
  sz = max (sz, 0);

  ## Check that parameters match requested dimensions in size
  ## Use 'size (ones (sz))' to ignore any trailing singleton dimensions in SZ
  if (! isscalar (k) && ! isequal (size (k), size (ones (sz))))
    error ("stdrrnd: K and DF must be scalars or of size SZ.");
  endif

  ## Check for class type
  if (isa (k, 'single') || isa (df, 'single'))
    cls = 'single';
  else
    cls = 'double';
  endif

  ## One draw per element: the range of k standard normals over the root of
  ## a chi-square on df degrees of freedom divided by df
  n = prod (sz);
  k = double (k(:)) .* ones (n, 1);
  df = double (df(:)) .* ones (n, 1);
  r = NaN (n, 1);
  valid = (k >= 2) & (k == fix (k)) & isfinite (k) & (df > 0);
  for kk = unique (k(valid))'
    idx = find (valid & k == kk);
    z = randn (kk, numel (idx));
    w = max (z, [], 1)' - min (z, [], 1)';
    s = ones (numel (idx), 1);
    fin = isfinite (df(idx));
    s(fin) = sqrt (2 * randg (df(idx(fin)) / 2) ./ df(idx(fin)));
    r(idx) = w ./ s;
  endfor
  r = cast (reshape (r, sz), cls);

endfunction

%!demo
%! ## Compare random samples with the density they are drawn from
%! rng (42);
%! r = stdrrnd (4, 10, 1, 10000);
%! x = 0:0.05:10;
%! hist (r, x, 1 / 0.05)
%! hold on
%! plot (x, stdrpdf (x, 4, 10), '-r', 'linewidth', 2)
%! hold off
%! xlim ([0, 10])
%! legend ({'10000 samples', 'k = 4, df = 10'}, 'location', 'northeast')
%! title ('Studentized range random samples')
%! xlabel ('values in r')
%! ylabel ('density')

## Test output
%!assert_equal (size (stdrrnd (3, 10)), [1, 1])
%!assert_equal (size (stdrrnd (3 * ones (2, 1), 10)), [2, 1])
%!assert_equal (size (stdrrnd (3, 10 * ones (2, 2))), [2, 2])
%!assert_equal (size (stdrrnd (3, 10, 3)), [3, 3])
%!assert_equal (size (stdrrnd (3, 10, [4, 1])), [4, 1])
%!assert_equal (size (stdrrnd (3, 10, 4, 1)), [4, 1])
%!assert_equal (size (stdrrnd (3, 10, 4, 1, 5)), [4, 1, 5])
%!assert_equal (size (stdrrnd (3, 10, 0, 1)), [0, 1])
%!assert_equal (size (stdrrnd (3, 10, [])), [0, 0])
%!assert_equal (size (stdrrnd (3, 10, [2, -1, 2])), [2, 0, 2])
%!assert_equal (stdrrnd ([1, 2.5, NaN], 10), [NaN, NaN, NaN])
%!assert_equal (stdrrnd (3, [0, -1, NaN]), [NaN, NaN, NaN])
%!assert_equal (all (stdrrnd (3, Inf, 1, 100) > 0), true)

## Test class of input preserved
%!assert_equal (class (stdrrnd (3, 10)), 'double')
%!assert_equal (class (stdrrnd (single (3), 10)), 'single')
%!assert_equal (class (stdrrnd (3, single ([10, 10]))), 'single')

## Test input validation
%!error<stdrrnd: function called with too few input arguments.> stdrrnd ()
%!error<stdrrnd: function called with too few input arguments.> stdrrnd (3)
%!error<stdrrnd: K and DF must be of common size or scalars.> ...
%! stdrrnd (ones (3), ones (2))
%!error<stdrrnd: K and DF must not be complex.> stdrrnd (i, 10)
%!error<stdrrnd: K and DF must not be complex.> stdrrnd (3, i)
%!error<stdrrnd: SZ must be a scalar or a row vector of integers.> ...
%! stdrrnd (3, 10, 1.2)
%!error<stdrrnd: SZ must be a scalar or a row vector of integers.> ...
%! stdrrnd (3, 10, ones (2))
%!error<stdrrnd: dimensions must be integers.> stdrrnd (3, 10, 2, 1.5)
%!error<stdrrnd: K and DF must be scalars or of size SZ.> ...
%! stdrrnd (3 * ones (2), 10, 3)
