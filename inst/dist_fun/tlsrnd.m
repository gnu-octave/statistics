## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2022-2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{r} =} tlsrnd (@var{mu}, @var{sigma}, @var{nu})
## @deftypefnx {statistics} {@var{r} =} tlsrnd (@var{mu}, @var{sigma}, @var{nu}, @var{rows})
## @deftypefnx {statistics} {@var{r} =} tlsrnd (@var{mu}, @var{sigma}, @var{nu}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} {@var{r} =} tlsrnd (@var{mu}, @var{sigma}, @var{nu}, [@var{sz}])
##
## Random arrays from the location-scale Student's T distribution.
##
## Return a matrix of random samples from the location-scale Student's T
## distribution with location parameter @var{mu}, scale parameter @var{sigma},
## and @var{nu} degrees of freedom.
##
## @code{@var{r} = tlsrnd (@var{nu})} returns an array of random numbers chosen
## from the location-scale Student's T distribution with location parameter
## @var{mu}, scale parameter @var{sigma}, and @var{nu} degrees of freedom.  The
## size of @var{r} is the common size of @var{mu}, @var{sigma}, and @var{nu}.  A
## scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## When called with a single size argument, @code{tlsrnd} returns a square
## matrix with the dimension specified.  When called with more than one scalar
## argument, the first two arguments are taken as the number of rows and columns
## and any further arguments specify additional matrix dimensions.  The size may
## also be specified with a row vector of dimensions, @var{sz}.
##
## Further information about the location-scale Student's T distribution can be
## found at @url{https://en.wikipedia.org/wiki/Student%27s_t-distribution#Location-scale_t_distribution}
##
## @seealso{tlscdf, tlsinv, tlspdf, tlsfit, tlslike, tlsstat}
## @end deftypefn

function r = tlsrnd (mu, sigma, nu, varargin)

  ## Check for valid number of input arguments
  if (nargin < 3)
    error ("tlsrnd: function called with too few input arguments.");
  endif

  ## Check for common size of MU, SIGMA, and NU
  if (! isscalar (mu) || ! isscalar (sigma) || ! isscalar (nu))
    [retval, mu, sigma, nu] = common_size (mu, sigma, nu);
    if (retval > 0)
      error ("tlsrnd: MU, SIGMA, and NU must be of common size or scalars.");
    endif
  endif

  ## Check for NU being real
  if (iscomplex (mu) || iscomplex (sigma) || iscomplex (nu))
    error ("tlsrnd: MU, SIGMA, and NU must not be complex.");
  endif

  ## Parse and check SIZE arguments
  if (nargin == 3)
    sz = size (nu);
  elseif (nargin == 4)
    if (isscalar (varargin{1}) && varargin{1} >= 0 ...
                               && varargin{1} == fix (varargin{1}))
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0) ...
                                && all (varargin{1} == fix (varargin{1})))
      sz = varargin{1};
    elseif
      error (strcat (["tlsrnd: SZ must be a scalar or a row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 4)
    posint = cellfun (@(x) (! isscalar (x) || x < 0 || x != fix (x)), varargin);
    if (any (posint))
      error ("tlsrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (nu) && ! isequal (size (nu), sz))
    error ("tlsrnd: MU, SIGMA, and NU must be scalar or of size SZ.");
  endif

  ## Check for class type
  if (isa (mu, "single") || isa (sigma, "single") || isa (nu, "single"))
    cls = "single";
  else
    cls = "double";
  endif

  ## Call trnd to do the work
  r = mu + sigma .* trnd (nu, sz);

  ## Force class type
  r = cast (r, cls);

endfunction

## Test output
%!assert (size (tlsrnd (1, 2, 3)), [1, 1])
%!assert (size (tlsrnd (ones (2,1), 2, 3)), [2, 1])
%!assert (size (tlsrnd (ones (2,2), 2, 3)), [2, 2])
%!assert (size (tlsrnd (1, 2, 3, 3)), [3, 3])
%!assert (size (tlsrnd (1, 2, 3, [4 1])), [4, 1])
%!assert (size (tlsrnd (1, 2, 3, 4, 1)), [4, 1])
%!assert (size (tlsrnd (1, 2, 3, 4, 1)), [4, 1])
%!assert (size (tlsrnd (1, 2, 3, 4, 1, 5)), [4, 1, 5])
%!assert (size (tlsrnd (1, 2, 3, 0, 1)), [0, 1])
%!assert (size (tlsrnd (1, 2, 3, 1, 0)), [1, 0])
%!assert (size (tlsrnd (1, 2, 3, 1, 2, 0, 5)), [1, 2, 0, 5])
%!assert (tlsrnd (1, 2, 0, 1, 1), NaN)
%!assert (tlsrnd (1, 2, [0, 0, 0], [1, 3]), [NaN, NaN, NaN])

## Test class of input preserved
%!assert (class (tlsrnd (1, 2, 3)), "double")
%!assert (class (tlsrnd (single (1), 2, 3)), "single")
%!assert (class (tlsrnd (single ([1, 1]), 2, 3)), "single")
%!assert (class (tlsrnd (1, single (2), 3)), "single")
%!assert (class (tlsrnd (1, single ([2, 2]), 3)), "single")
%!assert (class (tlsrnd (1, 2, single (3))), "single")
%!assert (class (tlsrnd (1, 2, single ([3, 3]))), "single")

## Test input validation
%!error<tlsrnd: function called with too few input arguments.> tlsrnd ()
%!error<tlsrnd: function called with too few input arguments.> tlsrnd (1)
%!error<tlsrnd: function called with too few input arguments.> tlsrnd (1, 2)
%!error<tlsrnd: MU, SIGMA, and NU must be of common size or scalars.> ...
%! tlsrnd (ones (3), ones (2), 1)
%!error<tlsrnd: MU, SIGMA, and NU must be of common size or scalars.> ...
%! tlsrnd (ones (2), 1, ones (3))
%!error<tlsrnd: MU, SIGMA, and NU must be of common size or scalars.> ...
%! tlsrnd (1, ones (2), ones (3))
%!error<tlsrnd: MU, SIGMA, and NU must not be complex.> tlsrnd (i, 2, 3)
%!error<tlsrnd: MU, SIGMA, and NU must not be complex.> tlsrnd (1, i, 3)
%!error<tlsrnd: MU, SIGMA, and NU must not be complex.> tlsrnd (1, 2, i)
%!error<tlsrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! tlsrnd (1, 2, 3, -1)
%!error<tlsrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! tlsrnd (1, 2, 3, 1.2)
%!error<tlsrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! tlsrnd (1, 2, 3, ones (2))
%!error<tlsrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! tlsrnd (1, 2, 3, [2 -1 2])
%!error<tlsrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! tlsrnd (1, 2, 3, [2 0 2.5])
%!error<tlsrnd: SZ must be a scalar or a row vector of non-negative integers.> ...
%! tlsrnd (ones (2), 2, 3, ones (2))
%!error<tlsrnd: dimensions must be non-negative integers.> ...
%! tlsrnd (1, 2, 3, 2, -1, 5)
%!error<tlsrnd: dimensions must be non-negative integers.> ...
%! tlsrnd (1, 2, 3, 2, 1.5, 5)
%!error<tlsrnd: MU, SIGMA, and NU must be scalar or of size SZ.> ...
%! tlsrnd (ones (2,2), 2, 3, 3)
%!error<tlsrnd: MU, SIGMA, and NU must be scalar or of size SZ.> ...
%! tlsrnd (1, ones (2,2), 3, 3)
%!error<tlsrnd: MU, SIGMA, and NU must be scalar or of size SZ.> ...
%! tlsrnd (1, 2, ones (2,2), 3)
%!error<tlsrnd: MU, SIGMA, and NU must be scalar or of size SZ.> ...
%! tlsrnd (1, 2, ones (2,2), [3, 3])
%!error<tlsrnd: MU, SIGMA, and NU must be scalar or of size SZ.> ...
%! tlsrnd (1, 2, ones (2,2), 2, 3)
