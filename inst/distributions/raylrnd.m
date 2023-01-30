## Copyright (C) 2006, 2007 Arno Onken <asnelt@asnelt.org>
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
## @deftypefn  {statistics} @var{r} = raylrnd (@var{sigma})
## @deftypefnx {statistics} @var{r} = raylrnd (@var{sigma}, @var{rows})
## @deftypefnx {statistics} @var{r} = raylrnd (@var{sigma}, @var{rows}, @var{cols}, @dots{})
## @deftypefnx {statistics} @var{r} = raylrnd (@var{sigma}, [@var{sz}])
##
## Random arrays from the Rayleigh distribution.
##
## @code{@var{r} = raylrnd (@var{sigma})} returns an array of random numbers
## chosen from the Rayleigh distribution with scale parameter @var{sigma}.  The
## size of @var{r} is the size of @var{sigma}.  @var{sigma} must be a finite
## real number greater than 0, otherwise NaN is returned.
##
## When called with a single size argument, return a square matrix with
## the dimension specified.  When called with more than one scalar argument the
## first two arguments are taken as the number of rows and columns and any
## further arguments specify additional matrix dimensions.  The size may also
## be specified with a vector of dimensions @var{sz}.
##
## @subheading References
##
## @enumerate
## @item
## Wendy L. Martinez and Angel R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB}. Appendix E, pages 547-557, Chapman & Hall/CRC,
## 2001.
##
## @item
## Athanasios Papoulis. @cite{Probability, Random Variables, and Stochastic
## Processes}. pages 104 and 148, McGraw-Hill, New York, second edition,
## 1984.
## @end enumerate
##
## @seealso{raylcdf, raylinv, raylrnd, raylstat}
## @end deftypefn

function r = raylrnd (sigma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1)
    print_us
  endif

  ## Check for SIGMA being real
  if (iscomplex (sigma))
    error ("raylrnd: SIGMA must not be complex.");
  endif

  ## Check for SIZE vector or DIMENSION input arguments
  if (nargin == 1)
    sz = size (sigma);
  elseif (nargin == 2)
    if (isscalar (varargin{1}) && varargin{1} >= 0)
      sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
      sz = varargin{1};
    else
      error (strcat (["raylrnd: dimension vector must be row vector"], ...
                     [" of non-negative integers."]));
    endif
  elseif (nargin > 2)
    if (any (cellfun (@(x) (! isscalar (x) || x < 0), varargin)))
      error ("raylrnd: dimensions must be non-negative integers.");
    endif
    sz = [varargin{:}];
  endif

  ## Check that parameters match requested dimensions in size
  if (! isscalar (sigma) && ! isequal (size (sigma), sz))
    error ("raylrnd: SIGMA must be scalar or of size SZ.");
  endif

  ## Generate random sample from Rayleigh distribution
  r = sqrt (-2 .* log (1 - rand (sz)) .* sigma .^ 2);

  ## Check for valid parameter
  k = find (! (sigma > 0));
  if (any (k))
    r(k) = NaN;
  endif

  ## Cast into appropriate class
  if (isa (sigma, "single"))
    r = cast (r, "single");
  endif

endfunction

%!test
%! sigma = 1:6;
%! r = raylrnd (sigma);
%! assert (size (r), size (sigma));
%! assert (all (r >= 0));

%!test
%! sigma = 0.5;
%! sz = [2, 3];
%! r = raylrnd (sigma, sz);
%! assert (size (r), sz);
%! assert (all (r >= 0));

%!test
%! sigma = 0.5;
%! rows = 2;
%! cols = 3;
%! r = raylrnd (sigma, rows, cols);
%! assert (size (r), [rows, cols]);
%! assert (all (r >= 0));


## Test input validation
%!error poissrnd ()
%!error poissrnd (1, -1)
%!error poissrnd (1, ones (2))
%!error poissrnd (1, 2, ones (2))
%!error poissrnd (i)
%!error poissrnd (1, 2, -1)
%!error poissrnd (1, [2 -1 2])
%!error poissrnd (ones (2,2), 3)
%!error poissrnd (ones (2,2), [3, 2])
%!error poissrnd (ones (2,2), 2, 3)
