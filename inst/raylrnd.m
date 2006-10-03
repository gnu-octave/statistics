## Copyright (C) 2006 Arno Onken
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

## -*- texinfo -*-
## @deftypefn {Function File} {@var{x} =} raylrnd (@var{sigma})
## @deftypefnx {Function File} {@var{x} =} raylrnd (@var{sigma}, @var{sz})
## @deftypefnx {Function File} {@var{x} =} raylrnd (@var{sigma}, @var{r}, @var{c})
## Returns a matrix of random samples from the Rayleigh distribution.
##
## Arguments are
##
## @itemize @bullet
## @item
## @var{sigma} is the parameter of the Rayleigh distribution. The elements
## of @var{sigma} must be positive.
##
## @item
## @var{sz} is the size of the matrix to be generated. @var{sz} must be a
## vector of non-negative integers.
##
## @item
## @var{r} is the number of rows of the matrix to be generated. @var{r} must
## be a non-negative integer.
##
## @item
## @var{c} is the number of columns of the matrix to be generated. @var{c}
## must be a non-negative integer.
## @end itemize
##
## Return values are
##
## @itemize @bullet
## @item
## @var{x} is a matrix of random samples from the Rayleigh distribution with
## corresponding parameter @var{sigma}. If neither @var{sz} nor @var{r} and
## @var{c} are specified, then @var{x} is of the same size as @var{sigma}.
## @end itemize
##
## Examples:
##
## @example
## sigma = 1:6;
## x = raylrnd (sigma)
##
## sz = [2, 3];
## x = raylrnd (0.5, sz)
##
## r = 2;
## c = 3;
## x = raylrnd (0.5, r, c)
## @end example
##
## References:
##
## @enumerate
## @item
## W. L. Martinez and A. R. Martinez. @cite{Computational Statistics
## Handbook with MATLAB.} Chapman & Hall/CRC, pages 547-557, 2001.
##
## @item
## Wikipedia contributors. Rayleigh distribution. @cite{Wikipedia, The Free
## Encyclopedia.}
## @uref{http://en.wikipedia.org/w/index.php?title=Rayleigh_distribution&oldid=69294908},
## August 2006.
## @end enumerate
## @end deftypefn

## Author: Arno Onken <whyly@gmx.net>
## Description: Random samples from the Rayleigh distribution

function x = raylrnd (sigma, r, c)

  # Check arguments
  if (nargin == 1)
    sz = size (sigma);
  elseif (nargin == 2)
    if (! isvector (r) || any ((r < 0) | round (r) != r))
      error ("raylrnd: sz must be a vector of non-negative integers")
    end
    sz = r(:)';
    if (! isscalar (sigma) && ! isempty (sigma) && (length (size (sigma)) != length (sz) || any (size (sigma) != sz)))
      error ("raylrnd: sigma must be scalar or of size sz");
    endif
  elseif (nargin == 3)
    if (! isscalar (r) || any ((r < 0) | round (r) != r))
      error ("raylrnd: r must be a non-negative integer")
    end
    if (! isscalar (c) || any ((c < 0) | round (c) != c))
      error ("raylrnd: c must be a non-negative integer")
    end
    sz = [r, c];
    if (! isscalar (sigma) && ! isempty (sigma) && (length (size (sigma)) != length (sz) || any (size (sigma) != sz)))
      error ("raylrnd: sigma must be scalar or of size [r, c]");
    endif
  else
    usage ("x = raylrnd (sigma [, sz |, r, c])");
  endif

  if (! isempty (sigma) && ! ismatrix (sigma))
    error ("raylrnd: sigma must be a numeric matrix");
  endif

  if (isempty (sigma))
    x = [];
  elseif (isscalar (sigma) && ! (sigma > 0))
    x = NaN .* ones (sz); 
  else
    # Draw random samples
    x = sqrt (-2 .* log (1 - rand (sz)) .* sigma .^ 2);

    # Continue argument check
    k = find (! (sigma > 0));
    if (any (k))
      x (k) = NaN;
    endif
  endif

endfunction

%!test
%! sigma = 1:6;
%! x = raylrnd (sigma);
%! assert (size (x), size (sigma));
%! assert (all (x >= 0));

%!test
%! sigma = 0.5;
%! sz = [2, 3];
%! x = raylrnd (sigma, sz);
%! assert (size (x), sz);
%! assert (all (x >= 0));

%!test
%! sigma = 0.5;
%! r = 2;
%! c = 3;
%! x = raylrnd (sigma, r, c);
%! assert (size (x), [r c]);
%! assert (all (x >= 0));
