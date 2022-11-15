## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} @var{r} = evrnd (@var{mu}, @var{sigma})
## @deftypefnx {Function File} @var{r} = evrnd (@var{mu}, @var{sigma}, @var{m}, @var{n}, @dots{})
## @deftypefnx {Function File} @var{r} = evrnd (@var{mu}, @var{sigma}, [@var{m}, @var{n}, @dots{}])
##
## Random arrays from the extreme value distribution.
##
## @code{@var{r} = evrnd (@var{mu}, @var{sigma})} returns an array of random
## numbers chosen from the type 1 extreme value distribution with location
## parameter @var{mu} and scale parameter @var{sigma}.  The size of @var{r} is
## the common size of @var{mu} and @var{sigma}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## @code{@var{r} = evrnd (@var{mu}, @var{sigma}, @var{m}, @var{n}, @dots{})} or
## @code{@var{r} = evrnd (@var{mu}, @var{sigma}, [@var{m}, @var{n}, @dots{}])}
## returns an M-by-N-by-... array.
##
## The type 1 extreme value distribution is also known as the Gumbel
## distribution.  The version used here is suitable for modeling minima; the
## mirror image of this distribution can be used to model maxima by negating
## @var{x}.  If @var{y} has a Weibull distribution, then
## @code{@var{x} = log (@var{y})} has the type 1 extreme value distribution.
##
## @seealso{evcdf, evinv, evpdf, evfit, evlike, evstat}
## @end deftypefn

function r = evrnd (mu, sigma, varargin)

  ## Check for valid number of input arguments
  if (nargin < 2)
    error ("evrnd: too few input arguments.");
  endif
  ## Check for appropriate class
  if (isa (mu, "single") || isa (sigma, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif
  ## Check for additional dimensions in varargin and get their size
  dim_vec = 1;
  if (nargin > 2)
    extra_varargin = numel (varargin(:));
    if (extra_varargin == 1)
      size_dim = varargin{1};
      ## Check for empty input argument
      if (isempty (size_dim))
        error ("evrnd: extra argument for size of output array cannot be empty.");
      endif
      dim_vec = zeros (size_dim, is_class);
    elseif (extra_varargin > 1)
      for i = 1:extra_varargin
        size_dim(i) = varargin{i};
      endfor
      dim_vec = zeros (size_dim, is_class);
    endif
  endif
  ## Check for common size of MU, SIGMA, and output based on given dimensions
  if (! isscalar (mu) || ! isscalar (sigma) || ! isscalar (dim_vec))
    [err, mu, sigma, dim_vec] = common_size (mu, sigma, dim_vec);
    if (err > 0)
      error ("evrnd: MU, SIGMA, and DIM vector must be of common size or scalars.");
    endif
  endif
  ## Get final dimensions of returning random array
  size_out = size (mu);
  ## Return NaNs for out of range values of SIGMA
  sigma(sigma < 0) = NaN;
  ## Generate uniform random values, and apply the extreme value inverse CDF.
  r = log (-log (rand (size_out, is_class))) .* sigma + mu;

endfunction

## Test input validation
%!error<evrnd: too few input arguments.> evrnd ()
%!error<evrnd: extra argument for size of output array cannot be empty.> ...
%! evrnd (ones (3), ones (2), [])
%!error<evrnd: MU, SIGMA, and DIM vector must be of common size or scalars.> ...
%! evrnd (ones (3), ones (2))
%!error<evrnd: MU, SIGMA, and DIM vector must be of common size or scalars.> ...
%! evrnd (ones (2), ones (2), 3, 2)
%!error<evrnd: MU, SIGMA, and DIM vector must be of common size or scalars.> ...
%! evrnd (ones (2), ones (2), 1, 2)
