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
## @deftypefn  {statistics} {@var{p} =} raylcdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} raylcdf (@var{x}, @var{sigma})
## @deftypefnx {statistics} {@var{p} =} raylcdf (@dots{}, "upper")
##
## Rayleigh cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the lognormal distribution with scale parameter
## @var{sigma}.  The size of @var{p} is the common size of @var{x} and
## @var{sigma}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.
##
## Default value is @var{sigma} = 1.
##
## @code{[@dots{}] = logncdf (@dots{}, "upper")} computes the upper tail
## probability of the lognormal distribution.
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
## @seealso{raylinv, raylpdf, raylrnd, raylstat}
## @end deftypefn

function p = raylcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 3)
    error ("raylcdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1  && ischar (varargin{end}) && ...
                         ! strcmpi (varargin{end}, "upper"))
    error ("raylcdf: invalid argument for upper tail.");
  elseif (nargin > 2  && ! strcmpi (varargin{end}, "upper"))
    error ("raylcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Get extra arguments (if they exist) or add defaults
  if (numel (varargin) > 0)
    sigma = varargin{1};
  else
    sigma = 1;
  endif

  ## Check for common size of X and SIGMA
  if (! isscalar (x) || ! isscalar (sigma))
    [retval, x, sigma] = common_size (x, sigma);
    if (retval > 0)
      error ("raylcdf: X and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X and SIGMA being reals
  if (iscomplex (x) || iscomplex (sigma))
    error ("raylcdf: X and SIGMA must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (sigma, "single"));
    p = zeros (size (x), "single");
  else
    p = zeros (size (x));
  endif

  ## Force 1 for upper flag and X <= 0
  k0 = sigma > 0 & x <= 0;
  if (uflag && any (k0(:)))
    p(k0) = 1;
  end

  ## Calculate Rayleigh CDF for valid parameter and data range
  k = sigma > 0 & x > 0;
  if (any (k(:)))
    if (uflag)
        p(k) = exp (-x(k) .^ 2 ./ (2 * sigma(k) .^ 2));
    else
        p(k) = - expm1 (-x(k) .^ 2 ./ (2 * sigma(k) .^ 2));
    endif
  endif

  ## Continue argument check
  p(! (k0 | k)) = NaN;

endfunction

%!test
%! x = 0:0.5:2.5;
%! sigma = 1:6;
%! p = raylcdf (x, sigma);
%! expected_p = [0.0000, 0.0308, 0.0540, 0.0679, 0.0769, 0.0831];
%! assert (p, expected_p, 0.001);

%!test
%! x = 0:0.5:2.5;
%! p = raylcdf (x, 0.5);
%! expected_p = [0.0000, 0.3935, 0.8647, 0.9889, 0.9997, 1.0000];
%! assert (p, expected_p, 0.001);

%!shared x, p
%! x = [-1, 0, 1, 2, Inf];
%! p = [0, 0, 0.39346934028737, 0.86466471676338, 1];
%!assert (raylcdf (x), p, 1e-14)
%!assert (raylcdf (x, 1), p, 1e-14)
%!assert (raylcdf (x, "upper"), 1 - p, 1e-14)

## Test input validation
%!error<raylcdf: invalid number of input arguments.> raylcdf ()
%!error<raylcdf: invalid number of input arguments.> raylcdf (1, 2, 3, 4)
%!error<raylcdf: invalid argument for upper tail.> raylcdf (1, 2, "uper")
%!error<raylcdf: invalid argument for upper tail.> raylcdf (1, 2, 3)
%!error<raylcdf: X and SIGMA must be of common size or scalars.> ...
%! raylcdf (ones (3), ones (2))
%!error<raylcdf: X and SIGMA must not be complex.> raylcdf (i, 2)
%!error<raylcdf: X and SIGMA must not be complex.> raylcdf (2, i)
