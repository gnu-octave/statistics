## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{p} =} gumbelcdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} gumbelcdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} gumbelcdf (@var{x}, @var{mu}, @var{beta})
## @deftypefnx {statistics} {@var{p} =} gumbelcdf (@dots{}, @qcode{"upper"})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} gumbelcdf (@var{x}, @var{mu}, @var{beta}, @var{pcov})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} gumbelcdf (@var{x}, @var{mu}, @var{beta}, @var{pcov}, @var{alpha})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} gumbelcdf (@dots{}, @qcode{"upper"})
##
## Gumbel cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the Gumbel distribution (also known as the extreme value or the type
## I generalized extreme value distribution) with location parameter @var{mu}
## and scale parameter @var{beta}.  The size of @var{p} is the common size of
## @var{x}, @var{mu} and @var{beta}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Default values are @var{mu} = 0, @var{beta} = 1.
##
## When called with three output arguments, i.e. @code{[@var{p}, @var{plo},
## @var{pup}]}, it computes the confidence bounds for @var{p} when the input
## parameters @var{mu} and @var{sigma} are estimates.  In such case, @var{pcov},
## a 2-by-2 matrix containing the covariance matrix of the estimated parameters,
## is necessary.  Optionally, @var{alpha}, which has a default value of 0.05,
## specifies the @qcode{100 * (1 - @var{alpha})} percent confidence bounds.
## @var{plo} and @var{pup} are arrays of the same size as @var{p} containing the
## lower and upper confidence bounds.
##
## @code{[@dots{}] = gumbelcdf (@dots{}, "upper")} computes the upper tail
## probability of the Gumbel distribution with parameters @var{mu} and
## @var{beta}, at the values in @var{x}.
##
## The Gumbel distribution is used to model the distribution of the maximum (or
## the minimum) of a number of samples of various distributions.  This version
## is suitable for modeling maxima.  For modeling minima, use the alternative
## extreme value CDF, @code{evcdf}.
##
## @code{[@dots{}] = gumbelcdf (@dots{}, "upper")} computes the upper tail
## probability of the extreme value (Gumbel) distribution.
##
## Further information about the Gumbel distribution can be found at
## @url{https://en.wikipedia.org/wiki/Gumbel_distribution}
##
## @seealso{gumbelinv, gumbelpdf, gumbelrnd, gumbelfit, gumbellike, gumbelstat,
## evcdf}
## @end deftypefn

function [varargout] = gumbelcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 6)
    error ("gumbelcdf: invalid number of input arguments.");
  endif

  ## Check for 'upper' flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1 && ischar (varargin{end}) && ...
          ! strcmpi (varargin{end}, "upper"))
    error ("gumbelcdf: invalid argument for upper tail.");
  else
    uflag = false;
  endif

  ## Get extra arguments (if they exist) or add defaults
  if (numel (varargin) > 0)
    mu = varargin{1};
  else
    mu = 0;
  endif
  if (numel (varargin) > 1)
    beta = varargin{2};
  else
    beta = 1;
  endif
  if (numel (varargin) > 2)
    pcov = varargin{3};
    ## Check for valid covariance matrix 2x2
    if (! isequal (size (pcov), [2, 2]))
      error ("gumbelcdf: invalid size of covariance matrix.");
    endif
  else
    ## Check that cov matrix is provided if 3 output arguments are requested
    if (nargout > 1)
      error ("gumbelcdf: covariance matrix is required for confidence bounds.");
    endif
    pcov = [];
  endif
  if (numel (varargin) > 3)
    alpha = varargin{4};
    ## Check for valid alpha value
    if (! isnumeric (alpha) || numel (alpha) !=1 || alpha <= 0 || alpha >= 1)
      error ("gumbelcdf: invalid value for alpha.");
    endif
  else
    alpha = 0.05;
  endif

  ## Check for common size of X, MU, and BETA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (beta))
    [err, x, mu, beta] = common_size (x, mu, beta);
    if (err > 0)
      error ("gumbelcdf: X, MU, and BETA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and BETA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (beta))
    error ("gumbelcdf: X, MU, and BETA must not be complex.");
  endif

  ## Return NaNs for out of range parameters.
  beta(beta <= 0) = NaN;
  ## Compute extreme value cdf
  z = (x - mu) ./ beta;
  if (uflag)
    p = -expm1 (-exp (-z));
  else
    p = exp (-exp (-z));
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (beta, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  varargout{1} = cast (p, is_class);
  if (nargout > 1)
    plo = NaN (size (z), is_class);
    pup = NaN (size (z), is_class);
  endif

  ## Check beta
  if (isscalar (beta))
    if (beta > 0)
      sigma_p = true (size (z));
    else
      if (nargout == 3)
        varargout{2} = plo;
        varargout{3} = pup;
      endif
      return;
    endif
  else
    sigma_p = beta > 0;
  endif

  ## Compute confidence bounds (if requested)
  if (nargout >= 2)
   zvar = (pcov(1,1) + 2 * pcov(1,2) * z(sigma_p) + ...
           pcov(2,2) * z(sigma_p) .^ 2) ./ (beta .^ 2);
   if (any (zvar < 0))
      error ("gumbelcdf: bad covariance matrix.");
   endif
   normz = -norminv (alpha / 2);
   halfwidth = normz * sqrt (zvar);
   zlo = z(sigma_p) - halfwidth;
   zup = z(sigma_p) + halfwidth;
   if (uflag)
     plo(sigma_p) = -expm1 (-exp (-zup));
     pup(sigma_p) = -expm1 (-exp (-zlo));
   else
     plo(sigma_p) = exp (-exp (-zlo));
     pup(sigma_p) = exp (-exp (-zup));
   endif
   varargout{2} = plo;
   varargout{3} = pup;
  endif

endfunction

%!demo
%! ## Plot various CDFs from the Gumbel distribution
%! x = -5:0.01:20;
%! p1 = gumbelcdf (x, 0.5, 2);
%! p2 = gumbelcdf (x, 1.0, 2);
%! p3 = gumbelcdf (x, 1.5, 3);
%! p4 = gumbelcdf (x, 3.0, 4);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c")
%! grid on
%! legend ({"μ = 0.5, β = 2", "μ = 1.0, β = 2", ...
%!          "μ = 1.5, β = 3", "μ = 3.0, β = 4"}, "location", "southeast")
%! title ("Gumbel CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-Inf, 1, 2, Inf];
%! y = [0, 0.3679, 0.6922, 1];
%!assert (gumbelcdf (x, ones (1,4), ones (1,4)), y, 1e-4)
%!assert (gumbelcdf (x, 1, ones (1,4)), y, 1e-4)
%!assert (gumbelcdf (x, ones (1,4), 1), y, 1e-4)
%!assert (gumbelcdf (x, [0, -Inf, NaN, Inf], 1), [0, 1, NaN, NaN], 1e-4)
%!assert (gumbelcdf (x, 1, [Inf, NaN, -1, 0]), [NaN, NaN, NaN, NaN], 1e-4)
%!assert (gumbelcdf ([x(1:2), NaN, x(4)], 1, 1), [y(1:2), NaN, y(4)], 1e-4)
%!assert (gumbelcdf (x, "upper"), [1, 0.3078, 0.1266, 0], 1e-4)

## Test class of input preserved
%!assert (gumbelcdf ([x, NaN], 1, 1), [y, NaN], 1e-4)
%!assert (gumbelcdf (single ([x, NaN]), 1, 1), single ([y, NaN]), 1e-4)
%!assert (gumbelcdf ([x, NaN], single (1), 1), single ([y, NaN]), 1e-4)
%!assert (gumbelcdf ([x, NaN], 1, single (1)), single ([y, NaN]), 1e-4)

## Test input validation
%!error<gumbelcdf: invalid number of input arguments.> gumbelcdf ()
%!error<gumbelcdf: invalid number of input arguments.> gumbelcdf (1,2,3,4,5,6,7)
%!error<gumbelcdf: invalid argument for upper tail.> gumbelcdf (1, 2, 3, 4, "uper")
%!error<gumbelcdf: X, MU, and BETA must be of common size or scalars.> ...
%! gumbelcdf (ones (3), ones (2), ones (2))
%!error<gumbelcdf: invalid size of covariance matrix.> gumbelcdf (2, 3, 4, [1, 2])
%!error<gumbelcdf: covariance matrix is required for confidence bounds.> ...
%! [p, plo, pup] = gumbelcdf (1, 2, 3)
%!error<gumbelcdf: invalid value for alpha.> [p, plo, pup] = ...
%! gumbelcdf (1, 2, 3, [1, 0; 0, 1], 0)
%!error<gumbelcdf: invalid value for alpha.> [p, plo, pup] = ...
%! gumbelcdf (1, 2, 3, [1, 0; 0, 1], 1.22)
%!error<gumbelcdf: invalid value for alpha.> [p, plo, pup] = ...
%! gumbelcdf (1, 2, 3, [1, 0; 0, 1], "alpha", "upper")
%!error<gumbelcdf: X, MU, and BETA must not be complex.> gumbelcdf (i, 2, 2)
%!error<gumbelcdf: X, MU, and BETA must not be complex.> gumbelcdf (2, i, 2)
%!error<gumbelcdf: X, MU, and BETA must not be complex.> gumbelcdf (2, 2, i)
%!error<gumbelcdf: bad covariance matrix.> ...
%! [p, plo, pup] = gumbelcdf (1, 2, 3, [1, 0; 0, -inf], 0.04)
