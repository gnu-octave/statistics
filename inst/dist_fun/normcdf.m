## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2012 Rik Wehbring
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
## @deftypefn  {statistics} {@var{p} =} normcdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} normcdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} normcdf (@var{x}, @var{mu}, @var{sigma})
## @deftypefnx {statistics} {@var{p} =} normcdf (@dots{}, @qcode{"upper"})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} normcdf (@var{x}, @var{mu}, @var{sigma}, @var{pcov})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} normcdf (@var{x}, @var{mu}, @var{sigma}, @var{pcov}, @var{alpha})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} normcdf (@dots{}, @qcode{"upper"})
##
## Normal cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the normal distribution with mean @var{mu} and standard deviation
## @var{sigma}.  The size of @var{p} is the common size of @var{x}, @var{mu} and
## @var{sigma}.  A scalar input functions as a constant matrix of the same size
## as the other inputs.
##
## Default values are @var{mu} = 0, @var{sigma} = 1.
##
## When called with three output arguments, i.e. @qcode{[@var{p}, @var{plo},
## @var{pup}]}, @code{normcdf} computes the confidence bounds for @var{p} when
## the input parameters @var{mu} and @var{sigma} are estimates.  In such case,
## @var{pcov}, a @math{2x2} matrix containing the covariance matrix of the
## estimated parameters, is necessary.  Optionally, @var{alpha}, which has a
## default value of 0.05, specifies the @qcode{100 * (1 - @var{alpha})} percent
## confidence bounds.  @var{plo} and @var{pup} are arrays of the same size as
## @var{p} containing the lower and upper confidence bounds.
##
## @code{[@dots{}] = normcdf (@dots{}, "upper")} computes the upper tail
## probability of the normal distribution with parameters @var{mu} and
## @var{sigma}, at the values in @var{x}.  This can be used to compute a
## right-tailed p-value.  To compute a two-tailed p-value, use
## @code{2 * normcdf (-abs (@var{x}), @var{mu}, @var{sigma})}.
##
## Further information about the normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Normal_distribution}
##
## @seealso{norminv, normpdf, normrnd, normfit, normlike, normstat}
## @end deftypefn

function [varargout] = normcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 6)
    error ("normcdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1  && ischar (varargin{end}) && ...
          ! strcmpi (varargin{end}, "upper"))
    error ("normcdf: invalid argument for upper tail.");
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
    sigma = varargin{2};
  else
    sigma = 1;
  endif
  if (numel (varargin) > 2)
    pcov = varargin{3};
    ## Check for valid covariance matrix 2x2
    if (! isequal (size (pcov), [2, 2]))
      error ("normcdf: invalid size of covariance matrix.");
    endif
  else
    ## Check that cov matrix is provided if 3 output arguments are requested
    if (nargout > 1)
      error ("normcdf: covariance matrix is required for confidence bounds.");
    endif
    pcov = [];
  endif
  if (numel (varargin) > 3)
    alpha = varargin{4};
    ## Check for valid alpha value
    if (! isnumeric (alpha) || numel (alpha) !=1 || alpha <= 0 || alpha >= 1)
      error ("normcdf: invalid value for alpha.");
   end
  else
    alpha = 0.05;
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (mu) || ! isscalar (sigma))
    [err, x, mu, sigma] = common_size (x, mu, sigma);
    if (err > 0)
      error ("normcdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("normcdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Compute normal CDF
  z = (x - mu) ./ sigma;
  if (uflag)
    z = -z;
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  p = NaN (size (z), is_class);
  if (nargout > 1)
    plo = NaN (size (z), is_class);
    pup = NaN (size (z), is_class);
  endif

  ## Check SIGMA
  if (isscalar (sigma))
    if (sigma > 0)
      sigma_p = true (size (z));
      sigma_z = false (size (z));
    elseif (sigma == 0)
      sigma_z = true (size (z));
      sigma_p = false (size (z));
    else
      if (nargout <= 1)
        varargout{1} = p;
      elseif (nargout == 3)
        varargout{1} = p;
        varargout{2} = plo;
        varargout{3} = pup;
      endif
      return;
    endif
  else
    sigma_p = sigma > 0;
    sigma_z = sigma == 0;
  endif

  ## Set edge cases when SIGMA = 0
  if (uflag)
    p(sigma_z & x < mu) = 1;
    p(sigma_z & x >= mu) = 0;
    if (nargout > 1)
      plo(sigma_z & x < mu) = 1;
      plo(sigma_z & x >= mu) = 0;
      pup(sigma_z & x < mu) = 1;
      pup(sigma_z & x >= mu) = 0;
    endif
  else
    p(sigma_z & x < mu) = 0;
    p(sigma_z & x >= mu) = 1;
    if (nargout >= 2)
      plo(sigma_z & x < mu) = 0;
      plo(sigma_z & x >= mu) = 1;
      pup(sigma_z & x < mu) = 0;
      pup(sigma_z & x >= mu) = 1;
    endif
  endif

  ## Compute cases when SIGMA > 0
  p(sigma_p) = 0.5 * erfc (-z(sigma_p) ./ sqrt (2));
  varargout{1} = p;

  ## Compute confidence bounds (if requested)
  if (nargout >= 2)
    zvar = (pcov(1,1) + 2 * pcov(1,2) * z(sigma_p) + ...
           pcov(2,2) * z(sigma_p) .^ 2) ./ (sigma .^ 2);
    if (any (zvar < 0))
      error ("normcdf: bad covariance matrix.");
    endif
    normz = -norminv (alpha / 2);
    halfwidth = normz * sqrt (zvar);
    zlo = z(sigma_p) - halfwidth;
    zup = z(sigma_p) + halfwidth;
    plo(sigma_p) = 0.5 * erfc (-zlo ./ sqrt (2));
    pup(sigma_p) = 0.5 * erfc (-zup ./ sqrt (2));
    varargout{2} = plo;
    varargout{3} = pup;
  endif

endfunction

%!demo
%! ## Plot various CDFs from the normal distribution
%! x = -5:0.01:5;
%! p1 = normcdf (x, 0, 0.5);
%! p2 = normcdf (x, 0, 1);
%! p3 = normcdf (x, 0, 2);
%! p4 = normcdf (x, -2, 0.8);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r", x, p4, "-c")
%! grid on
%! xlim ([-5, 5])
%! legend ({"μ = 0, σ = 0.5", "μ = 0, σ = 1", ...
%!          "μ = 0, σ = 2", "μ = -2, σ = 0.8"}, "location", "southeast")
%! title ("Normal CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-Inf 1 2 Inf];
%! y = [0, 0.5, 1/2*(1+erf(1/sqrt(2))), 1];
%!assert (normcdf (x, ones (1,4), ones (1,4)), y)
%!assert (normcdf (x, 1, ones (1,4)), y)
%!assert (normcdf (x, ones (1,4), 1), y)
%!assert (normcdf (x, [0, -Inf, NaN, Inf], 1), [0, 1, NaN, NaN])
%!assert (normcdf (x, 1, [Inf, NaN, -1, 0]), [NaN, NaN, NaN, 1])
%!assert (normcdf ([x(1:2), NaN, x(4)], 1, 1), [y(1:2), NaN, y(4)])
%!assert (normcdf (x, "upper"), [1, 0.1587, 0.0228, 0], 1e-4)

## Test class of input preserved
%!assert (normcdf ([x, NaN], 1, 1), [y, NaN])
%!assert (normcdf (single ([x, NaN]), 1, 1), single ([y, NaN]), eps ("single"))
%!assert (normcdf ([x, NaN], single (1), 1), single ([y, NaN]), eps ("single"))
%!assert (normcdf ([x, NaN], 1, single (1)), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<normcdf: invalid number of input arguments.> normcdf ()
%!error<normcdf: invalid number of input arguments.> normcdf (1,2,3,4,5,6,7)
%!error<normcdf: invalid argument for upper tail.> normcdf (1, 2, 3, 4, "uper")
%!error<normcdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! normcdf (ones (3), ones (2), ones (2))
%!error<normcdf: invalid size of covariance matrix.> normcdf (2, 3, 4, [1, 2])
%!error<normcdf: covariance matrix is required for confidence bounds.> ...
%! [p, plo, pup] = normcdf (1, 2, 3)
%!error<normcdf: invalid value for alpha.> [p, plo, pup] = ...
%! normcdf (1, 2, 3, [1, 0; 0, 1], 0)
%!error<normcdf: invalid value for alpha.> [p, plo, pup] = ...
%! normcdf (1, 2, 3, [1, 0; 0, 1], 1.22)
%!error<normcdf: invalid value for alpha.> [p, plo, pup] = ...
%! normcdf (1, 2, 3, [1, 0; 0, 1], "alpha", "upper")
%!error<normcdf: X, MU, and SIGMA must not be complex.> normcdf (i, 2, 2)
%!error<normcdf: X, MU, and SIGMA must not be complex.> normcdf (2, i, 2)
%!error<normcdf: X, MU, and SIGMA must not be complex.> normcdf (2, 2, i)
%!error<normcdf: bad covariance matrix.> ...
%! [p, plo, pup] =normcdf (1, 2, 3, [1, 0; 0, -inf], 0.04)
