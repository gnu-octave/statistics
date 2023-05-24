## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
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
## @deftypefn  {statistics} {@var{p} =} logncdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} logncdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} logncdf (@var{x}, @var{mu}, @var{sigma})
## @deftypefnx {statistics} {@var{p} =} logncdf (@dots{}, @qcode{"upper"})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} logncdf (@var{x}, @var{mu}, @var{sigma}, @var{pcov})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} logncdf (@var{x}, @var{mu}, @var{sigma}, @var{pcov}, @var{alpha})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} logncdf (@dots{}, @qcode{"upper"})
##
## Log-normal cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the log-normal distribution with mean @var{mu} and standard
## deviation @var{sigma} corresponding to the associated normal distribution.
## The size of @var{p} is the common size of @var{x}, @var{mu} and @var{sigma}.
## A scalar input functions as a constant matrix of the same size as the other
## inputs.
##
## If a random variable follows this distribution, its logarithm is normally
## distributed with mean @var{mu} and standard deviation @var{sigma}.
##
## Default parameter values are @qcode{@var{mu} = 0} and
## @qcode{@var{sigma} = 1}.  Both parameters must be reals and
## @qcode{@var{sigma} > 0}.  For @qcode{@var{sigma} <= 0}, @qcode{NaN} is
## returned.
##
## When called with three output arguments, i.e. @qcode{[@var{p}, @var{plo},
## @var{pup}]}, @code{logncdf} computes the confidence bounds for @var{p} when
## the input parameters @var{mu} and @var{sigma} are estimates.  In such case,
## @var{pcov}, a @math{2x2} matrix containing the covariance matrix of the
## estimated parameters, is necessary.  Optionally, @var{alpha}, which has a
## default value of 0.05, specifies the @qcode{100 * (1 - @var{alpha})} percent
## confidence bounds.  @var{plo} and @var{pup} are arrays of the same size as
## @var{p} containing the lower and upper confidence bounds.
##
## @code{[@dots{}] = logncdf (@dots{}, "upper")} computes the upper tail
## probability of the log-normal distribution with parameters @var{mu} and
## @var{sigma}, at the values in @var{x}.
##
## Further information about the log-normal distribution can be found at
## @url{https://en.wikipedia.org/wiki/Log-normal_distribution}
##
## @seealso{logninv, lognpdf, lognrnd, lognfit, lognlike, lognstat}
## @end deftypefn

function [varargout] = logncdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 6)
    error ("logncdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1  && ischar (varargin{end}) && ...
          ! strcmpi (varargin{end}, "upper"))
    error ("logncdf: invalid argument for upper tail.");
  elseif (nargin > 2 && isempty (varargin{end}))
    uflag = false;
    varargin(end) = [];
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
      error ("logncdf: invalid size of covariance matrix.");
    endif
  else
    ## Check that cov matrix is provided if 3 output arguments are requested
    if (nargout > 1)
      error ("logncdf: covariance matrix is required for confidence bounds.");
    endif
    pcov = [];
  endif
  if (numel (varargin) > 3)
    alpha = varargin{4};
    ## Check for valid alpha value
    if (! isnumeric (alpha) || numel (alpha) !=1 || alpha <= 0 || alpha >= 1)
      error ("logncdf: invalid value for alpha.");
    endif
  else
    alpha = 0.05;
  endif

  ## Check for common size of X, MU, and SIGMA
  if (! isscalar (x) || ! isscalar (mu) || ! isscalar (sigma))
    [err, x, mu, sigma] = common_size (x, mu, sigma);
    if (err > 0)
      error ("logncdf: X, MU, and SIGMA must be of common size or scalars.");
    endif
  endif

  ## Check for X, MU, and SIGMA being reals
  if (iscomplex (x) || iscomplex (mu) || iscomplex (sigma))
    error ("logncdf: X, MU, and SIGMA must not be complex.");
  endif

  ## Return NaN for out of range parameters.
  sigma(sigma <= 0) = NaN;

  ## Negative data would create complex values, which erfc cannot handle.
  x(x < 0) = 0;

  ## Compute lognormal cdf
  z = (log (x) - mu) ./ sigma;
  if (uflag)
    z = -z;
  endif
  p = 0.5 * erfc (-z ./ sqrt(2));

  ## Compute confidence bounds (if requested)
  if (nargout >= 2)
    zvar = (pcov(1,1) + 2 * pcov(1,2) * z + pcov(2,2) * z .^ 2) ./ (sigma .^ 2);
    if (any (zvar(:) < 0))
      error ("logncdf: bad covariance matrix.");
    end
    normz = -norminv (alpha / 2);
    halfwidth = normz * sqrt (zvar);
    zlo = z - halfwidth;
    zup = z + halfwidth;

    plo = 0.5 * erfc (-zlo ./ sqrt (2));
    pup = 0.5 * erfc (-zup ./ sqrt (2));
  endif

  ## Check for class type
  if (isa (x, "single") || isa (mu, "single") || isa (sigma, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Prepare output
  varargout{1} = cast (p, is_class);
  if (nargout > 1)
    varargout{2} = cast (plo, is_class);
    varargout{3} = cast (pup, is_class);
  endif

endfunction

%!demo
%! ## Plot various CDFs from the log-normal distribution
%! x = 0:0.01:3;
%! p1 = logncdf (x, 0, 1);
%! p2 = logncdf (x, 0, 0.5);
%! p3 = logncdf (x, 0, 0.25);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"μ = 0, σ = 1", "μ = 0, σ = 0.5", "μ = 0, σ = 0.25"}, ...
%!         "location", "southeast")
%! title ("Log-normal CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1, 0, 1, e, Inf];
%! y = [0, 0, 0.5, 1/2+1/2*erf(1/2), 1];
%!assert (logncdf (x, zeros (1,5), sqrt(2)*ones (1,5)), y, eps)
%!assert (logncdf (x, zeros (1,5), sqrt(2)*ones (1,5), []), y, eps)
%!assert (logncdf (x, 0, sqrt(2)*ones (1,5)), y, eps)
%!assert (logncdf (x, zeros (1,5), sqrt(2)), y, eps)
%!assert (logncdf (x, [0 1 NaN 0 1], sqrt(2)), [0 0 NaN y(4:5)], eps)
%!assert (logncdf (x, 0, sqrt(2)*[0 NaN Inf 1 1]), [NaN NaN y(3:5)], eps)
%!assert (logncdf ([x(1:3) NaN x(5)], 0, sqrt(2)), [y(1:3) NaN y(5)], eps)

## Test class of input preserved
%!assert (logncdf ([x, NaN], 0, sqrt(2)), [y, NaN], eps)
%!assert (logncdf (single ([x, NaN]), 0, sqrt(2)), single ([y, NaN]), eps ("single"))
%!assert (logncdf ([x, NaN], single (0), sqrt(2)), single ([y, NaN]), eps ("single"))
%!assert (logncdf ([x, NaN], 0, single (sqrt(2))), single ([y, NaN]), eps ("single"))

## Test input validation
%!error<logncdf: invalid number of input arguments.> logncdf ()
%!error<logncdf: invalid number of input arguments.> logncdf (1,2,3,4,5,6,7)
%!error<logncdf: invalid argument for upper tail.> logncdf (1, 2, 3, 4, "uper")
%!error<logncdf: X, MU, and SIGMA must be of common size or scalars.> ...
%! logncdf (ones (3), ones (2), ones (2))
%!error<logncdf: invalid size of covariance matrix.> logncdf (2, 3, 4, [1, 2])
%!error<logncdf: covariance matrix is required for confidence bounds.> ...
%! [p, plo, pup] = logncdf (1, 2, 3)
%!error<logncdf: invalid value for alpha.> [p, plo, pup] = ...
%! logncdf (1, 2, 3, [1, 0; 0, 1], 0)
%!error<logncdf: invalid value for alpha.> [p, plo, pup] = ...
%! logncdf (1, 2, 3, [1, 0; 0, 1], 1.22)
%!error<logncdf: invalid value for alpha.> [p, plo, pup] = ...
%! logncdf (1, 2, 3, [1, 0; 0, 1], "alpha", "upper")
%!error<logncdf: X, MU, and SIGMA must not be complex.> logncdf (i, 2, 2)
%!error<logncdf: X, MU, and SIGMA must not be complex.> logncdf (2, i, 2)
%!error<logncdf: X, MU, and SIGMA must not be complex.> logncdf (2, 2, i)
%!error<logncdf: bad covariance matrix.> ...
%! [p, plo, pup] =logncdf (1, 2, 3, [1, 0; 0, -inf], 0.04)
