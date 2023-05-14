## Copyright (C) 2012 Rik Wehbring
## Copyright (C) 1995-2016 Kurt Hornik
## Copyright (C) 2023 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {@var{x} =} expinv (@var{p})
## @deftypefnx {statistics} {@var{x} =} expinv (@var{p}, @var{mu})
## @deftypefnx {statistics} {[@var{x}, @var{xlo}, @var{xup}] =} expinv (@var{p}, @var{mu}, @var{pcov})
## @deftypefnx {statistics} {[@var{x}, @var{xlo}, @var{xup}] =} expinv (@var{p}, @var{mu}, @var{pcov}, @var{alpha})
##
## Inverse of the exponential cumulative distribution function (iCDF).
##
## For each element of @var{p}, compute the quantile (the inverse of the CDF) of
## the exponential distribution with mean @var{mu}.  The size of @var{x} is the
## common size of @var{p} and @var{mu}.  A scalar input functions as a constant
## matrix of the same size as the other inputs.
##
## Default value is @var{mu} = 1.
##
## A common alternative parameterization of the exponential distribution is to
## use the parameter @math{λ} defined as the mean number of events in an
## interval as opposed to the parameter @math{μ}, which is the mean wait time
## for an event to occur. @math{λ} and @math{μ} are reciprocals,
## i.e. @math{μ = 1 / λ}.
##
## When called with three output arguments, i.e. @qcode{[@var{x}, @var{xlo},
## @var{xup}]}, @code{expinv} computes the confidence bounds for @var{x} when
## the input parameter @var{mu} is an estimate.  In such case, @var{pcov}, a
## scalar value with the variance of the estimated parameter @var{mu}, is
## necessary.  Optionally, @var{alpha}, which has a default value of 0.05,
## specifies the @qcode{100 * (1 - @var{alpha})} percent confidence bounds.
## @var{xlo} and @var{xup} are arrays of the same size as @var{x} containing the
## lower and upper confidence bounds.
##
## Further information about the exponential distribution can be found at
## @url{https://en.wikipedia.org/wiki/Exponential_distribution}
##
## @seealso{expcdf, exppdf, exprnd, expfit, explike, expstat}
## @end deftypefn

function [varargout] = expinv (p, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 4)
    error ("expinv: invalid number of input arguments.");
  endif

  ## Get extra arguments (if they exist) or add defaults
  if (numel (varargin) > 0)
    mu = varargin{1};
  else
    mu = 1;
  endif
  if (numel (varargin) > 1)
    pcov = varargin{2};
    ## Check for variance being a scalar
    if (! isscalar (pcov))
      error ("expinv: invalid size of variance, PCOV must be a scalar.");
    endif
    if (pcov < 0)
      error ("expinv: variance, PCOV, cannot be negative.");
    endif
  else
    ## Check that cov matrix is provided if 3 output arguments are requested
    if (nargout > 1)
      error ("expinv: variance, PCOV, is required for confidence bounds.");
    endif
    pcov = [];
  endif
  if (numel (varargin) > 2)
    alpha = varargin{3};
    ## Check for valid alpha value
    if (! isnumeric (alpha) || numel (alpha) != 1 || alpha <= 0 || alpha >= 1)
      error ("expinv: invalid value for alpha.");
    endif
  else
    alpha = 0.05;
  endif

  ## Check for common size of P and MU
  if (! isscalar (p) || ! isscalar (mu))
    [retval, p, mu] = common_size (p, mu);
    if (retval > 0)
      error ("expinv: P and MU must be of common size or scalars.");
    endif
  endif

  ## Check for P and MU being reals
  if (iscomplex (p) || iscomplex (mu))
    error ("expinv: P and MU must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (p, "single") || isa (mu, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Create output matrix
  if (isa (p, "single") || isa (mu, "single"))
    x = NaN (size (p), "single");
  else
    x = NaN (size (p));
  endif

  ## Handle edge cases
  k = (p == 1) & (mu > 0);
  x(k) = Inf;

  ## Handle valid cases
  k = (p >= 0) & (p < 1) & (mu > 0);
  if (isscalar (mu))
    x(k) = - mu * log (1 - p(k));
  else
    x(k) = - mu(k) .* log (1 - p(k));
  endif

  ## Prepare output
  varargout{1} = cast (x, is_class);
  if (nargout > 1)
    xlo = NaN (size (z), is_class);
    xup = NaN (size (z), is_class);
  endif

  ## Compute confidence bounds (if requested)
  if (nargout >= 2)

    ## Convert to log scale
    log_x = log (x);
    z = -probit (alpha / 2);
    halfwidth = z * sqrt (pcov ./ (mu.^2));

    ## Convert to original scale
    xlo = exp (log_x - halfwidth);
    xup = exp (log_x + halfwidth);

    ## Prepare output
    varargout{2} = plo;
    varargout{3} = pup;
  endif

endfunction

%!demo
%! ## Plot various iCDFs from the exponential distribution
%! p = 0.001:0.001:0.999;
%! x1 = expinv (p, 2/3);
%! x2 = expinv (p, 1.0);
%! x3 = expinv (p, 2.0);
%! plot (p, x1, "-b", p, x2, "-g", p, x3, "-r")
%! grid on
%! ylim ([0, 5])
%! legend ({"μ = 2/3", "μ = 1", "μ = 2"}, "location", "northwest")
%! title ("Exponential iCDF")
%! xlabel ("probability")
%! ylabel ("values in x")

## Test output
%!shared p
%! p = [-1 0 0.3934693402873666 1 2];
%!assert (expinv (p, 2*ones (1,5)), [NaN 0 1 Inf NaN], eps)
%!assert (expinv (p, 2), [NaN 0 1 Inf NaN], eps)
%!assert (expinv (p, 2*[1 0 NaN 1 1]), [NaN NaN NaN Inf NaN], eps)
%!assert (expinv ([p(1:2) NaN p(4:5)], 2), [NaN 0 NaN Inf NaN], eps)

## Test class of input preserved
%!assert (expinv ([p, NaN], 2), [NaN 0 1 Inf NaN NaN], eps)
%!assert (expinv (single ([p, NaN]), 2), single ([NaN 0 1 Inf NaN NaN]), eps)
%!assert (expinv ([p, NaN], single (2)), single ([NaN 0 1 Inf NaN NaN]), eps)

## Test input validation
%!error<expinv: invalid number of input arguments.> expinv ()
%!error<expinv: invalid number of input arguments.> expinv (1, 2 ,3 ,4 ,5)
%!error<expinv: P and MU must be of common size or scalars.> ...
%! expinv (ones (3), ones (2))
%!error<expinv: invalid size of variance, PCOV must be a scalar.> ...
%! expinv (2, 3, [1, 2])
%!error<expinv: variance, PCOV, is required for confidence bounds.> ...
%! [x, xlo, xup] = expinv (1, 2)
%!error<expinv: invalid value for alpha.> [x, xlo, xup] = ...
%! expinv (1, 2, 3, 0)
%!error<expinv: invalid value for alpha.> [x, xlo, xup] = ...
%! expinv (1, 2, 3, 1.22)
%!error<expinv: invalid value for alpha.> [x, xlo, xup] = ...
%! expinv (1, 2, 3, [0.05, 0.1])
%!error<expinv: P and MU must not be complex.> expinv (i, 2)
%!error<expinv: P and MU must not be complex.> expinv (2, i)
%!error<expinv: variance, PCOV, cannot be negative.> ...
%! [x, xlo, xup] = expinv (1, 2, -1, 0.04)
