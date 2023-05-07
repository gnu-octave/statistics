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
## You should have received lambda copy of the GNU General Public License
## along with this program; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {@var{p} =} wblcdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} wblcdf (@var{x}, @var{lambda})
## @deftypefnx {statistics} {@var{p} =} wblcdf (@var{x}, @var{lambda}, @var{k})
## @deftypefnx {statistics} {@var{p} =} wblcdf (@dots{}, @qcode{"upper"})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} wblcdf (@var{x}, @var{lambda}, @var{k}, @var{pcov})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} wblcdf (@var{x}, @var{lambda}, @var{k}, @var{pcov}, @var{alpha})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} wblcdf (@dots{}, @qcode{"upper"})
##
## Weibull cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) at @var{x} of the Weibull distribution with scale parameter
## @var{lambda} and shape parameter @var{k}.  The size of @var{p} is the common
## size of @var{x}, @var{lambda} and @var{k}.  A scalar input functions as a
## constant matrix of the same size as the other inputs.
##
## Default values are @var{lambda} = 0, @var{k} = 1.
##
## When called with three output arguments, @code{[@var{p}, @var{plo},
## @var{pup}]} it computes the confidence bounds for @var{p} when the input
## parameters @var{lambda} and @var{k} are estimates.  In such case, @var{pcov},
## a 2-by-2 matrix containing the covariance matrix of the estimated parameters,
## is necessary.  Optionally, @var{alpha} has a default value of 0.05, and
## specifies 100 * (1 - @var{alpha})% confidence bounds. @var{plo} and @var{pup}
## are arrays of the same size as @var{p} containing the lower and upper
## confidence bounds.
##
## @code{[@dots{}] = wblcdf (@dots{}, "upper")} computes the upper tail
## probability of the lognormal distribution.
##
## Further information about the Weibull distribution can be found at
## @url{https://en.wikipedia.org/wiki/Weibull_distribution}
##
## @seealso{wblinv, wblpdf, wblrnd, wblstat, wblplot}
## @end deftypefn

function [varargout] = wblcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 6)
    error ("wblcdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1  && ischar (varargin{end}) && ...
          ! strcmpi (varargin{end}, "upper"))
    error ("wblcdf: invalid argument for upper tail.");
  elseif (nargin > 1 && isempty (varargin{end}))
    uflag = false;
    varargin(end) = [];
  else
    uflag = false;
  endif

  ## Get extra arguments (if they exist) or add defaults
  if (numel (varargin) > 0)
    lambda = varargin{1};
  else
    lambda = 1;
  endif
  if (numel (varargin) > 1)
    k = varargin{2};
  else
    k = 1;
  endif
  if (numel (varargin) > 2)
    pcov = varargin{3};
    ## Check for valid covariance matrix 2x2
    if (! isequal (size (pcov), [2, 2]))
      error ("wblcdf: invalid size of covariance matrix.");
    endif
  else
    ## Check that cov matrix is provided if 3 output arguments are requested
    if (nargout > 1)
      error ("wblcdf: covariance matrix is required for confidence bounds.");
    endif
    pcov = [];
  endif
  if (numel (varargin) > 3)
    alpha = varargin{4};
    ## Check for valid alpha value
    if (! isnumeric (alpha) || numel (alpha) !=1 || alpha <= 0 || alpha >= 1)
      error ("wblcdf: invalid value for alpha.");
    endif
  else
    alpha = 0.05;
  endif

  ## Check for common size of X, LAMBDA, and K
  if (! isscalar (x) || ! isscalar (lambda) || ! isscalar (k))
    [err, x, lambda, k] = common_size (x, lambda, k);
    if (err > 0)
      error ("wblcdf: X, LAMBDA, and K must be of common size or scalars.");
    endif
  endif

  ## Check for X, LAMBDA, and K being reals
  if (iscomplex (x) || iscomplex (lambda) || iscomplex (k))
    error ("wblcdf: X, LAMBDA, and K must not be complex.");
  endif

  ## Return NaN for out of range parameters.
  lambda(lambda <= 0) = NaN;
  k(k <= 0) = NaN;

  ## Force 0 for negative data
  x(x < 0) = 0;

  ## Compute z
  z = (x ./ lambda) .^ k;
  if (uflag)
    p = exp(-z);
  else
    p = -expm1(-z);
  endif

  ## Compute confidence bounds (if requested)
  if (nargout >= 2)
    ## Work on log scale
    log_z = log (z);
    d_lambda = 1 ./ lambda;
    d_k = -1 ./ (k .^ 2);
    log_zvar = (pcov(1,1) .* d_lambda .^ 2 + ...
               2 * pcov(1,2) .* d_lambda .* d_k .* log_z + ...
               pcov(2,2) .* (d_k .* log_z) .^ 2) .* (k .^ 2);
    if (any(log_zvar < 0))
      error ("wblcdf: bad covariance matrix.");
    endif
    normz = -norminv (alpha / 2);
    halfwidth = normz * sqrt (log_zvar);
    zlo = log_z - halfwidth;
    zup = log_z + halfwidth;
    ## Convert back from log scale
    if uflag == true
      plo = exp (-exp (zup));
      pup = exp (-exp (zlo));
    else
      plo = -expm1 (-exp (zlo));
      pup = -expm1 (-exp (zup));
    endif
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (lambda, "single") || isa (k, "single"));
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
%! ## Plot various CDFs from the Weibull distribution
%! x = 0:0.001:2.5;
%! p1 = wblcdf (x, 1, 0.5);
%! p2 = wblcdf (x, 1, 1);
%! p3 = wblcdf (x, 1, 1.5);
%! p4 = wblcdf (x, 1, 5);
%! plot (x, p1, "-b", x, p2, "-r", x, p3, "-m", x, p4, "-g")
%! grid on
%! legend ({"λ = 1, k = 0.5", "λ = 1, k = 1", ...
%!          "λ = 1, k = 1.5", "λ = 1, k = 5"}, "location", "southeast")
%! title ("Weibull CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, y
%! x = [-1 0 0.5 1 Inf];
%! y = [0, 1-exp(-x(2:4)), 1];
%!assert (wblcdf (x, ones (1,5), ones (1,5)), y)
%!assert (wblcdf (x, ones (1,5), ones (1,5), "upper"), 1 - y)
%!assert (wblcdf (x, "upper"), 1 - y)
%!assert (wblcdf (x, 1, ones (1,5)), y)
%!assert (wblcdf (x, ones (1,5), 1), y)
%!assert (wblcdf (x, [0 1 NaN Inf 1], 1), [NaN 0 NaN 0 1])
%!assert (wblcdf (x, [0 1 NaN Inf 1], 1, "upper"), 1 - [NaN 0 NaN 0 1])
%!assert (wblcdf (x, 1, [0 1 NaN Inf 1]), [NaN 0 NaN y(4:5)])
%!assert (wblcdf (x, 1, [0 1 NaN Inf 1], "upper"), 1 - [NaN 0 NaN y(4:5)])
%!assert (wblcdf ([x(1:2) NaN x(4:5)], 1, 1), [y(1:2) NaN y(4:5)])
%!assert (wblcdf ([x(1:2) NaN x(4:5)], 1, 1, "upper"), 1 - [y(1:2) NaN y(4:5)])

## Test class of input preserved
%!assert (wblcdf ([x, NaN], 1, 1), [y, NaN])
%!assert (wblcdf (single ([x, NaN]), 1, 1), single ([y, NaN]))
%!assert (wblcdf ([x, NaN], single (1), 1), single ([y, NaN]))
%!assert (wblcdf ([x, NaN], 1, single (1)), single ([y, NaN]))

## Test input validation
%!error<wblcdf: invalid number of input arguments.> wblcdf ()
%!error<wblcdf: invalid number of input arguments.> wblcdf (1,2,3,4,5,6,7)
%!error<wblcdf: invalid argument for upper tail.> wblcdf (1, 2, 3, 4, "uper")
%!error<wblcdf: X, LAMBDA, and K must be of common size or scalars.> ...
%! wblcdf (ones (3), ones (2), ones (2))
%!error<wblcdf: invalid size of covariance matrix.> wblcdf (2, 3, 4, [1, 2])
%!error<wblcdf: covariance matrix is required for confidence bounds.> ...
%! [p, plo, pup] = wblcdf (1, 2, 3)
%!error<wblcdf: invalid value for alpha.> [p, plo, pup] = ...
%! wblcdf (1, 2, 3, [1, 0; 0, 1], 0)
%!error<wblcdf: invalid value for alpha.> [p, plo, pup] = ...
%! wblcdf (1, 2, 3, [1, 0; 0, 1], 1.22)
%!error<wblcdf: invalid value for alpha.> [p, plo, pup] = ...
%! wblcdf (1, 2, 3, [1, 0; 0, 1], "alpha", "upper")
%!error<wblcdf: X, LAMBDA, and K must not be complex.> wblcdf (i, 2, 2)
%!error<wblcdf: X, LAMBDA, and K must not be complex.> wblcdf (2, i, 2)
%!error<wblcdf: X, LAMBDA, and K must not be complex.> wblcdf (2, 2, i)
%!error<wblcdf: bad covariance matrix.> ...
%! [p, plo, pup] =wblcdf (1, 2, 3, [1, 0; 0, -inf], 0.04)
