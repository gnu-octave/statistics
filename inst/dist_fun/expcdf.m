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
## @deftypefn  {statistics} {@var{p} =} expcdf (@var{x})
## @deftypefnx {statistics} {@var{p} =} expcdf (@var{x}, @var{mu})
## @deftypefnx {statistics} {@var{p} =} expcdf (@dots{}, @qcode{"upper"})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} expcdf (@var{x}, @var{mu}, @var{pcov})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} expcdf (@var{x}, @var{mu}, @var{pcov}, @var{alpha})
## @deftypefnx {statistics} {[@var{p}, @var{plo}, @var{pup}] =} expcdf (@dots{}, @qcode{"upper"})
##
## Exponential cumulative distribution function (CDF).
##
## For each element of @var{x}, compute the cumulative distribution function
## (CDF) of the exponential distribution with mean parameter @var{mu}.  The size
## of @var{p} is the common size of @var{x} and @var{mu}.  A scalar input
## functions as a constant matrix of the same size as the other inputs.
##
## Default value for @var{mu} = 1.
##
## A common alternative parameterization of the exponential distribution is to
## use the parameter @math{λ} defined as the mean number of events in an
## interval as opposed to the parameter @math{μ}, which is the mean wait time
## for an event to occur. @math{λ} and @math{μ} are reciprocals,
## i.e. @math{μ = 1 / λ}.
##
## When called with three output arguments, i.e. @code{[@var{p}, @var{plo},
## @var{pup}]}, it computes the confidence bounds for @var{p} when the input
## parameter @var{mu} is an estimate.  In such case, @var{pcov}, a scalar value
## with the variance of the estimated parameter @var{mu}, is necessary.
## Optionally, @var{alpha}, which has a default value of 0.05, specifies the
## @qcode{100 * (1 - @var{alpha})} percent confidence bounds.  @var{plo} and
## @var{pup} are arrays of the same size as @var{p} containing the lower and
## upper confidence bounds.
##
## @code{[@dots{}] = expcdf (@dots{}, "upper")} computes the upper tail
## probability of the exponential distribution.
##
## Further information about the exponential distribution can be found at
## @url{https://en.wikipedia.org/wiki/Exponential_distribution}
##
## @seealso{expinv, exppdf, exprnd, expfit, explike, expstat}
## @end deftypefn

function [varargout] = expcdf (x, varargin)

  ## Check for valid number of input arguments
  if (nargin < 1 || nargin > 5)
    error ("expcdf: invalid number of input arguments.");
  endif

  ## Check for "upper" flag
  if (nargin > 1 && strcmpi (varargin{end}, "upper"))
    uflag = true;
    varargin(end) = [];
  elseif (nargin > 1 && ischar (varargin{end}) && ...
          ! strcmpi (varargin{end}, "upper"))
    error ("expcdf: invalid argument for upper tail.");
  else
    uflag = false;
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
      error ("expcdf: invalid size of variance, PCOV must be a scalar.");
    endif
    if (pcov < 0)
      error ("expcdf: variance, PCOV, cannot be negative.");
    endif
  else
    ## Check that cov matrix is provided if 3 output arguments are requested
    if (nargout > 1)
      error ("expcdf: variance, PCOV, is required for confidence bounds.");
    endif
    pcov = [];
  endif
  if (numel (varargin) > 2)
    alpha = varargin{3};
    ## Check for valid alpha value
    if (! isnumeric (alpha) || numel (alpha) != 1 || alpha <= 0 || alpha >= 1)
      error ("expcdf: invalid value for alpha.");
    endif
  else
    alpha = 0.05;
  endif

  ## Check for common size of X and MU
  if (! isscalar (x) || ! isscalar (mu))
    [err, x, mu] = common_size (x, mu);
    if (err > 0)
      error ("expcdf: X and MU must be of common size or scalars.");
    endif
  endif

  ## Check for X and MU being reals
  if (iscomplex (x) || iscomplex (mu))
    error ("expcdf: X and MU must not be complex.");
  endif

  ## Check for appropriate class
  if (isa (x, "single") || isa (mu, "single"));
    is_class = "single";
  else
    is_class = "double";
  endif

  ## Return NaNs for out of range parameters.
  mu(mu <= 0) = NaN;

  ## Compute P value for exponential cdf
  z = x ./ mu;

  ## Force 0 for negative X
  z(z < 0) = 0;

  ## Check uflag
  if (uflag)
    p = exp (-z);
  else
    p = -expm1 (-z);
  endif

  ## Prepare output
  varargout{1} = cast (p, is_class);
  if (nargout > 1)
    plo = NaN (size (z), is_class);
    pup = NaN (size (z), is_class);
  endif

  ## Compute confidence bounds (if requested)
  if (nargout >= 2)

    ## Convert to log scale
    log_z = log (z);
    norm_z = -probit (alpha / 2);
    halfwidth = norm_z * sqrt (pcov ./ (mu .^ 2));
    zlo = log_z - halfwidth;
    zup = log_z + halfwidth;

    ## Convert to original scale
    if (uflag)
      plo = exp (-exp (zup));
      pup = exp (-exp (zlo));
    else
       plo = - expm1 (-exp (zlo));
       pup = - expm1 (-exp (zup));
    endif

    ## Prepare output
    varargout{2} = plo;
    varargout{3} = pup;
  endif

endfunction

%!demo
%! ## Plot various CDFs from the exponential distribution
%! x = 0:0.01:5;
%! p1 = expcdf (x, 2/3);
%! p2 = expcdf (x, 1.0);
%! p3 = expcdf (x, 2.0);
%! plot (x, p1, "-b", x, p2, "-g", x, p3, "-r")
%! grid on
%! legend ({"μ = 2/3", "μ = 1", "μ = 2"}, "location", "southeast")
%! title ("Exponential CDF")
%! xlabel ("values in x")
%! ylabel ("probability")

## Test output
%!shared x, p
%! x = [-1 0 0.5 1 Inf];
%! p = [0, 1 - exp(-x(2:end)/2)];
%!assert (expcdf (x, 2*ones (1,5)), p)
%!assert (expcdf (x, 2), p)
%!assert (expcdf (x, 2*[1 0 NaN 1 1]), [0 NaN NaN p(4:5)])
## Test class of input preserved
%!assert (expcdf ([x, NaN], 2), [p, NaN])
%!assert (expcdf (single ([x, NaN]), 2), single ([p, NaN]))
%!assert (expcdf ([x, NaN], single (2)), single ([p, NaN]))

## Test values against MATLAB output
%!test
%! [p, plo, pup] = expcdf (1, 2, 3);
%! assert (p, 0.39346934028737, 1e-14);
%! assert (plo, 0.08751307220484, 1e-14);
%! assert (pup, 0.93476821257933, 1e-14);
%!test
%! [p, plo, pup] = expcdf (1, 2, 2, 0.1);
%! assert (p, 0.39346934028737, 1e-14);
%! assert (plo, 0.14466318041675, 1e-14);
%! assert (pup, 0.79808291849140, 1e-14);
%!test
%! [p, plo, pup] = expcdf (1, 2, 2, 0.1, "upper");
%! assert (p, 0.60653065971263, 1e-14);
%! assert (plo, 0.20191708150860, 1e-14);
%! assert (pup, 0.85533681958325, 1e-14);

## Test input validation
%!error<expcdf: invalid number of input arguments.> expcdf ()
%!error<expcdf: invalid number of input arguments.> expcdf (1, 2 ,3 ,4 ,5, 6)
%!error<expcdf: invalid argument for upper tail.> expcdf (1, 2, 3, 4, "uper")
%!error<expcdf: X and MU must be of common size or scalars.> ...
%! expcdf (ones (3), ones (2))
%!error<expcdf: invalid size of variance, PCOV must be a scalar.> ...
%! expcdf (2, 3, [1, 2])
%!error<expcdf: variance, PCOV, is required for confidence bounds.> ...
%! [p, plo, pup] = expcdf (1, 2)
%!error<expcdf: invalid value for alpha.> [p, plo, pup] = ...
%! expcdf (1, 2, 3, 0)
%!error<expcdf: invalid value for alpha.> [p, plo, pup] = ...
%! expcdf (1, 2, 3, 1.22)
%!error<expcdf: invalid value for alpha.> [p, plo, pup] = ...
%! expcdf (1, 2, 3, "alpha", "upper")
%!error<expcdf: X and MU must not be complex.> expcdf (i, 2)
%!error<expcdf: X and MU must not be complex.> expcdf (2, i)
%!error<expcdf: variance, PCOV, cannot be negative.> ...
%! [p, plo, pup] = expcdf (1, 2, -1, 0.04)
