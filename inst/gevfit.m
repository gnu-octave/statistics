## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2012-2021 Nir Krakauer <nkrakauer@ccny.cuny.edu>
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
## @deftypefn {Function File} [@var{paramhat}, @var{paramci}] = gevfit (@var{data})
## @deftypefnx {Function File} [@var{paramhat}, @var{paramci}] = gevfit (@var{data}, @var{paramguess})
## @deftypefnx {Function File} [@var{paramhat}, @var{paramci}] = gevfit (@var{data}, @var{alpha})
## @deftypefnx {Function File} [@var{paramhat}, @var{paramci}] = gevfit (@dots{}, @var{options})
##
## Find the maximum likelihood estimator @var{paramhat} of the generalized
## extreme value (GEV) distribution to fit @var{data}.
##
## @subheading Arguments
##
## @itemize @bullet
## @item
## @var{data} is the vector of given values.
## @item
## @var{paramguess} is an initial guess for the maximum likelihood parameter
## vector. If not given, this defaults to @var{k_0}=0 and @var{sigma}, @var{mu}
## determined by matching the data mean and standard deviation to their expected
## values.
## @item
## @var{alpha} returns 100(1-ALPHA) percent confidence intervals for the
## parameter estimates.  Pass in [] for ALPHA to use the default values.
## @item
## @var{options} a structure that specifies control parameters for the iterative
## algorithm used to compute ML estimates.  The structure must contain the
## following fields:
##
## 'Display' = "off"; 'MaxFunEvals' = 400; 'MaxIter' = 200;
## 'TolFun' = 1e-6; 'TolX' = 1e-6.
##
## If not provided, the aforementioned values are used by default.
## @end itemize
##
## @subheading Return values
##
## @itemize @bullet
## @item
## @var{paramhat} is the 3-parameter maximum-likelihood parameter vector
## @code{[@var{k_0}, @var{sigma}, @var{mu}]}, where @var{k_0} is the shape
## parameter of the GEV distribution, @var{sigma} is the scale parameter of the
## GEV distribution, and @var{mu} is the location parameter of the GEV
## distribution.
## @item
## @var{paramci} has the approximate 95% confidence intervals of the parameter
## values based on the Fisher information matrix at the maximum-likelihood
## position.
##
## @end itemize
##
## When K < 0, the GEV is the type III extreme value distribution.  When K > 0,
## the GEV distribution is the type II, or Frechet, extreme value distribution.
## If W has a Weibull distribution as computed by the WBLFIT function, then -W
## has a type III extreme value distribution and 1/W has a type II extreme value
## distribution.  In the limit as K approaches 0, the GEV is the mirror image of
## the type I extreme value distribution as computed by the EVFIT function.
##
## The mean of the GEV distribution is not finite when K >= 1, and the variance
## is not finite when PSI >= 1/2.  The GEV distribution is defined
## for K*(X-MU)/SIGMA > -1.
##
## @subheading Examples
##
## @example
## @group
## data = 1:50;
## [pfit, pci] = gevfit (data);
## p1 = gevcdf (data, pfit(1), pfit(2), pfit(3));
## plot (data, p1)
## @end group
## @end example
## @seealso{gevcdf, gevinv, gevlike, gevpdf, gevrnd, gevstat}
## @end deftypefn

function [paramhat, paramci] = gevfit (data, varargin)

  ## Check arguments
  if (nargin < 1)
    print_usage;
  endif
  ## Check data is vector
  if (! isvector (data))
    error ("gevfit: DATA must be a vector");
  endif
  ## Get data type and convert to double for computation
  is_type = class (data);
  if (strcmpi (is_type, "single"))
    data = double (data);
  endif
  ## Check that data is not constant and does not contain NaNs
  sample_size = length (data);
  if (sample_size == 0 || any (isnan (data)))
    paramhat = NaN (1,3, is_type);
    paramci = NaN (2,3, is_type);
    warning ("gevfit: data contain NaNs");
    return
  elseif (numel (unique (data)) == 1)
    paramhat = cast ([0, 0, unique(data)], is_type);
    if (length (data) == 1)
      paramci = cast ([-Inf, 0, -Inf; Inf, Inf, Inf], is_type);
    else
      paramci = [paramhat; paramhat];
    endif
    warning ("gevfit: data is a constant vector");
    return
  endif
  ## Check second input argument
  ## If it is a scalar, then it is ALPHA value,
  ## otherwise it is an initial parameter vector
  paramguess = [];
  alpha = [];
  if (nargin > 1 && isequal (size (varargin{1}), [1, 3]))
    paramguess = varargin{1};
  elseif (nargin > 1 && isscalar (varargin{1}) && varargin{1} > 0 ...
                     && varargin{1} < 1)
    alpha = varargin{1};
  endif
  ## Compute initial parameters if not parsed as an input argument
  if (isempty (paramguess))
    F = (0.5:1:(sample_size - 0.5))' ./ sample_size;
    k_0 = fminsearch (@(k) 1 - corr (data, gevinv (F, k, 1, 0)), 0);
    paramguess = [k_0, polyfit(gevinv(F,k_0,1,0),data,1)];
    #paramguess = [k_0, tmp(1), tmp(2)];
    ## Check if data support initial parameters or fall back to unbounded evfit
    if (k_0 < 0 && (max (data) > - paramguess(2) / k_0 + paramguess(3)) || ...
        k_0 > 0 && (min (data) < - paramguess(2) / k_0 + paramguess(3)))
      paramguess = [evfit(data), 0];
      paramguess = flip (paramguess);
    endif
  endif
  if (isempty (alpha))
    alpha = 0.05;
  endif
  ## Get options structure or add defaults
  if (nargin == 3)
    options = varargin{3};
  else
    options.Display = "off";
    options.MaxFunEvals = 400;
    options.MaxIter = 200;
    options.TolFun = 1e-6;
    options.TolX = 1e-6;
  endif
  ## Minimize the negative log-likelihood according to initial parameters
  paramguess(2) = log (paramguess(2));
  fhandle = @(paramguess) nll (paramguess, data);
  [paramhat, ~, exitflag, output] = fminsearch (fhandle, paramguess, options);
  paramhat(2) = exp (paramhat(2));
  ## Display errors and warnings if any
  if (exitflag == 0)
    if (output.funcCount >= output.iterations)
      warning ("gevfit: maximum number of evaluations reached");
    else
      warning ("gevfit: reached iteration limit");
    endif
  elseif (exitflag == -1)
    error ("gevfit: No solution");
  endif
  ## Return a row vector for Matlab compatibility
  paramhat = paramhat(:)';
  ## Check for second output argument
  if (nargout > 1)
  	[~, ~, ACOV] = gevlike (paramhat, data);
  	param_se = sqrt (diag (ACOV))';
    if (any (iscomplex (param_se)))
      warning (["gevfit: Fisher information matrix not positive definite;", ...
                " parameter optimization likely did not converge"]);
      paramci = NaN (2, 3, is_type);
    else
      p_vals = [alpha/2; 1-alpha/2];
      k_ci = norminv (p_vals, paramhat(1), param_se(1));
      s_ci = exp (norminv (p_vals, log (paramhat(2)), param_se(2) ./ paramhat(2)));
      m_ci = norminv (p_vals, paramhat(3), param_se(3));
      paramci = [k_ci, s_ci, m_ci];
    endif
  endif
endfunction

## Negative log-likelihood for the GEV (log(sigma) parameterization)
function out = nll (parms, data)
  k_0 = parms(1);
  log_sigma = parms(2);
  sigma = exp (log_sigma);
  mu = parms(3);
  n = numel (data);
  z = (data - mu) ./ sigma;
  if abs(k_0) > eps
      u = 1 + k_0.*z;
      if min(u) > 0
          lnu = log1p (k_0 .* z);
          out = n * log_sigma + sum (exp ((-1 / k_0) * lnu)) + ...
                                (1 + 1 / k_0) * sum (lnu);
      else
          out = Inf;
      endif
  else
      out = n * log_sigma + sum (exp (-z) + z);
  endif
endfunction

%!demo
%! data = 1:50;
%! [pfit, pci] = gevfit (data);
%! p1 = gevcdf (data, pfit(1), pfit(2), pfit(3));
%! plot (data, p1);

%!test
%! data = 1:50;
%! [pfit, pci] = gevfit (data);
%! pfit_out = [-0.4407, 15.1923, 21.5309];
%! pci_out = [-0.7532, 11.5878, 16.5686; -0.1282, 19.9183, 26.4926];
%! assert (pfit, pfit_out, 1e-3);
%! assert (pci, pci_out, 1e-3);
%!error [pfit, pci] = gevfit (ones (2,5));

%!test
%! data = 1:2:50;
%! [pfit, pci] = gevfit (data);
%! pfit_out = [-0.4434, 15.2024, 21.0532];
%! pci_out = [-0.8904, 10.3439, 14.0168; 0.0035, 22.3429, 28.0896];
%! assert (pfit, pfit_out, 1e-3);
%! assert (pci, pci_out, 1e-3);
