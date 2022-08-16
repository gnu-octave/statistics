## Copyright (C) 2022 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Function File} [@var{f}, @var{x}] = ecdf (@var{y})
## @deftypefnx {Function File} [@var{f}, @var{x}, @var{flo}, @var{fup}] = ecdf (@var{y})
## @deftypefnx {Function File} ecdf (@dots{})
## @deftypefnx {Function File} ecdf (@var{ax}, @dots{})
## @deftypefnx {Function File} [@dots{}] = ecdf (@var{y}, @var{name}, @var{value}, @dots{})
## @deftypefnx {Function File} [@dots{}] = ecdf (@var{ax}, @var{y}, @var{name}, @var{value}, @dots{})
##
## Empirical (Kaplan-Meier) cumulative distribution function.
##
## @code{[@var{f}, @var{x}] = ecdf (@var{y})} calculates the Kaplan-Meier
## estimate of the cumulative distribution function (cdf), also known as the
## empirical cdf.  @var{y} is a vector of data values.  @var{f} is a vector of
## values of the empirical cdf evaluated at @var{x}.
##
## @code{[@var{f}, @var{x}, @var{flo}, @var{fup}] = ecdf (@var{y})} also returns
## lower and upper confidence bounds for the cdf.  These bounds are calculated
## using Greenwood's formula, and are not simultaneous confidence bounds.
##
## @code{ecdf (@dots{})} without output arguments produces a plot of the
## empirical cdf.
##
## @code{ecdf (@var{ax}, @dots{})} plots into existing axes @var{ax}.
##
## @code{[@dots{}] = ecdf (@var{y}, @var{name}, @var{value}, @dots{})} specifies
## additional parameter name/value pairs chosen from the following:
##
## @multitable @columnfractions 0.20 0.8
## @headitem Name @tab Value
## @item "censoring" @tab A boolean vector of the same size as Y that is 1 for
## observations that are right-censored and 0 for observations that are observed
## exactly.  Default is all observations observed exactly.
##
## @item "frequency" @tab A vector of the same size as Y containing non-negative
## integer counts.  The jth element of this vector gives the number of times the
## jth element of Y was observed.  Default is 1 observation per Y element.
##
## @item "alpha" @tab A value @var{alpha} between 0 and 1 specifying the
## significance level.  Default is 0.05 for 5% significance.
##
## @item "function" @tab The type of function returned as the F output argument,
## chosen from "cdf" (the default), "survivor", or "cumulative hazard".
##
## @item "bounds" @tab Either "on" to include bounds or "off" (the default) to
## omit them.  Used only for plotting.
## @end multitable
##
## Type @code{demo ecdf} to see examples of usage.
##
## @seealso{cdfplot, ecdfhist}
## @end deftypefn

function [Fout, x, Flo, Fup] = ecdf (y, varargin)

  ## Check for valid input arguments
  narginchk (1, Inf);
  ## Parse input arguments
  if (nargin == 1)
    ax = [];
    if (! isvector (y) || ! isreal (y))
      error ("ecdf: Y must be a vector of real numbers.");
    endif
    ##x = varargin{1};
  else
    ## ax = y;
    ## Check that ax is a valid axis handle
    try
      isstruct (get (y));
      ax = y;
      y = varargin{1};
      varargin{1} = [];
    catch
      ##error ("ecdf: invalid handle %f.", ax);
      ax = [];
    end_try_catch
    #y = varargin{1};
    #varargin{1} = [];
  endif
  ## Make y a column vector
  x = y(:);
  ## Add defaults
  cens = zeros (size (x));
  freq = ones (size (x));
  alpha = 0.05;
  fname = "cdf";
  bound = "off";
  ## Check for remaining varargins and parse extra parameters
  if (length (varargin) > 0 && mod (numel (varargin), 2) == 0)
    [~, prop] = parseparams (varargin);
    while (!isempty (prop))
      switch (lower (prop{1}))
        case "censoring"
          cens = prop{2};
          ## Check for valid input
          if (! isequal (size (x), size (cens)))
            error ("ecdf: censoring data mismatch Y vector.");
          endif
          ## Make double in case censoring data is logical
          if (islogical (cens))
            cens = double (cens);
          endif
        case "frequency"
          freq = prop{2};
          ## Check for valid input
          if (! isequal (size (x), size (freq)))
            error ("ecdf: frequency data mismatch Y vector.");
          endif
          ## Make double in case frequency data is logical
          if (islogical (freq))
            freq = double (freq);
          endif
        case "alpha"
          alpha = prop{2};
          ## Check for valid alpha value
          if (numel (alpha) != 1 || ! isnumeric (alpha) || alpha <= 0 || alpha >= 1)
            error ("ecdf: alpha must be a numeric scalar in the range (0,1).");
          endif
        case "function"
          fname = prop{2};
          ## Check for valid function name option
          if (sum (strcmpi (fname, {"cdf", "survivor", "cumulative hazard"})) < 1)
            error ("ecdf: wrong function name.");
          endif
        case "bounds"
          bound = prop{2};
          ## Check for valid bounds option
          if (! (strcmpi (bound, "off") || strcmpi (bound, "on")))
            error ("ecdf: wrong bounds.");
          endif
        otherwise
          error ("ecdf: unknown option %s", prop{1});
      endswitch
      prop = prop(3:end);
    endwhile
  elseif nargin > 2
    error ("ecdf: optional parameters must be in name/value pairs.");
  endif
  ## Remove NaNs from data
  rn = ! isnan (x) & ! isnan (cens) & freq > 0;
  x = x(rn);
  if (length (x) == 0)
    error ("ecdf: not enought data.");
  endif
  cens = cens(rn);
  freq = freq(rn);
  ## Sort data in ascending order
  [x, sr] = sort (x);
  cens = cens(sr);
  freq = freq(sr);
  ## Keep class for data (single | double)
  if (isa (x, "single"))
    freq = single (freq);
  endif
  ## Calculate cumulative frequencies
  tcfreq = cumsum (freq);
  ocfreq = cumsum (freq .* ! cens);
  x_diff = (diff (x) == 0);
  if (any (x_diff))
    x(x_diff) = [];
    tcfreq(x_diff) = [];
    ocfreq(x_diff) = [];
  endif
  max_count = tcfreq(end);
  ## Get Deaths and Number at Risk for each unique X
  Death = [ocfreq(1); diff(ocfreq)];
  NRisk = max_count - [0; tcfreq(1:end-1)];
  ## Remove no Death observations
  x = x(Death > 0);
  NRisk = NRisk(Death > 0);
  Death = Death(Death > 0);
  ## Estimate function
  switch (fname)
    case "cdf"
      S = cumprod (1 - Death ./ NRisk);
      Fun_x = 1 - S;
      Fzero = 0;
      fdisp = "F(x)";
    case "survivor"
      S = cumprod (1 - Death ./ NRisk);
      Fun_x = S;
      Fzero = 1;
      fdisp = "S(x)";
    case "cumulative hazard"
      Fun_x = cumsum (Death ./ NRisk);
      Fzero = 0;
      fdisp = "H(x)";
  endswitch
  ## Check for remaining non-censored data and add a starting value
  if (! isempty (Death))
    x = [min(y); x];
    F = [Fzero; Fun_x];
  else
    warning("ecdf: No Death in data");
    F = Fun_x;
  endif
  ## Calculate lower and upper confidence bounds if requested
  if (nargout > 2 || (nargout == 0 && strcmpi (bound, "on")))
    switch (fname)
      case {"cdf", "survivor"}
        se = NaN (size (Death));
        if (! isempty (Death))
          if (NRisk(end) == Death(end))
            t = 1:length (NRisk) - 1;
          else
            t = 1:length (NRisk);
          endif
          se(t) = S(t) .* sqrt (cumsum (Death(t) ./ ...
                  (NRisk(t) .* (NRisk(t) - Death(t)))));
        endif
      case "cumulative hazard"
        se = sqrt (cumsum (Death ./ (NRisk .* NRisk)));
    endswitch
    ## Calculate confidence limits
    if (! isempty (se))
      z_a = - norminv (alpha / 2);
      h_w = z_a * se;
      Flo = max (0, Fun_x - h_w);
      Flo(isnan (h_w)) = NaN;
      switch (fname)
        case {"cdf", "survivor"}
          Fup = min (1, Fun_x + h_w);
          Fup(isnan (h_w)) = NaN;
        case "cumulative hazard"
          Fup = Fun_x + h_w;
      endswitch
      Flo = [NaN; Flo];
      Fup = [NaN; Fup];
    else
      Flo = [];
      Fup = [];
    endif
  else
    Flo = [];
    Fup = [];
  endif
  ## Plot stairs if no output is requested
  if (nargout == 0)
    if (isempty (ax))
      ax = newplot();
    end
    h = stairs(ax, x , [F, Flo, Fup]);
    xlabel (ax, "x");
    ylabel (ax, fdisp);
    title ("ecdf");
  else
    Fout = F;
  endif
endfunction

%!demo
%! y = exprnd (10, 50, 1);    ## random failure times are exponential(10)
%! d = exprnd (20, 50, 1);    ## drop-out times are exponential(20)
%! t = min (y, d);            ## we observe the minimum of these times
%! censored = (y > d);        ## we also observe whether the subject failed
%!
%! ## Calculate and plot the empirical cdf and confidence bounds
%! [f, x, flo, fup] = ecdf (t, "censoring", censored);
%! stairs (x, f);
%! hold on;
%! stairs (x, flo, "r:"); stairs (x, fup, "r:");
%!
%! ## Superimpose a plot of the known true cdf
%! xx = 0:.1:max (t); yy = 1 - exp (-xx / 10); plot (xx, yy, "g-");
%! hold off;

%!demo
%! R = wblrnd (100, 2, 100, 1);
%! ecdf (R, "Function", "survivor", "Alpha", 0.01, "Bounds", "on");
%! hold on
%! x = 1:1:250;
%! wblsurv = 1 - cdf ("weibull", x, 100, 2);
%! plot (x, wblsurv, "g-", "LineWidth", 2)
%! legend ("Empirical survivor function", "Lower confidence bound", ...
%!         "Upper confidence bound", "Weibull survivor function", ...
%!         "Location", "northeast");
%! hold off

## Test input
%!error ecdf ();
%!error ecdf (randi (15,2));
%!error ecdf ([3,2,4,3+2i,5]);
%!error kstest ([2,3,4,5,6],"tail");
%!error kstest ([2,3,4,5,6],"tail", "whatever");
%!error kstest ([2,3,4,5,6],"function", "");
%!error kstest ([2,3,4,5,6],"badoption", 0.51);
%!error kstest ([2,3,4,5,6],"tail", 0);
%!error kstest ([2,3,4,5,6],"alpha", 0);
%!error kstest ([2,3,4,5,6],"alpha", NaN);
%!error kstest ([NaN,NaN,NaN,NaN,NaN],"tail", "unequal");
%!error kstest ([2,3,4,5,6],"alpha", 0.05, "CDF", [2,3,4;1,3,4;1,2,1]);
## Test output against MATLAB results
%!test
%! x = [2, 3, 4, 3, 5, 4, 6, 5, 8, 3, 7, 8, 9, 0];
%! [F, x, Flo, Fup] = ecdf (x);
%! F_out = [0; 0.0714; 0.1429; 0.3571; 0.5; 0.6429; 0.7143; 0.7857; 0.9286; 1];
%! assert (F, F_out, ones (10,1) * 1e-4);
%! x_out = [0 0 2 3 4 5 6 7 8 9]';
%! assert (x, x_out);
%! Flo_out = [NaN, 0, 0, 0.1061, 0.2381, 0.3919, 0.4776, 0.5708, 0.7937, NaN]';
%! assert (Flo, Flo_out, ones (10,1) * 1e-4);
%! Fup_out = [NaN, 0.2063, 0.3262, 0.6081, 0.7619, 0.8939, 0.9509, 1, 1, NaN]';
%! assert (Fup, Fup_out, ones (10,1) * 1e-4);
%!shared visibility_setting
%! visibility_setting = get (0, "DefaultFigureVisible");
%!test
%! x = [2, 3, 4, 3, 5, 4, 6, 5, 8, 3, 7, 8, 9, 0];
%! ecdf (x);
%! set (0, "DefaultFigureVisible", visibility_setting);
