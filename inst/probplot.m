## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {} probplot (@var{y})
## @deftypefnx {statistics} {} probplot (@var{dist}, @var{y})
## @deftypefnx {statistics} {} probplot (@var{dist}, @var{y}, @var{cens})
## @deftypefnx {statistics} {} probplot (@var{dist}, @var{y}, @var{cens}, @var{freq})
## @deftypefnx {statistics} {} probplot (@var{ax}, @dots{})
## @deftypefnx {statistics} {} probplot (@dots{}, @qcode{"noref"})
## @deftypefnx {statistics} {@var{h} =} probplot (@dots{})
##
## Produce a probability plot of the data in @var{y} against the distribution
## @var{dist}.
##
## On a probability plot the ordered data are drawn against a nonlinear
## probability axis chosen so that a sample from the reference distribution
## @var{dist} falls approximately along a straight line.  Systematic departures
## from the reference line indicate departures from the distribution.
##
## @var{dist} is one of @qcode{"normal"} (the default when @var{dist} is
## omitted), @qcode{"lognormal"}, @qcode{"exponential"}, @qcode{"extreme value"},
## @qcode{"weibull"}, @qcode{"rayleigh"}, @qcode{"logistic"}, or
## @qcode{"loglogistic"}.  For @qcode{"lognormal"}, @qcode{"weibull"}, and
## @qcode{"loglogistic"} the data axis is logarithmic.
##
## @var{y} is a numeric vector, or a matrix in which case each column is plotted
## as a separate sample.  @code{NaN} values are ignored.
##
## @var{cens} is a logical vector the same size as @var{y} that is true for
## right-censored observations; censored points are not plotted and the plotting
## positions of the remaining points follow the Kaplan-Meier estimate.  @var{freq}
## is a vector of nonnegative integer frequencies (counts) the same size as
## @var{y}.  Pass @code{[]} to omit either one.
##
## @code{probplot (@var{ax}, @dots{})} plots into the axes @var{ax} instead of
## the current axes.  The trailing option @qcode{"noref"} suppresses the
## reference line.
##
## @code{@var{h} = probplot (@dots{})} returns a column vector of handles to the
## plotted line objects (the data markers, followed by the reference line unless
## @qcode{"noref"} was given).
##
## The reference line is a robust fit through the first and third quartiles of the
## data on the transformed scale.  For censored data the quartiles are taken from
## the Kaplan-Meier plotting positions; when heavy censoring prevents the data
## from reaching a quartile the position is linearly extrapolated, which may
## deviate slightly from @sc{matlab}.
##
## @seealso{normplot, wblplot, qqplot, cdfplot, ecdf}
## @end deftypefn

function h = probplot (varargin)

  if (nargin < 1)
    print_usage ();
  endif

  args = varargin;

  ## Optional leading axes handle.
  ax = [];
  if (isscalar (args{1}) && ishghandle (args{1}) ...
      && strcmp (get (args{1}, "type"), "axes"))
    ax = args{1};
    args(1) = [];
  endif
  if (isempty (args))
    print_usage ();
  endif

  ## Distribution name (default "normal") followed by the data.
  if (ischar (args{1}))
    dist = lower (args{1});
    args(1) = [];
  else
    dist = "normal";
  endif
  if (isempty (args))
    error ("probplot: missing data vector Y.");
  endif
  y = args{1};
  args(1) = [];

  ## Remaining arguments: numeric CENS / FREQ and the trailing "noref" flag.
  noref = false;
  numargs = {};
  for k = 1:numel (args)
    a = args{k};
    if (ischar (a))
      if (strcmpi (a, "noref"))
        noref = true;
      else
        error ("probplot: unknown option '%s'.", a);
      endif
    else
      numargs{end+1} = a;
    endif
  endfor
  if (numel (numargs) > 2)
    error ("probplot: too many input arguments.");
  endif
  cens = [];
  freq = [];
  if (numel (numargs) >= 1)
    cens = numargs{1};
  endif
  if (numel (numargs) >= 2)
    freq = numargs{2};
  endif

  ## Look up the distribution transform and axis scale.
  [invfun, logx, prettyname] = dist_transform (dist);

  ## Validate the data.
  if (! (isnumeric (y) && isreal (y)))
    error ("probplot: Y must be real numeric.");
  endif
  if (isrow (y))
    y = y(:);
  endif
  if (ndims (y) > 2)
    error ("probplot: Y must be a vector or a 2-D matrix.");
  endif
  ncol = columns (y);
  if (ncol > 1 && (! isempty (cens) || ! isempty (freq)))
    error ("probplot: CENS and FREQ are only supported for a vector Y.");
  endif

  if (isempty (ax))
    ax = newplot ();
  endif

  ## Probability grid for the y-axis ticks (as in a normal probability plot).
  pgrid = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.25, 0.5, ...
           0.75, 0.9, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999];

  hmark = [];
  href = [];
  hold (ax, "on");
  for j = 1:ncol
    yj = y(:,j);
    cj = cens;
    fj = freq;
    [xu, pp] = plot_positions (yj, cj, fj);
    if (isempty (xu))
      continue;
    endif
    zz = invfun (pp);
    hmark(end+1) = line (ax, xu, zz, "linestyle", "none", "marker", "x", ...
                                     "color", "b");
    if (! noref)
      href(end+1) = ref_line (ax, xu, pp, invfun, logx);
    endif
  endfor
  hold (ax, "off");

  ## Data (x) axis scale.
  if (logx)
    set (ax, "xscale", "log");
  else
    set (ax, "xscale", "linear");
  endif

  ## Probability (y) axis: tick at invfun(pgrid), labelled with the probability.
  ok = isfinite (invfun (pgrid));
  yt = invfun (pgrid(ok));
  ylabels = arrayfun (@(p) num2str (p), pgrid(ok), "UniformOutput", false);
  set (ax, "ytick", yt, "yticklabel", ylabels);

  title (ax, sprintf ("Probability plot for %s distribution", prettyname));
  xlabel (ax, "Data");
  ylabel (ax, "Probability");
  grid (ax, "on");
  box (ax, "off");

  if (nargout > 0)
    h = [hmark(:); href(:)];
  endif

endfunction

## Standardized inverse cdf (transform), log-x flag, and pretty name per dist.
function [invfun, logx, prettyname] = dist_transform (dist)
  switch (dist)
    case "normal"
      invfun = @(p) norminv (p); logx = false; prettyname = "normal";
    case "lognormal"
      invfun = @(p) norminv (p); logx = true;  prettyname = "lognormal";
    case "exponential"
      invfun = @(p) -log (1 - p); logx = false; prettyname = "exponential";
    case {"extreme value", "ev"}
      invfun = @(p) log (-log (1 - p)); logx = false;
      prettyname = "extreme value";
    case {"weibull", "wbl"}
      invfun = @(p) log (-log (1 - p)); logx = true; prettyname = "weibull";
    case {"rayleigh", "rayl"}
      invfun = @(p) sqrt (-2 .* log (1 - p)); logx = false;
      prettyname = "rayleigh";
    case "logistic"
      invfun = @(p) log (p ./ (1 - p)); logx = false; prettyname = "logistic";
    case "loglogistic"
      invfun = @(p) log (p ./ (1 - p)); logx = true;
      prettyname = "loglogistic";
    otherwise
      error ("probplot: unrecognized distribution '%s'.", dist);
  endswitch
endfunction

## Sorted uncensored data and Kaplan-Meier survival-midpoint plotting positions.
function [xu, pp] = plot_positions (y, cens, freq)

  keep = ! isnan (y);
  y = y(keep);
  if (isempty (cens))
    cens = false (size (y));
  else
    cens = logical (cens(:));
    cens = cens(keep);
  endif
  if (isempty (freq))
    freq = ones (size (y));
  else
    freq = freq(:);
    freq = freq(keep);
  endif

  ## Expand by integer frequency.
  y = repelem (y(:), freq);
  cens = repelem (cens(:), freq);

  [ys, ord] = sort (y);
  cs = cens(ord);

  ## Kaplan-Meier survival, sampled to the midpoint of each uncensored jump.
  n = numel (ys);
  S = 1;
  nrisk = n;
  xu = [];
  pp = [];
  for i = 1:n
    if (! cs(i))
      Safter = S * (nrisk - 1) / nrisk;
      xu(end+1) = ys(i);
      pp(end+1) = 1 - (S + Safter) / 2;
      S = Safter;
    endif
    nrisk -= 1;
  endfor
  xu = xu(:);
  pp = pp(:);

endfunction

## Robust quartile reference line drawn across the data range.
function hl = ref_line (ax, xu, pp, invfun, logx)

  ## Quartiles of the data from the plotting positions (transformed scale).
  if (logx)
    u = log (xu);
  else
    u = xu;
  endif
  ## interp1 needs strictly increasing, unique sample sites.
  [ppu, iu] = unique (pp);
  q = interp1 (ppu, u(iu), [0.25, 0.75], "linear", "extrap");
  z = invfun ([0.25, 0.75]);
  slope = (z(2) - z(1)) / (q(2) - q(1));
  intercept = z(1) - slope * q(1);

  ## Sample the straight line across the data range on the transformed scale.
  uu = linspace (min (u), max (u), 100)';
  zz = intercept + slope * uu;
  if (logx)
    xx = exp (uu);
  else
    xx = uu;
  endif
  hl = line (ax, xx, zz, "linestyle", "--", "marker", "none", "color", "r");

endfunction

%!demo
%! probplot ([1:20]);

%!demo
%! probplot ("weibull", [1 2 3 4 5 6 7 8 9 10 15 20]);

## shared data
%!shared y, c, f
%! y = [2.1 3.4 1.8 5.2 2.9 4.1 3.3 2.7 6.0 3.8 4.5 2.2];
%! c = logical ([0 0 0 1 0 0 0 0 1 0 1 0]);
%! f = [1 2 1 1 3 1 1 2 1 1 1 2];

## normal: marker positions and reference line, verified against MATLAB
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   h = probplot (y);
%!   assert_equal (numel (h), 2);
%!   xd = get (h(1), "xdata");
%!   yd = get (h(1), "ydata");
%!   assert_equal (xd(:)', sort (y), 1e-12);
%!   ymatlab = [-1.7317 -1.1503 -0.8122 -0.5485 -0.3186 -0.1046 ...
%!               0.1046 0.3186 0.5485 0.8122 1.1503 1.7317];
%!   assert_equal (yd(:)', ymatlab, 1e-4);
%!   rl = get (h(2), "ydata");
%!   rx = get (h(2), "xdata");
%!   ## reference line slope matches MATLAB robust quartile fit
%!   slope = (rl(end) - rl(1)) / (rx(end) - rx(1));
%!   assert_equal (slope, 0.72923, 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## weibull uses a logarithmic data axis and the SEV transform
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   ax = axes (hf);
%!   h = probplot (ax, "weibull", y);
%!   assert_equal (get (ax, "xscale"), "log");
%!   yd = get (h(1), "ydata");
%!   ymatlab = [-3.1568 -2.0134 -1.4541 -1.0647 -0.7550 -0.4892 ...
%!              -0.2483 -0.0194 0.2088 0.4502 0.7321 1.1563];
%!   assert_equal (yd(:)', ymatlab, 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## remaining distributions: transform and data-axis scale vs MATLAB
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   dists = {"exponential", "lognormal", "rayleigh", "logistic", ...
%!            "loglogistic", "extreme value"};
%!   gold = {[0.0426 0.1335 0.2336 0.3448 0.4700 0.6131 0.7802 0.9808 ...
%!            1.2321 1.5686 2.0794 3.1781], ...
%!           [-1.7317 -1.1503 -0.8122 -0.5485 -0.3186 -0.1046 0.1046 ...
%!             0.3186 0.5485 0.8122 1.1503 1.7317], ...
%!           [0.2918 0.5168 0.6835 0.8305 0.9695 1.1073 1.2491 1.4006 ...
%!            1.5698 1.7712 2.0393 2.5211], ...
%!           [-3.1355 -1.9459 -1.3350 -0.8873 -0.5108 -0.1671 0.1671 ...
%!             0.5108 0.8873 1.3350 1.9459 3.1355], ...
%!           [-3.1355 -1.9459 -1.3350 -0.8873 -0.5108 -0.1671 0.1671 ...
%!             0.5108 0.8873 1.3350 1.9459 3.1355], ...
%!           [-3.1568 -2.0134 -1.4541 -1.0647 -0.7550 -0.4892 -0.2483 ...
%!            -0.0194 0.2088 0.4502 0.7321 1.1563]};
%!   logx = {false, true, false, false, true, false};
%!   for k = 1:numel (dists)
%!     ax = axes (hf);
%!     h = probplot (ax, dists{k}, y);
%!     assert_equal (get (h(1), "ydata")(:)', gold{k}, 1e-4);
%!     if (logx{k})
%!       assert_equal (get (ax, "xscale"), "log");
%!     else
%!       assert_equal (get (ax, "xscale"), "linear");
%!     endif
%!     delete (ax);
%!   endfor
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## noref suppresses the reference line
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   h = probplot (y, "noref");
%!   assert_equal (numel (h), 1);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## censoring: only uncensored points plotted, at Kaplan-Meier positions
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   h = probplot ("normal", y, c);
%!   xd = get (h(1), "xdata");
%!   yd = get (h(1), "ydata");
%!   assert_equal (xd(:)', [1.8 2.1 2.2 2.7 2.9 3.3 3.4 3.8 4.1], 1e-12);
%!   ymatlab = [-1.7317 -1.1503 -0.8122 -0.5485 -0.3186 -0.1046 ...
%!               0.1046 0.3186 0.5485];
%!   assert_equal (yd(:)', ymatlab, 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## interspersed censoring: KM redistributes the plotting positions
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   cm = logical ([0 0 0 0 0 1 1 1 0 0 0 0]);
%!   h = probplot ("normal", y, cm);
%!   yd = get (h(1), "ydata");
%!   ymatlab = [-1.7317 -1.1503 -0.8122 -0.5334 -0.2574 0.0196 ...
%!               0.3462 0.7764 1.4544];
%!   assert_equal (yd(:)', ymatlab, 1e-4);
%!   rl = get (h(2), "ydata");
%!   rx = get (h(2), "xdata");
%!   slope = (rl(end) - rl(1)) / (rx(end) - rx(1));
%!   assert_equal (slope, 0.53519, 1e-4);
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## frequency: each observation expanded into repeated markers
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   h = probplot ("normal", y, [], f);
%!   xd = get (h(1), "xdata");
%!   assert_equal (numel (xd), sum (f));
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect

## Test input validation
%!error <Invalid call to probplot> probplot ()
%!error <unrecognized distribution 'foo'.> probplot ("foo", [1 2 3 4])
%!error <unknown option 'bar'.> probplot ([1 2 3 4], "bar")
