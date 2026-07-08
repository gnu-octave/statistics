## Copyright (C) 2026 Aman Behera <aman.behera.systesms@gmail.com>
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
## @deftypefn  {statistics} {@var{h} =} __anova_plot_diagnostics__ (@var{residuals}, @var{fitted}, @var{leverage}, @var{cooksd}, @var{dfe})
## @deftypefnx {statistics} {@var{h} =} __anova_plot_diagnostics__ (@dots{}, @var{name}, @var{value})
##
## Shared private helper that produces the four standard regression
## diagnostic panels (Normal Q-Q, Spread-Location, Residual-Leverage,
## Cook's distance) used by both @code{anovan} and the @code{anova}
## classdef.  Mirrors the existing implementation in @code{anovan.m}
## (this skeleton is finalized in Week 5).
##
## @end deftypefn

function h = __anova_plot_diagnostics__ (residuals, fitted, leverage, ...
                                         cooksd, dfe, varargin)

  if (nargin < 5)
    error ("__anova_plot_diagnostics__: too few input arguments.");
  endif

  if (isempty (residuals) || isempty (fitted) || isempty (leverage) ...
      || isempty (cooksd))
    error ("__anova_plot_diagnostics__: diagnostic inputs must be non-empty.");
  endif

  residuals = residuals(:);
  fitted    = fitted(:);
  leverage  = leverage(:);
  cooksd    = cooksd(:);
  n = numel (residuals);
  if (numel (fitted) != n || numel (leverage) != n || numel (cooksd) != n)
    error ("__anova_plot_diagnostics__: diagnostic inputs must match.");
  endif

  fig_name = "Diagnostic Plots: Model Residuals";
  visible = "on";
  if (mod (numel (varargin), 2) != 0)
    error ("__anova_plot_diagnostics__: name-value pairs must come in pairs.");
  endif
  for k = 1:2:numel (varargin)
    switch (lower (varargin{k}))
      case "figurename"
        fig_name = varargin{k + 1};
      case "visible"
        visible = varargin{k + 1};
      otherwise
        error ("__anova_plot_diagnostics__: unknown option '%s'.", varargin{k});
    endswitch
  endfor

  mse = sum (residuals .^ 2) / max (dfe, 1);
  t = residuals ./ sqrt (mse * max (1 - leverage, eps));
  [~, DI] = sort (cooksd, "descend");
  nk = min (4, n);

  h = figure ("Name", fig_name, "Visible", visible);

  ## Normal Q-Q
  subplot (2, 2, 1);
  x = ((1:n)' - 0.5) / n;
  [ts, I] = sort (t);
  q = norminv (x);
  plot (q, ts, "ok", "markersize", 3);
  box off; grid on;
  xlabel ("Theoretical quantiles");
  ylabel ("Studentized residuals");
  title ("Normal Q-Q Plot");
  arrayfun (@(i) text (q(I == DI(i)), t(DI(i)), ...
                       sprintf ("  %u", DI(i))), 1:nk);
  iqr = [0.25; 0.75];
  yl = quantile (t, iqr, 1, 6);
  xl = norminv (iqr);
  slope = diff (yl) / diff (xl);
  int = yl(1) - slope * xl(1);
  ax1_xlim = get (gca, "XLim");
  hold on;
  plot (ax1_xlim, slope * ax1_xlim + int, "k-");
  hold off;
  set (gca, "Xlim", ax1_xlim);

  ## Spread-Location
  subplot (2, 2, 2);
  plot (fitted, sqrt (abs (t)), "ko", "markersize", 3);
  box off;
  xlabel ("Fitted values");
  ylabel ("sqrt ( | Studentized residuals | )");
  title ("Spread-Location Plot");
  ax2_xlim = get (gca, "XLim");
  hold on;
  plot (ax2_xlim, ones (1, 2) * sqrt (2), "k:");
  plot (ax2_xlim, ones (1, 2) * sqrt (3), "k-.");
  plot (ax2_xlim, ones (1, 2) * sqrt (4), "k--");
  hold off;
  arrayfun (@(i) text (fitted(DI(i)), sqrt (abs (t(DI(i)))), ...
                       sprintf ("  %u", DI(i))), 1:nk);
  xlim (ax2_xlim);

  ## Residual-Leverage
  subplot (2, 2, 3);
  plot (leverage, t, "ko", "markersize", 3);
  box off;
  xlabel ("Leverage");
  ylabel ("Studentized residuals");
  title ("Residual-Leverage Plot");
  ax3_xlim = get (gca, "XLim");
  ax3_ylim = get (gca, "YLim");
  hold on;
  plot (ax3_xlim, zeros (1, 2), "k-");
  hold off;
  arrayfun (@(i) text (leverage(DI(i)), t(DI(i)), ...
                       sprintf ("  %u", DI(i))), 1:nk);
  set (gca, "ygrid", "on");
  xlim (ax3_xlim);
  ylim (ax3_ylim);

  ## Cook's distance
  subplot (2, 2, 4);
  stem (cooksd, "ko", "markersize", 3);
  box off;
  xlabel ("Obs. number");
  ylabel ("Cook's distance");
  title ("Cook's Distance Stem Plot");
  xlim ([0, n]);
  ax4_xlim = get (gca, "XLim");
  ax4_ylim = get (gca, "YLim");
  hold on;
  plot (ax4_xlim, ones (1, 2) * 4 / max (dfe, eps), "k:");
  plot (ax4_xlim, ones (1, 2) * 0.5, "k-.");
  plot (ax4_xlim, ones (1, 2), "k--");
  hold off;
  arrayfun (@(i) text (DI(i), cooksd(DI(i)), ...
                       sprintf ("  %u", DI(i))), 1:nk);
  xlim (ax4_xlim);
  ylim (ax4_ylim);

  set (findall (gcf, "-property", "FontSize"), "FontSize", 7);

endfunction

%!demo
%! y  = [10; 12; 11; 14; 16; 15; 9; 8; 10];
%! g  = [1;1;1;2;2;2;3;3;3];
%! [~, ~, stats] = anovan (y, {g}, "display", "off");
%! h = full (stats.X) * pinv (full (stats.X' * stats.X)) * stats.X';
%! lev = diag (h);
%! mse = stats.mse;
%! __anova_plot_diagnostics__ (stats.resid, full (stats.X) * stats.coeffs(:,1), ...
%!                             lev, stats.CooksD, stats.dfe);

%!error <__anova_plot_diagnostics__: too few input arguments.> ...
%!  __anova_plot_diagnostics__ ([], [], [], [])
%!error <diagnostic inputs must match> ...
%!  __anova_plot_diagnostics__ ([1;2], [1], [0;0], [0;0], 1)
