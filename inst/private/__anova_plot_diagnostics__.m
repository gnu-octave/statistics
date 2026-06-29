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

  fig_name = "Diagnostic Plots: Model Residuals";
  for k = 1:2:numel (varargin)
    if (strcmpi (varargin{k}, "FigureName"))
      fig_name = varargin{k + 1};
    endif
  endfor

  n = numel (residuals);
  mse = sum (residuals .^ 2) / max (dfe, 1);
  t = residuals ./ sqrt (mse * max (1 - leverage, eps));
  [~, DI] = sort (cooksd, "descend");
  nk = min (4, n);

  h = figure ("Name", fig_name);

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

  ## Spread-Location
  subplot (2, 2, 2);
  plot (fitted, sqrt (abs (t)), "ko", "markersize", 3);
  box off;
  xlabel ("Fitted values");
  ylabel ("sqrt ( | Studentized residuals | )");
  title ("Spread-Location Plot");

  ## Residual-Leverage
  subplot (2, 2, 3);
  plot (leverage, t, "ko", "markersize", 3);
  box off;
  xlabel ("Leverage");
  ylabel ("Studentized residuals");
  title ("Residual-Leverage Plot");

  ## Cook's distance
  subplot (2, 2, 4);
  stem (cooksd, "ko", "markersize", 3);
  box off;
  xlabel ("Observation");
  ylabel ("Cook's distance");
  title ("Cook's Distance");
  xlim ([0, n]);

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
