## Copyright (C) 2014 Bj{\"o}rn Vennberg
##
## This file is part of the statistics package for GNU Octave.
##
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn  {statistics} {} wblplot (@var{data}, @dots{})
## @deftypefnx {statistics} {@var{handle} =} wblplot (@var{data}, @dots{})
## @deftypefnx {statistics} {[@var{handle}, @var{param}] =} wblplot (@var{data})
## @deftypefnx {statistics} {[@var{handle}, @var{param}] =} wblplot (@var{data}, @var{censor})
## @deftypefnx {statistics} {[@var{handle}, @var{param}] =} wblplot (@var{data}, @var{censor}, @var{freq})
## @deftypefnx {statistics} {[@var{handle}, @var{param}] =} wblplot (@var{data}, @var{censor}, @var{freq}, @var{confint})
## @deftypefnx {statistics} {[@var{handle}, @var{param}] =} wblplot (@var{data}, @var{censor}, @var{freq}, @var{confint}, @var{fancygrid})
## @deftypefnx {statistics} {[@var{handle}, @var{param}] =} wblplot (@var{data}, @var{censor}, @var{freq}, @var{confint}, @var{fancygrid}, @var{showlegend})
##
## Plot a column vector @var{data} on a Weibull probability plot using rank
## regression.
##
## @var{censor}: optional parameter is a column vector of same size as
## @var{data} with 1 for right censored data and 0 for exact observation.
## Pass [] when no censor data are available.
##
## @var{freq}: optional vector same size as @var{data} with the number of
## occurences for corresponding data.
## Pass [] when no frequency data are available.
##
## @var{confint}: optional confidence limits for plotting upper and lower
## confidence bands using beta binomial confidence bounds.  If a single
## value is given this will be used such as LOW = a and HIGH = 1 - a.
## Pass [] if confidence bounds is not requested.
##
## @var{fancygrid}: optional parameter which if set to anything but 1 will turn
## off the fancy gridlines.
##
## @var{showlegend}: optional parameter that when set to zero(0) turns off the
## legend.
##
## If one output argument is given, a @var{handle} for the data marker and
## plotlines is returned, which can be used for further modification of line and
## marker style.
##
## If a second output argument is specified, a @var{param} vector with scale,
## shape and correlation factor is returned.
##
## @seealso{normplot, wblpdf}
## @end deftypefn

function [handle, param] = wblplot (data, censor = [], freq = [], ...
                                    confint = [], fancygrid = 1, showlegend = 1)
  [mm, nn] = size (data);
  if (mm > 1 && nn > 1)
    error ("wblplot: can only handle a single data vector")
  elseif (mm == 1 && nn > 1)
    data = data(:);
    mm = nn;
  endif
  if (any (data <= 0))
    error ("wblplot: data vector must be positive and non zero")
  endif

  if (isempty (freq))
    freq = ones (mm, 1);
    N = mm;
  else
    [mmf nnf] = size (freq);
    if ((mmf == mm && nnf == 1) || (mmf == 1 && nnf == mm))
      freq = freq(:);
      N = sum (freq);	## Total number of samples
      if (any (freq <= 0))
	      error ("wblplot: frequency vector must be positive non zero integers")
      endif
    else
      error ("wblplot: frequency must be vector of same length as data")
    endif
  endif

  if (isempty (censor))
      censor = zeros(mm,1);
  else
    [mmc, nnc] = size(censor);
    if ((mmc == mm && nnc == 1) || (mmc == 1 && nnc == mm))
      censor = censor(:);
    else
      error ("wblplot: censor must be a vector of same length as data")
    endif
    ## Make sure censored data is sorted corectly so that no censored samples
    ## are processed before failures if they have the same time.
    if (any (censor > 0))
       ind = find (censor > 0);
       ind2 = find (data(1:end-1) == data(2:end));
       if ((! isempty (ind)) && (! isempty (ind2)))
         if (any (ind == ind2))
           tmp = censor(ind2);
           censor(ind2) = censor(ind2 + 1);
           censor(ind2+1) = tmp;
           tmp = freq(ind2);
           freq(ind2) = freq(ind2 + 1);
           freq(ind2 + 1) = tmp;
         endif
       endif
    endif
  endif

  ## Determine the order number
  wbdat = zeros (length (find (censor == 0)), 3);
  Op = 0;
  Oi = 0;
  c = N;
  nf = 0;
  for k = 1 : mm
    if (censor(k, 1) == 0)
      nf = nf + 1;
      wbdat(nf, 1) = data(k, 1);
      for s = 1 : freq(k, 1);
        Oi = Op + ((N + 1) - Op) / (1 + c);
        Op = Oi;
        c = c - 1;
      endfor
      wbdat(nf, 3) = Oi;
    else
      c = c - freq(k, 1);
    endif
  endfor
  ## Compute median rank
  a = wbdat(:, 3) ./ (N - wbdat(:, 3) + 1);
  f = finv(0.5, 2 * (N - wbdat(:, 3) + 1), 2 * wbdat(:, 3));

  wbdat(:, 2) = a ./ (f+a);

  datx = log (wbdat(:,1));
  daty = log (log (1 ./ (1 - wbdat(:,2))));

  ## Rank regression
  poly = polyfit (datx, daty, 1);
  ## Shape factor
  beta_rry = poly(1);
  ## Scale factor
  eta_rry = exp (-(poly(2) / beta_rry));

  ## Determine min-max values of view port
  aa = ceil (log10 (max (wbdat(:,1))));
  bb = log10 (max (wbdat(:,1)));
  if ((aa - bb) < 0.2)
    aa = ceil (log10 (max (wbdat(:,1)))) + 1;
  endif
  xmax = 10 ^ aa;

  if ((log10 (min (wbdat(:,1))) - floor (log10 (min (wbdat(:,1))))) < 0.2)
    xmin = 10 ^ (floor (log10 (min (wbdat(:,1)))) - 1);
  else
    xmin = 10 ^ floor (log10 (min (wbdat(:,1))));
  endif

  if (min (wbdat(:,2)) > 0.20)
    ymin = log (log (1 / (1 - 0.1)));
  elseif (min (wbdat(:,2)) > 0.02)
    ymin = log (log (1 / (1 - 0.01)));
  elseif (min (wbdat(:,2)) > 0.002)
    ymin = log (log (1 / (1 - 0.001)));
  else
    ymin = log (log (1 / (1 - 0.0001)));
  endif

  ymax= log (log (1 / (1 - 0.999)));
  x = [0;0];
  y = [0;0];

  label = char('0.10', '1.00', '10.00', '99.00');
  prob = [0.001 0.01 0.1 0.99];
  tick  = log (log (1 ./ (1 - prob)));
  xbf = [xmin; xmax];
  ybf = polyval (poly, log (xbf));

  newplot();
  x(1, 1) = xmin;
  x(2, 1) = xmax;
  if (fancygrid == 1)
    for k = 1 : 4
      ## Y major grids
      x(1, 1) = xmin;
      x(2, 1) = xmax*10;
      y(1, 1) = log (log (1 / (1 - 10 ^ (-k))));
      y(2, 1) = y(1, 1);
      ymajorgrid(k) = line (x, y, 'LineStyle', '-', 'Marker', 'none', ...
                      'Color', [1 0.75 0.75], 'LineWidth', 0.1);
    endfor
    ## Y Minor grids 2 - 9
    x(1, 1) = xmin;
    x(2, 1) = xmax * 10;
    for m = 1 : 4
      for k = 1 : 8
        y(1, 1) = log (log (1 / (1 - ((k + 1) / (10 ^ m)))));
        y(2, 1) = y(1, 1);
        yminorgrid(k) = line (x, y, 'LineStyle', '-', 'Marker', 'none', ...
                        'Color', [0.75 1 0.75], 'LineWidth', 0.1);
      endfor
    endfor
    ## X-axis grid
    y(1, 1) = ymin;
    y(2, 1) = ymax;
    for m = log10 (xmin) : log10 (xmax)
      x(1, 1) = 10 ^ m;
      x(2, 1) =  x(1, 1);
      y(1, 1) = ymin;
      y(2, 1) = ymax;
      xmajorgrid(k) = line (x, y, 'LineStyle', '-', 'Marker', 'none', ...
                      'Color', [1 0.75 0.75]);
      for k = 1 : 8
        ## X Minor grids - 2 - 9
        x(1, 1) = (k + 1) * (10 ^ m);
        x(2, 1) = (k + 1) * (10 ^ m);
        xminorgrid(k) = line (x, y, 'LineStyle', '-', 'Marker', 'none', ...
                        'Color', [0.75 1 0.75], 'LineWidth', 0.1);
      endfor
    endfor
  endif

  set (gca, 'XScale', 'log');
  set (gca, 'YTick', tick, 'YTickLabel', label);

  xlabel ('Data', 'FontSize', 12);
  ylabel ('Unreliability, F(t)=1-R(t)', 'FontSize', 12);
  title ('Weibull Probability Plot', 'FontSize', 12);
  set (gcf, 'Color', [0.9, 0.9, 0.9]);
  set (gcf, 'name', 'WblPlot');
  hold on

  h = plot (wbdat(:,1), daty, 'o');
  set (h, 'markerfacecolor', [0, 0, 1]);
  set (h, 'markersize', 8);
  h2 = line (xbf, ybf, 'LineStyle', '-', 'Marker', 'none', ...
            'Color', [0.25 0.25 1], 'LineWidth', 1);
  ## If requested plot beta binomial confidens bounds
  if (! isempty (confint))
    cb_high = [];
    cb_low = [];
    if (length (confint) == 1)
      if (confint > 0.5)
        cb_high = confint;
        cb_low = 1 - confint;
      else
        cb_high = 1 - confint;
        cb_low = confint;
      endif
    else
      cb_high = confint(2);
      cb_low = confint(1);
    endif
    conf = zeros (N + 4, 3);
    betainv = 1 / beta_rry;
    N2 = [1:N]';
    N2 = [0.3; 0.7; N2; N2(end) + 0.5; N2(end) + 0.8]; ## Extend the ends a bit
    ypos = medianranks (0.5, N, N2);
    conf(:, 1) = eta_rry * log (1 ./ (1 - ypos)) .^ betainv;
    conf(:, 2) = medianranks (cb_low, N, N2);
    conf(:, 3) = medianranks (cb_high, N, N2);
    confy = log (log (1 ./ (1 - conf(:,2:3))));
    confu = [conf(:,1) confy];

    if (conf(1,1) > xmin)  ## It looks better to extend the lines.
      p1 = polyfit (log (conf(1:2,1)), confy(1:2,1), 1);
      y1 = polyval (p1, log (xmin));
      p2 = polyfit (log (conf(1:2,1)), confy(1:2,2), 1);
      y2 = polyval (p2, log (xmin));
      confu = [xmin y1 y2; confu];
    endif

    if (conf(end,1) < xmax)
      p3 = polyfit (log (conf(end-1:end,1)), confy(end-1:end,1), 1);
      y3 = polyval (p3, log (xmax));
      p4 = polyfit (log (conf(end-1:end,1)), confy(end-1:end,2), 1);
      y4 = polyval (p4, log (xmax));
      confu = [confu; xmax, y3, y4];
    endif
    h3 = plot (confu(:,1), confu(:,2:3), 'LineStyle', '-' ,'Marker', 'none', ...
               'Color', [1 0.25 0.25], 'LineWidth', 1);
  endif
  ## Correlation coefficient
  rsq = corr (datx, daty);

  if (showlegend == 1)
    s1 = sprintf (' RRY\n \\beta=%.3f \n \\eta=%.2f  \n \\rho=%.4f', ...
                  beta_rry, eta_rry, rsq);
    if (! isempty (confint))
      s2 = sprintf ('CB_H=%.2f', cb_high);
      s3 = sprintf ('CB_L=%.2f', cb_low);
      legend ([h; h2; h3], "Data", s1, s2, s3, "location", "northeastoutside");
    else
      legend ([h; h2], "Data", s1, "location", "northeastoutside");
    endif
    legend ("boxoff");
  endif
  axis ([xmin, xmax, ymin, (log (log (1 / (1 - 0.99))))]);
  hold off

  if (nargout >= 2)
    param = [eta_rry, beta_rry, rsq];
    if (! isempty (confint))
      handle = [h; h2; h3];
    else
      handle = [h; h2];
    endif
  endif
  if (nargout == 1)
    if (! isempty (confint))
      handle = [h; h2; h3];
    else
      handle = [h; h2];
    endif
  endif
endfunction


function ret = medianranks (alpha, n, ii)
  a = ii ./ (n - ii + 1);
  f = finv (alpha, 2 * (n - ii + 1), 2 * ii);
  ret = a ./ (f + a);
endfunction

%!demo
%! x = [16 34 53 75 93 120];
%! wblplot (x);

%!demo
%! x = [2 3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61 67]';
%! c = [0 1 0 1 0 1 1 1 0 0 1 0 1 0 1 1 0 1 1]';
%! [h, p] = wblplot (x, c);
%! p

%!demo
%! x = [16, 34, 53, 75, 93, 120, 150, 191, 240 ,339];
%! [h, p] = wblplot (x, [], [], 0.05);
%! p
%! ## Benchmark Reliasoft eta = 146.2545 beta 1.1973 rho = 0.9999

%!demo
%! x = [46 64 83 105 123 150 150];
%! c = [0 0 0 0 0 0 1];
%! f = [1 1 1 1 1 1 4];
%! wblplot (x, c, f, 0.05);

%!demo
%! x = [46 64 83 105 123 150 150];
%! c = [0 0 0 0 0 0 1];
%! f = [1 1 1 1 1 1 4];
%! ## Subtract 30.92 from x to simulate a 3 parameter wbl with gamma = 30.92
%! wblplot (x - 30.92, c, f, 0.05);

## Test plotting
%!test
%! hf = figure ("visible", "off");
%! unwind_protect
%!   x = [16, 34, 53, 75, 93, 120, 150, 191, 240 ,339];
%!   [h, p] = wblplot (x, [], [], 0.05);
%!   assert (numel (h), 4)
%!   assert (p(1), 146.2545, 1E-4)
%!   assert (p(2), 1.1973, 1E-4)
%!   assert (p(3), 0.9999, 5E-5)
%! unwind_protect_cleanup
%!   close (hf);
%! end_unwind_protect
