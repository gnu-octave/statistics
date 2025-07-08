## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Private Function} [@var{h}] = __plot__ (@var{pd}, @var{Discrete}, @var{varargin})
##
## Plot a probability distribution object.
##
## @end deftypefn

function h = __plot__ (pd, DistType, varargin)

  ## Add defaults (Discrete is passed by the calling method)
  ax = [];
  PlotType = "pdf";
  Discrete = DistType;

  ## Parse optional arguments
  if (mod (numel (varargin), 2) != 0)
    error ("plot: optional arguments must be in NAME-VALUE pairs.");
  endif
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))

      case "plottype"
        ValidTypes = {"pdf", "cdf", "probability"};
        try
          selected_T = strcmpi (varargin{2}, ValidTypes);
        catch
          error ("plot: invalid VALUE size for 'Parameter' argument.");
        end_try_catch
        if (! any (selected_T) || sum (selected_T) > 1)
          error ("plot: invalid VALUE for 'PlotType' argument.");
        endif
        PlotType = ValidTypes{selected_T};

      case "discrete"
        if (! (islogical (varargin{2}) && isscalar (varargin{2})))
          error ("plot: invalid VALUE for 'Discrete' argument.");
        endif
        ## Only for discrete distributions this can be changed by the user
        if (DistType)
          Discrete = varargin{2};
        endif

      case "parent"
        if (! isaxes (varargin{2}))
          error ("plot: invalid VALUE for 'Parent' argument.");
        endif
        ax = varargin{2};

      otherwise
        error ("plot: invalid NAME for optional argument.");

    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Check for invalid cases of probability type before creating new axes
  if (strcmpi (PlotType, "probability"))
    if (! isprop (pd, "InputData"))
      msg = "plot: 'probability' PlotType is not supported for '%s'.";
      error (sprintf (msg, pd.DistributionName));
    endif
    if (isempty (pd.InputData))
      error ("plot: no fitted DATA to plot a probability plot.");
    endif
  endif

  ## Get current axes or create new ones
  if (isempty (ax))
    ax = gca ();
  endif

  ## Switch to PlotType
  switch (PlotType)
    case "pdf"
      h = plot_pdf (pd, ax, DistType, Discrete);
    case "cdf"
      h = plot_cdf (pd, ax, DistType, Discrete);
    case "probability"
      h = plot_prob (pd, ax, DistType, Discrete);
  endswitch

endfunction

function x = expand_freq (data, freq)
  x = [];
  for i = 1:numel (freq)
    x = [x, repmat(data(i), 1, freq(i))];
  endfor
endfunction

function [lb, ub, xmin, xmax] = compute_boundaries (pd)
  ## Compute moments to determine plot boundaries
  m = mean (pd);
  s = std (pd);
  lb = m - 3 * s;
  ub = m + 3 * s;
  xmin = m - 3.5 * s;
  xmax = m + 3.5 * s;
  ## Fix boundaries for specific distributions
  PD = {"bino", "bisa", "exp", "gam", "invg", "logl", ...
        "logn", "naka", "nbin", "poiss", "rayl", "rice", "wbl"};
  if (strcmpi (pd.DistributionCode, "beta"))
    lb = xmin = 0;
    ub = xmax = 1;
  elseif (strcmpi (pd.DistributionCode, "burr"))
    lb = xmin = 0;
    ub = xmax = m + 3 * iqr (pd);
  elseif (any (strcmpi (pd.DistributionCode, PD)))
    lb = max (m - 3 * s, 0);
    xmin = max (m - 3.5 * s, 0);
  elseif (strcmpi (pd.DistributionCode, "gev"))

  elseif (strcmpi (pd.DistributionCode, "gp"))

  elseif (strcmpi (pd.DistributionCode, "hn"))
    lb = max (m - 3 * s, m);
    xmin = max (m - 3.5 * s, m);
  endif
endfunction

function h = plot_pdf (pd, ax, DistType, Discrete)

  ## Handle special case of multinomial distribution
  if (strcmpi (pd.DistributionCode, "mn"))
    y = pd.ParameterValues';
    x = [1:numel(y)]';
    if (Discrete)
      h = stem (ax, x, y, "color", "b");
    else
      h = plot (ax, x, y, ";;b-");
    endif
    xlim ([0.5, max(x)+0.5]);
    xlabel ("Data");
    ylabel ("Probability");
    return
  endif

  ## Handle special case of piecewise linear distribution
  if (strcmpi (pd.DistributionCode, "pl"))
    x = pd.ParameterValues(:,1);
    y = pd.ParameterValues(:,2);
    h = plot (ax, x, y, ";;b-");
    xgap = (x(end) - x(1)) * 0.1;
    xlim ([x(1)-xgap, x(end)+xgap]);
    xlabel ("Data");
    ylabel ("Probability");
    return
  endif

  ## Handle special case of triangular distribution
  if (strcmpi (pd.DistributionCode, "tri"))
    lb = pd.A;
    ub = pd.C;
    xmin = lb - (ub - lb) * 0.1;
    xmax = ub + (ub - lb) * 0.1;
    x = [lb:(ub-lb)/100:ub]';
    y = pdf (pd, x);
    h = plot (ax, x, y, ";;r-", "linewidth", 2);
    xlim ([xmin, xmax]);
    xlabel ("Data");
    ylabel ("PDF");
    return
  endif

  ## Handle special case of log-uniform and uniform distributions
  if (any (strcmpi (pd.DistributionCode, {"logu", "unif"})))
    lb = pd.Lower;
    ub = pd.Upper;
    xmin = lb - (ub - lb) * 0.1;
    xmax = ub + (ub - lb) * 0.1;
    x = [lb:(ub-lb)/100:ub]';
    y = pdf (pd, x);
    h = plot (ax, x, y, ";;r-", "linewidth", 2);
    xlim ([xmin, xmax]);
    xlabel ("Data");
    ylabel ("PDF");
    return
  endif

  ## Check for fitted distribution
  if (isempty (pd.InputData)) # fixed parameters, no data

    ## Compute plot boundaries
    [lb, ub, xmin, xmax] = compute_boundaries (pd);

    ## Compute stem or line for PDF
    if (DistType)
      x = [floor(lb):ceil(ub)]';
      y = pdf (pd, x);
    else
      x = [lb:(ub-lb)/100:ub]';
      y = pdf (pd, x);
    endif

    ## Plot
    if (Discrete)
      h = stem (ax, x, y, "color", "r");
      xlim ([min(y)-0.5, max(y)+0.5]);
      xlabel ("Data");
      ylabel ("Probability");
    else
      h = plot (ax, x, y, ";;r-", "linewidth", 2);
      xlim ([xmin, xmax]);
      xlabel ("Data");
      ylabel ("PDF");
    endif

  else # fitted distribution, data available

    ## Expand frequency vector (if necessary)
    if (any (pd.InputData.freq != 1))
      x = expand_freq (pd.InputData.data, pd.InputData.freq);
    else
      x = pd.InputData.data;
    endif

    ## Keep data within plotting boundaries
    [lb, ub, xmin, xmax] = compute_boundaries (pd);
    x(x < lb | x > ub) = [];

    ## Compute the patch or histogram for data
    xsize = numel (x);
    if (DistType)
      binwidth = 1;
      xmin = min (x) - 1;
      xmax = max (x) + 1;
      [binsize, bincenter] = hist (x, [xmin:xmax]);
    else
      nbins = ceil (sqrt (xsize));
      [binsize, bincenter] = hist (x, nbins);
      binwidth = max (diff (bincenter));
      xmin = min (x) - binwidth / 2;
      xmax = max (x) + binwidth / 2;
    endif

    ## Compute stem or line for PDF
    if (Discrete)
      x = [min(x):max(x)]';
      y = pdf (pd, x);
    else
      x = [xmin:(xmax-xmin)/100:xmax]';
      y = pdf (pd, x);
    endif

    ## Normalize density line
    y = xsize * y * binwidth;

    ## Plot
    if (DistType)
      h(2) = patch (ax, bincenter, binsize, 1, "facecolor", "b");
      hold on;
      if (Discrete)
        h(1) = stem (ax, x, y, "color", "r");
      else
        h(1) = plot (ax, x, y, ";;r-");
      endif
      xlim ([xmin, xmax]);
      xlabel ("Data");
      ylabel ("Probability");
      hold off;
    else
      h(2) = bar (ax, bincenter, binsize, 1, "facecolor", "b");
      hold on;
      h(1) = plot (ax, x, y, ";;r-", "linewidth", 2);
      xlim ([xmin, xmax]);
      xlabel ("Data");
      ylabel ("PDF");
      hold off;
    endif
  endif

endfunction

function h = plot_cdf (pd, ax, DistType, Discrete)

  ## Handle special case of multinomial distribution
  if (strcmpi (pd.DistributionCode, "mn"))
    y = pd.ParameterValues';
    x = [1:numel(y)]';
    if (Discrete)
      h = stem (ax, x, y, "color", "b");
    else
      h = plot (ax, x, y, ";;b-");
    endif
    xlim ([0.5, max(x)+0.5]);
    xlabel ("Data");
    ylabel ("Probability");
    return
  endif

  ## Handle special case of piecewise linear distribution
  if (strcmpi (pd.DistributionCode, "pl"))
    x = pd.ParameterValues(:,1);
    y = pd.ParameterValues(:,2);
    h = plot (ax, x, y, ";;b-");
    xgap = (x(end) - x(1)) * 0.1;
    xlim ([x(1)-xgap, x(end)+xgap]);
    xlabel ("Data");
    ylabel ("Probability");
    return
  endif

  ## Handle special case of triangular distribution
  if (strcmpi (pd.DistributionCode, "tri"))
    lb = pd.A;
    ub = pd.C;
    xmin = lb - (ub - lb) * 0.1;
    xmax = ub + (ub - lb) * 0.1;
    x = [lb:(ub-lb)/100:ub]';
    y = pdf (pd, x);
    h = plot (ax, x, y, ";;r-", "linewidth", 2);
    xlim ([xmin, xmax]);
    xlabel ("Data");
    ylabel ("PDF");
    return
  endif

  ## Handle special case of log-uniform and uniform distributions
  if (any (strcmpi (pd.DistributionCode, {"logu", "unif"})))
    lb = pd.Lower;
    ub = pd.Upper;
    xmin = lb - (ub - lb) * 0.1;
    xmax = ub + (ub - lb) * 0.1;
    x = [lb:(ub-lb)/100:ub]';
    y = pdf (pd, x);
    h = plot (ax, x, y, ";;r-", "linewidth", 2);
    xlim ([xmin, xmax]);
    xlabel ("Data");
    ylabel ("PDF");
    return
  endif

  ## Compute plot boundaries
  [lb, ub, xmin, xmax] = compute_boundaries (pd);

  ## Check for fitted distribution
  if (isempty (pd.InputData)) # fixed parameters, no data

    ## Compute stem or line for PDF
    if (DistType)
      x = [floor(lb):ceil(ub)]';
      p = cdf (pd, x);
    else
      x = [lb:(ub-lb)/100:ub]';
      p = cdf (pd, x);
    endif

    ## Plot
    if (Discrete)
      h = stairs (ax, x, p, "color", "r");
      xlim ([lb-0.5, ub+0.5]);
      ylim ([0, 1]);
      xlabel ("Data");
      ylabel ("CDF");
    else
      h = plot (ax, x, p, ";;r-", "linewidth", 2);
      xlim ([xmin, xmax]);
      ylim ([0, 1]);
      xlabel ("Data");
      ylabel ("CDF");
    endif

  else # fitted distribution, data available

    ## Expand frequency vector (if necessary)
    if (any (pd.InputData.freq != 1))
      x = expand_freq (pd.InputData.data, pd.InputData.freq);
    else
      x = pd.InputData.data;
    endif

    ## Compute the stairs for data
    [yy, xx, ~, ~, eid] = cdfcalc (x);
    n = length (xx);
    ## Create vectors for plotting
    nidx = reshape (repmat (1:n, 2, 1), 2*n, 1);
    xCDF = [-Inf; xx(nidx); Inf];
    yCDF = [0; 0; yy(1+nidx)];

    ## Compute stairs or line for CDF
    if (DistType)
      x = [min(x):max(x)]';
      p = cdf (pd, x);
    else
      x = [xmin:(xmax-xmin)/100:xmax]';
      p = cdf (pd, x);
    endif

    ## Plot
    if (DistType)
      h(2) = plot (ax, xCDF, yCDF, ";;b-");
      hold on;
      if (Discrete)
        h(1) = stem (ax, x, p, "color", "r");
      else
        h(1) = plot (ax, x, p, ";;r-");
      endif
      xlim ([xmin, xmax]);
      ylim ([0, 1]);
      xlabel ("Data");
      ylabel ("CDF");
      hold off;
    else
      h(2) = plot (ax, xCDF, yCDF, ";;b-");
      hold on;
      h(1) = plot (ax, x, p, ";;r-", "linewidth", 2);
      xlim ([xmin, xmax]);
      ylim ([0, 1]);
      xlabel ("Data");
      ylabel ("CDF");
      hold off;
    endif
  endif

endfunction

function h = plot_prob (pd, ax, DistType, Discrete)

  ## Expand frequency vector (if necessary)
  if (any (pd.InputData.freq != 1))
    x = expand_freq (pd.InputData.data, pd.InputData.freq);
  else
    x = pd.InputData.data;
  endif

  ## Compute the probabilities for data
  n = rows (x);
  y = icdf (pd, ([1:n]' - 0.5) / n);
  x = sort (x);

  ## Plot reference line
  X = Y = [x(1); x(end)];
  h(2) = line (ax, X, Y, "LineStyle", "-.", "Marker", "none", "color", "red");
  hold on;
  h(1) = plot (ax, x, y, "LineStyle", "none", "Marker", "+", "color", "blue");

  ## Plot labels
  ylabel "Probability"
  xlabel "Data"

  ## Plot grid
  p = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.25, 0.5, ...
       0.75, 0.90, 0.95, 0.98, 0.99, 0.995, 0.999];
  label = {"0.001", "0.005", "0.01", "0.02", "0.05", "0.10", "0.25", "0.50", ...
           "0.75", "0.90", "0.95", "0.98", "0.99", "0.995", "0.999"};
  tick = icdf (pd, p);
  set (ax, "ytick", tick, "yticklabel", label);

  ## Compute plot boundaries
  [~, ~, xmin, xmax] = compute_boundaries (pd);
  ## Set view range with a bit of space around data
  ymin = icdf (pd, 0.25 ./ n);
  ymax = icdf (pd, (n - 0.25) ./ n);
  set (ax, "ylim", [ymin, ymax], "xlim", [xmin, xmax]);
  grid (ax, "on");
  box (ax, "off");
  hold off;

endfunction
