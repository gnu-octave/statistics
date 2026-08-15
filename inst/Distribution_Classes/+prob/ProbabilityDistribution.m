## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftp {statistics} prob.ProbabilityDistribution
##
## Abstract base class of the probability distribution objects.
##
## It holds the behaviour every distribution object shares -- how it is
## displayed, how its parameter confidence intervals are computed, how it is
## plotted and how a profile likelihood is taken -- as protected methods, so
## the 29 distribution classes inherit one implementation and no part of it is
## reachable from outside them.
##
## These helpers were ordinary files in a @file{private} directory until the
## classes moved into the @code{prob} namespace.  Octave does not resolve a
## @file{private} directory from inside a package directory, for a classdef or
## for a plain function, so the only way to keep them out of the public
## interface is to make them protected methods of a shared base.
##
## @end deftp

classdef (Abstract) ProbabilityDistribution

  methods (Access = protected)

    function __disp__ (pd, distname)

      if (isscalar (pd))

        ## Handle special case of prob.PiecewiseLinearDistribution
        if (isa (pd, 'prob.PiecewiseLinearDistribution'))
          ## Print distribution header
          fprintf ("  %s\n\n", pd.DistributionName);
          ## Print parameter values
          for i = 1:numel (pd.x)
            fprintf ("F(%g) = %g\n", pd.x(i), pd.Fx(i));
          endfor
          ## Print truncation interval if applicable
          if (pd.IsTruncated)
            fprintf ("  Truncated to the interval [%g, %g]\n\n", pd.Truncation);
          else
            fprintf ("\n");
          endif

        ## Handle special case of prob.MultinomialDistribution
        elseif (isa (pd, 'prob.MultinomialDistribution'))
          ## Print distribution header
          fprintf ("  %s\n\n", pd.DistributionName);
          ## Print parameter values
          fprintf ("  Probabilities:\n");
          fprintf ("    %0.4f", pd.Probabilities);
          fprintf ("\n\n");
          ## Print truncation interval if applicable
          if (pd.IsTruncated)
            fprintf ("  Truncated to the interval [%g, %g]\n\n", pd.Truncation);
          endif

        ## Handle special case of prob.KernelDistribution
        elseif (isa (pd, 'prob.KernelDistribution'))
          ## Print distribution header
          fprintf ("  %s\n\n", pd.DistributionName);
          ## Print kernel, bandwidth, and support
          fprintf ("    Kernel = %s\n", pd.Kernel);
          fprintf ("    Bandwidth = %g\n", pd.Bandwidth);
          if (ischar (pd.Support.range))
            fprintf ("    Support = %s\n", pd.Support.range);
          else
            fprintf ("    Support = [%g, %g]\n", pd.Support.range);
          endif
          ## Print truncation interval if applicable
          if (pd.IsTruncated)
            fprintf ("  Truncated to the interval [%g, %g]\n\n", pd.Truncation);
          else
            fprintf ("\n");
          endif

        ## Handle all other cases
        else
          ## Get required length for parameter values
          PVlen = max (arrayfun (@(x) numel (sprintf ("%g", x)), ...
                                 pd.ParameterValues));
          PVstr = sprintf ("%%%dg", PVlen);

          ## Prepare template for fitted and not fitted distributions
          pat1 = ['  %+7s = ', PVstr, "   [%g, %g]\n"];
          pat2 = ['  %+7s = ', PVstr, "\n"];

          ## Grad distributions that are non fittable
          if (any (strcmpi (pd.DistributionCode, {'unif', 'tri', 'logu'})))
            fitted = false;
            ParameterIsFixed = true;
          elseif (all (pd.ParameterIsFixed))
            fitted = false;
            ParameterIsFixed = pd.ParameterIsFixed;
          else
            fitted = true;
            ParameterIsFixed = pd.ParameterIsFixed;
          endif

          ## Print distribution header
          fprintf ("  %s\n\n", pd.DistributionName);
          fprintf ("  %s\n", distname);
          ## Print parameter values
          for i = 1:pd.NumParameters
            if (fitted && ! ParameterIsFixed(i))
              fprintf (pat1, pd.ParameterNames{i}, pd.ParameterValues(i), ...
                             pd.ParameterCI(1,i), pd.ParameterCI(2,i));
            else
              fprintf (pat2, pd.ParameterNames{i}, pd.ParameterValues(i));
            endif
          endfor
          ## Print truncation interval if applicable
          if (pd.IsTruncated)
            fprintf ("  Truncated to the interval [%g, %g]\n\n", pd.Truncation);
          else
            fprintf ("\n");
          endif
        endif
      else
        fprintf ("%dx%d %s array\n", size (pd), class (pd));
      endif

    endfunction

    function ci = __paramci__ (pd, varargin)

      ## Get Distribution specific info
      distname = pd.DistributionCode;
      parnames = pd.ParameterNames;

      ## Add defaults and parse optional arguments
      alpha = 0.05;
      param = logical ([1:numel(parnames)]);
      if (mod (numel (varargin), 2) != 0)
        error ("paramci: optional arguments must be in NAME-VALUE pairs.");
      endif
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))

          case 'alpha'
            alpha = varargin{2};
            if (! isscalar (alpha) || ! isnumeric (alpha) || ...
                  alpha <= 0 || alpha >= 1)
              error ("paramci: invalid VALUE for 'Alpha' argument.");
            endif

          case 'parameter'
            if (! isvector (varargin{2}) ||
                ((iscellstr (varargin{2}) || isnumeric (varargin{2})) &&
                numel (varargin{2}) > numel (parnames)))
              error ("paramci: invalid VALUE size for 'Parameter' argument.");
            endif
            if (iscellstr (varargin{2}))
              tmp = cellfun (@(x) strcmpi (x, parnames), varargin{2}, ...
                            'UniformOutput', false);
              param = or (tmp{:});
            elseif (isnumeric (varargin{2}))
              param = ismember (parnames, varargin{2});
            else  # assume it is a character vector
              param = strcmpi (varargin{2}, parnames);
            endif
            if (! any (param))
              error ("paramci: unknown distribution parameter.");
            endif

          case {'type', 'logflag'}
            printf ("paramci: '%s' argument not supported yet.", varargin{1});

          otherwise
            error ("paramci: invalid NAME for optional argument.");

        endswitch
        varargin([1:2]) = [];
      endwhile

      ## Get confidence intervals for all parameters from selected distribution
      if (strcmpi (distname, 'bino'))
        ntrials = pd.N;
        [~, ci] = mle (pd.InputData.data, 'distribution', distname, ...
                       'alpha', alpha, 'ntrials', ntrials, ...
                       'frequency', pd.InputData.freq);

      elseif (strcmpi (distname, 'gp'))
        theta = pd.theta;
        [~, ci] = mle (pd.InputData.data, 'distribution', distname, ...
                       'alpha', alpha, 'theta', theta, ...
                       'frequency', pd.InputData.freq);

      elseif (strcmpi (distname, 'hn'))
        mu = pd.mu;
        [~, ci] = mle (pd.InputData.data, 'distribution', distname, ...
                       'alpha', alpha, 'mu', mu, 'frequency', pd.InputData.freq);

      elseif (! pd.CensoringAllowed)
        [~, ci] = mle (pd.InputData.data, 'distribution', distname, ...
                       'alpha', alpha, 'frequency', pd.InputData.freq);

      else
        [~, ci] = mle (pd.InputData.data, 'distribution', distname, ...
                       'alpha', alpha, 'censoring', pd.InputData.cens, ...
                       'frequency', pd.InputData.freq);

      endif

      ## MLE reports intervals only for the estimated parameters of some
      ## distributions.  Report one column per parameter, as MATLAB does, holding a
      ## fixed parameter at its own value.
      if (columns (ci) != numel (pd.ParameterValues))
        fullci = [pd.ParameterValues; pd.ParameterValues];
        fullci(:, ! pd.ParameterIsFixed) = ci;
        ci = fullci;
      endif

      ## Return ci only for requested parameters
      ci(:, ! param) = [];
    endfunction

    function h = __plot__ (pd, DistType, varargin)

      ## Add defaults (Discrete is passed by the calling method)
      ax = [];
      PlotType = 'pdf';
      Discrete = DistType;

      ## Parse optional arguments
      if (mod (numel (varargin), 2) != 0)
        error ("plot: optional arguments must be in NAME-VALUE pairs.");
      endif
      while (numel (varargin) > 0)
        switch (tolower (varargin{1}))

          case 'plottype'
            ValidTypes = {'pdf', 'cdf', 'probability'};
            try
              selected_T = strcmpi (varargin{2}, ValidTypes);
            catch
              error ("plot: invalid VALUE size for 'Parameter' argument.");
            end_try_catch
            if (! any (selected_T) || sum (selected_T) > 1)
              error ("plot: invalid VALUE for 'PlotType' argument.");
            endif
            PlotType = ValidTypes{selected_T};

          case 'discrete'
            if (! (islogical (varargin{2}) && isscalar (varargin{2})))
              error ("plot: invalid VALUE for 'Discrete' argument.");
            endif
            ## Only for discrete distributions this can be changed by the user
            if (DistType)
              Discrete = varargin{2};
            endif

          case 'parent'
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
      if (strcmpi (PlotType, 'probability'))
        if (! isprop (pd, 'InputData'))
          msg = 'plot: ''probability'' PlotType is not supported for ''%s''.';
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
        case 'pdf'
          h = plot_pdf (pd, ax, DistType, Discrete);
        case 'cdf'
          h = plot_cdf (pd, ax, DistType, Discrete);
        case 'probability'
          h = plot_prob (pd, ax, DistType, Discrete);
      endswitch

    endfunction

    function [nlogl, param, other] = __proflik__ (pd, pnum, varargin)

      ## Default to the first non-fixed parameter
      npvec = find (pd.ParameterIsFixed == false);
      if (nargin < 2 || isempty (pnum))
        pnum = npvec(1);
      endif

      ## Check for non-fixed pnum
      if (! (isnumeric (pnum) && isscalar (pnum) && ismember (pnum, npvec)))
        error (strcat ("proflik: PNUM must be a scalar number", ...
                       " indexing a non-fixed parameter."));
      endif

      ## Add defaults and parse optional arguments
      param = [];
      Display = false;
      while (numel (varargin) > 0)
        if (isnumeric (varargin{1}))
          if (! isvector (varargin{1}))
            error ("proflik: SETPARAM must be a numeric vector.");
          endif
          param = varargin{1};
          varargin(1) = [];
        elseif (ischar (varargin{1}))
          if (strcmpi (varargin{1}, 'display'))
            if (numel (varargin) < 2)
              error ("proflik: missing VALUE for 'Display' argument.");
            endif
            if (! ischar (varargin{2}))
              error ("proflik: invalid VALUE type for 'Display' argument.");
            endif
            if (size (varargin{2}, 1) != 1)
              error ("proflik: invalid VALUE size for 'Display' argument.");
            endif
            if (strcmpi (varargin{2}, 'off'))
              Display = false;
            elseif (strcmpi (varargin{2}, 'on'))
              Display = true;
            else
              error ("proflik: invalid VALUE for 'Display' argument.");
            endif
            varargin([1:2]) = [];
          else
            error ("proflik: invalid NAME for optional arguments.");
          endif
        else
          error ("proflik: invalid optional argument.");
        endif
      endwhile

      ## Optimal parameter values and the free parameters to profile out (the
      ## non-fixed parameters other than the selected one)
      optpar = pd.ParameterValues;
      fname = sprintf ("%slike", pd.DistributionCode);
      freeidx = npvec(npvec != pnum);

      ## Create parameter vector
      pname = pd.ParameterNames{pnum};
      if (isempty (param))
        ## Default range: equally spaced values over the 98% confidence interval,
        ## restricted to the non-fixed range.  MATLAB takes 101 values when the
        ## selected parameter is the only one estimated and 21 when the others must
        ## be profiled out at each of them.
        ci = paramci (pd, "Alpha", 0.02);
        if (any (isnan (ci(:, pnum))))
          error (strcat ("proflik: no confidence interval is defined for '%s',", ...
                         " so the default range cannot be built; supply", ...
                         " SETPARAM instead."), pd.ParameterNames{pnum});
        endif
        lower = max (ci(1, pnum), pd.ParameterRange(1, pnum));
        upper = min (ci(2, pnum), pd.ParameterRange(2, pnum));
        if (isempty (freeidx))
          param = linspace (lower, upper, 101);
        else
          param = linspace (lower, upper, 21);
        endif
      else
        ## Restrict user defined parameter range within non-fixed range
        param(param < pd.ParameterRange(1,pnum)) = [];
        param(param > pd.ParameterRange(2,pnum)) = [];
      endif

      ## Compute the profile log likelihood: at each value of the selected
      ## parameter, maximize the log likelihood over the remaining free parameters
      params = pd.ParameterValues;
      opts = optimset ("Display", "off", "TolX", 1e-6, "TolFun", 1e-6);
      nlogl = zeros (1, numel (param));

      ## Each row of OTHER holds every parameter but the selected one, at the
      ## values maximizing the likelihood; a fixed parameter keeps its own value
      otheridx = [1:numel(params)];
      otheridx(pnum) = [];
      other = zeros (numel (param), numel (otheridx));

      for i = 1:numel (param)
        p0 = params;
        p0(pnum) = param(i);
        if (isempty (freeidx))
          nlogl(i) = - like_value (fname, p0, pd);
        else
          objfun = @(pf) like_free (pf, p0, freeidx, fname, pd);
          [pfhat, fval] = fminsearch (objfun, params(freeidx), opts);
          nlogl(i) = - fval;
          p0(freeidx) = pfhat;
        endif
        other(i,:) = p0(otheridx);
      endfor
      optnll = - like_value (fname, optpar, pd);

      ## Plot the profile log likelihood against the selected parameter, marking
      ## the estimate and the 95% profile-likelihood confidence threshold
      if (Display)
        nll_conf = optnll - 0.5 * chi2inv (0.95, 1);
        plot (optpar(pnum), optnll, 'ok;Estimate;', ...
              param, nlogl, '-r;Profile log likelihood;', ...
              param, repmat (nll_conf, size (param)), ':b;95% confidence;');
        xlabel (pname);
        ylabel ('log likelihood');
        xlim ([param(1), param(end)]);
      endif

    endfunction

  endmethods

endclassdef

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
  PD = {'bino', 'bisa', 'exp', 'gam', 'invg', 'logl', 'logn', ...
        'naka', 'nbin', 'poiss', 'rayl', 'rice', 'wbl'};
  if (strcmpi (pd.DistributionCode, 'beta'))
    lb = xmin = 0;
    ub = xmax = 1;
  elseif (strcmpi (pd.DistributionCode, 'burr'))
    lb = xmin = 0;
    ub = xmax = m + 3 * iqr (pd);
  elseif (any (strcmpi (pd.DistributionCode, PD)))
    lb = max (m - 3 * s, 0);
    xmin = max (m - 3.5 * s, 0);
  elseif (strcmpi (pd.DistributionCode, 'gev'))

  elseif (strcmpi (pd.DistributionCode, 'gp'))

  elseif (strcmpi (pd.DistributionCode, 'hn'))
    lb = max (m - 3 * s, m);
    xmin = max (m - 3.5 * s, m);
  elseif (strcmpi (pd.DistributionCode, 'kernel'))
    ## Clamp the plotting range to a bounded kernel support
    if (ischar (pd.Support.range))
      if (strcmp (pd.Support.range, 'positive'))
        lb = xmin = max (m - 3 * s, 0);
      endif
    else
      lb = xmin = max (m - 3 * s, pd.Support.range(1));
      ub = xmax = min (m + 3 * s, pd.Support.range(2));
    endif
  endif
endfunction

function h = plot_pdf (pd, ax, DistType, Discrete)

  ## Handle special case of multinomial distribution
  if (strcmpi (pd.DistributionCode, 'mn'))
    y = pd.ParameterValues';
    x = [1:numel(y)]';
    if (Discrete)
      h = stem (ax, x, y, 'color', 'b');
    else
      h = plot (ax, x, y, ';;b-');
    endif
    xlim (ax, [0.5, max(x)+0.5]);
    xlabel ('Data');
    ylabel ('Probability');
    return
  endif

  ## Handle special case of piecewise linear distribution
  if (strcmpi (pd.DistributionCode, 'pl'))
    x = pd.ParameterValues(:,1);
    y = pd.ParameterValues(:,2);
    h = plot (ax, x, y, ';;b-');
    xgap = (x(end) - x(1)) * 0.1;
    xlim (ax, [x(1)-xgap, x(end)+xgap]);
    xlabel ('Data');
    ylabel ('Probability');
    return
  endif

  ## Handle special case of triangular distribution
  if (strcmpi (pd.DistributionCode, 'tri'))
    lb = pd.A;
    ub = pd.C;
    xmin = lb - (ub - lb) * 0.1;
    xmax = ub + (ub - lb) * 0.1;
    x = [lb:(ub-lb)/100:ub]';
    y = pdf (pd, x);
    h = plot (ax, x, y, ';;r-', 'linewidth', 2);
    xlim (ax, [xmin, xmax]);
    xlabel ('Data');
    ylabel ('PDF');
    return
  endif

  ## Handle special case of log-uniform and uniform distributions
  if (any (strcmpi (pd.DistributionCode, {'logu', 'unif'})))
    lb = pd.Lower;
    ub = pd.Upper;
    xmin = lb - (ub - lb) * 0.1;
    xmax = ub + (ub - lb) * 0.1;
    x = [lb:(ub-lb)/100:ub]';
    y = pdf (pd, x);
    h = plot (ax, x, y, ';;r-', 'linewidth', 2);
    xlim (ax, [xmin, xmax]);
    xlabel ('Data');
    ylabel ('PDF');
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
      xlim (ax, [min(x)-0.5, max(x)+0.5]);
      h = stem (ax, x, y, 'color', 'r');
      #xlim (ax, [min(x)-0.5, max(x)+0.5]); # before Octave 11 this emits an error
      xlabel ('Data');
      ylabel ('Probability');
    else
      h = plot (ax, x, y, ';;r-', 'linewidth', 2);
      xlim (ax, [xmin, xmax]);
      xlabel ('Data');
      ylabel ('PDF');
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
      h(2) = patch (ax, bincenter, binsize, 1, 'facecolor', 'b');
      xlim (ax, [xmin, xmax]);
      hold on;
      if (Discrete)
        h(1) = stem (ax, x, y, 'color', 'r');
      else
        h(1) = plot (ax, x, y, ';;r-');
      endif
      xlabel ('Data');
      ylabel ('Probability');
      hold off;
    else
      h(2) = bar (ax, bincenter, binsize, 1, 'facecolor', 'b');
      hold on;
      h(1) = plot (ax, x, y, ';;r-', 'linewidth', 2);
      xlim (ax, [xmin, xmax]);
      xlabel ('Data');
      ylabel ('PDF');
      hold off;
    endif
  endif

endfunction

function h = plot_cdf (pd, ax, DistType, Discrete)

  ## Handle special case of multinomial distribution
  if (strcmpi (pd.DistributionCode, 'mn'))
    y = pd.ParameterValues';
    x = [1:numel(y)]';
    xlim (ax, [0.5, max(x)+0.5]);
    if (Discrete)
      h = stem (ax, x, y, 'color', 'b');
    else
      h = plot (ax, x, y, ';;b-');
    endif
    xlabel ('Data');
    ylabel ('Probability');
    return
  endif

  ## Handle special case of piecewise linear distribution
  if (strcmpi (pd.DistributionCode, 'pl'))
    x = pd.ParameterValues(:,1);
    y = pd.ParameterValues(:,2);
    h = plot (ax, x, y, ';;b-');
    xgap = (x(end) - x(1)) * 0.1;
    xlim (ax, [x(1)-xgap, x(end)+xgap]);
    xlabel ('Data');
    ylabel ('Probability');
    return
  endif

  ## Handle special case of triangular distribution
  if (strcmpi (pd.DistributionCode, 'tri'))
    lb = pd.A;
    ub = pd.C;
    xmin = lb - (ub - lb) * 0.1;
    xmax = ub + (ub - lb) * 0.1;
    x = [lb:(ub-lb)/100:ub]';
    y = pdf (pd, x);
    h = plot (ax, x, y, ';;r-', 'linewidth', 2);
    xlim (ax, [xmin, xmax]);
    xlabel ('Data');
    ylabel ('PDF');
    return
  endif

  ## Handle special case of log-uniform and uniform distributions
  if (any (strcmpi (pd.DistributionCode, {'logu', 'unif'})))
    lb = pd.Lower;
    ub = pd.Upper;
    xmin = lb - (ub - lb) * 0.1;
    xmax = ub + (ub - lb) * 0.1;
    x = [lb:(ub-lb)/100:ub]';
    y = pdf (pd, x);
    h = plot (ax, x, y, ';;r-', 'linewidth', 2);
    xlim (ax, [xmin, xmax]);
    xlabel ('Data');
    ylabel ('PDF');
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
      h = stairs (ax, x, p, 'color', 'r');
      xlim (ax, [lb-0.5, ub+0.5]);
      ylim (ax, [0, 1]);
      xlabel ('Data');
      ylabel ('CDF');
    else
      h = plot (ax, x, p, ';;r-', 'linewidth', 2);
      xlim (ax, [xmin, xmax]);
      ylim (ax, [0, 1]);
      xlabel ('Data');
      ylabel ('CDF');
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
      h(2) = plot (ax, xCDF, yCDF, ';;b-');
      xlim (ax, [xmin, xmax]);
      ylim (ax, [0, 1]);
      hold on;
      if (Discrete)
        h(1) = stem (ax, x, p, 'color', 'r');
      else
        h(1) = plot (ax, x, p, ';;r-');
      endif
      xlabel ('Data');
      ylabel ('CDF');
      hold off;
    else
      h(2) = plot (ax, xCDF, yCDF, ';;b-');
      xlim (ax, [xmin, xmax]);
      ylim (ax, [0, 1]);
      hold on;
      h(1) = plot (ax, x, p, ';;r-', 'linewidth', 2);
      xlabel ('Data');
      ylabel ('CDF');
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
  h(2) = line (ax, X, Y, 'LineStyle', '-.', 'Marker', 'none', 'color', 'red');
  hold on;
  h(1) = plot (ax, x, y, 'LineStyle', 'none', 'Marker', '+', 'color', 'blue');

  ## Plot labels
  ylabel 'Probability'
  xlabel 'Data'

  ## Plot grid
  p = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.25, 0.5, ...
       0.75, 0.90, 0.95, 0.98, 0.99, 0.995, 0.999];
  label = {'0.001', '0.005', '0.01', '0.02', '0.05', '0.10', '0.25', '0.50', ...
           '0.75', '0.90', '0.95', '0.98', '0.99', '0.995', '0.999'};
  tick = icdf (pd, p);
  set (ax, 'ytick', tick, 'yticklabel', label);

  ## Compute plot boundaries
  [~, ~, xmin, xmax] = compute_boundaries (pd);
  ## Set view range with a bit of space around data
  ymin = icdf (pd, 0.25 ./ n);
  ymax = icdf (pd, (n - 0.25) ./ n);
  set (ax, 'ylim', [ymin, ymax], 'xlim', [xmin, xmax]);
  grid (ax, 'on');
  box (ax, 'off');
  hold off;

endfunction

## Evaluate the family negative log likelihood at a full parameter vector.
## The family function takes (params, x, censor, freq) or (params, x, freq);
## dispatch on what it declares rather than on the distribution's own
## CensoringAllowed, since Burr does not allow censoring yet still takes the
## four-argument form, and passing the frequencies positionally into the
## censoring slot marks every observation censored.
function nll = like_value (fname, params, pd)
  if (nargin (fname) > 3)
    nll = feval (fname, params, pd.InputData.data, ...
                 pd.InputData.cens, pd.InputData.freq);
  else
    nll = feval (fname, params, pd.InputData.data, pd.InputData.freq);
  endif
endfunction

## Negative log likelihood as a function of the free parameters only.  The
## search is unconstrained, so parameters outside their own range are rejected
## here: the family likelihood may well return a finite value there, and for
## some families it grows without bound as they run away.
function nll = like_free (pf, p0, freeidx, fname, pd)
  if (any (pf(:)' < pd.ParameterRange(1,freeidx)) ||
      any (pf(:)' > pd.ParameterRange(2,freeidx)))
    nll = Inf;
    return;
  endif
  p = p0;
  p(freeidx) = pf;
  nll = like_value (fname, p, pd);
  if (! isreal (nll) || ! isfinite (nll))
    nll = Inf;
  endif
endfunction
