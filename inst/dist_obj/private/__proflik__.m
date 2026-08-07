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
## @deftypefn {Private Function} {[@var{nlogl}, @var{param}, @var{other}] =} __proflik__ (@var{pd}, @var{pnum}, @var{varargin})
##
## Compute the negative log likelihood for a range of the selected parameter of
## a probability distribution object.
##
## @end deftypefn

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
