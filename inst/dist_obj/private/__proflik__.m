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
## @deftypefn {Private Function} [@var{nlogl}, @var{param}, @var{other}] = __proflik__ (@var{pd}, @var{pnum}, @var{varargin})
##
## Compute the negative log likelihood for a range of the selected parameter of
## a probability distribution object.
##
## @end deftypefn

function [nlogl, param] = __proflik__ (pd, pnum, varargin)

  ## Check for valid pnum
  npvec = [1:pd.NumParameters];
  if (! (isnumeric (pnum) && isscalar (pnum) && ismember (pnum, npvec)))
    error ("proflik: PNUM must be a scalar number indexing a valid parameter.");
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
      if (strcmpi (varargin{1}, "display"))
        if (numel (varargin) < 2)
          error ("proflik: missing VALUE for 'Display' argument.");
        endif
        if (! ischar (varargin{2}))
          error ("proflik: invalid VALUE type for 'Display' argument.");
        endif
        if (size (varargin{2}, 1) != 1)
          error ("proflik: invalid VALUE size for 'Display' argument.");
        endif
        if (strcmpi (varargin{2}, "off"))
          Display = false;
        elseif (strcmpi (varargin{2}, "on"))
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

  ## Get fitted parameter CI and restrict for valid range
  pname = pd.ParameterNames{pnum};
  lower = pd.ParameterCI(1, pnum);
  if (lower < pd.ParameterRange(1,pnum))
    lower = pd.ParameterRange(1,pnum);
  endif
  upper = pd.ParameterCI(2, pnum);
  if (upper > pd.ParameterRange(2,pnum))
    upper = pd.ParameterRange(2,pnum);
  endif
  ## Create parameter vector
  if (isempty (param))
    param = [lower:(upper-lower)/100:upper];
  else
    ## Restrict user defined parameter range within valid range
    param(param < pd.ParameterRange(1,pnum)) = [];
    param(param > pd.ParameterRange(2,pnum)) = [];
  endif

  ## Get nlogl estimate and optimal parameter values
  optnll = negloglik (pd);
  optpar = pd.ParameterValues;

  ## Compute negative log likelihood for the given parameter range by
  ## calling the appropriate likelihood function
  params = pd.ParameterValues;
  for i = 1:numel (param)
    params(pnum) = param(i);
    fname = sprintf ("%slike", pd.DistributionCode);
    nlogl(i) = - feval (fname, params, pd.InputData.data, ...
                        pd.InputData.cens, pd.InputData.freq);
  endfor

  ## Plot the loglikelihood values against the range of the selected parameter
  if (Display)
    ## Get 95% confidence
    params(pnum) = lower;
    nll_conf = - feval (fname, params, pd.InputData.data, ...
                        pd.InputData.cens, pd.InputData.freq);
    plot (optpar(pnum), optnll, "ok;Estimate;", ...
          param, nlogl, "-r;Exact log likelihood;", ...
          param, repmat (nll_conf, size (param)), ":b;95% confidence;");
    xlabel (pname);
    ylabel ("log likelihood");
    xlim ([param(1), param(end)]);
  endif

endfunction
