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
## @deftypefn {Private Function} [@var{nlogl}, @var{param}, @var{other}] = __paramci__ (@var{pd}, @var{varargin})
##
## Compute the cofindence intervals for the selected parameters of a
## probability distribution object.
##
## @end deftypefn

function ci = __paramci__ (pd, varargin)

  ## Get Distribution specific info
  distname = pd.DistributionCode;
  parnames = pd.ParameterNames(! pd.ParameterIsFixed);

  ## Add defaults and parse optional arguments
  alpha = 0.05;
  param = logical ([1:numel(parnames)]);
  if (mod (numel (varargin), 2) != 0)
    error ("paramci: optional arguments must be in NAME-VALUE pairs.");
  endif
  while (numel (varargin) > 0)
    switch (tolower (varargin{1}))

      case "alpha"
        alpha = varargin{2};
        if (! isscalar (alpha) || ! isnumeric (alpha) || ...
              alpha <= 0 || alpha >= 1)
          error ("paramci: invalid VALUE for 'Alpha' argument.");
        endif

      case "parameter"
        if (iscell (varargin{2}) && numel (varargin{2}) > numel (parnames))
          error ("paramci: invalid VALUE size for 'Parameter' argument.");
        endif
        param = strcmpi (varargin{2}, parnames);
        if (! any (param))
          error ("paramci: unknown distribution parameter.");
        endif

      case {"type", "logflag"}
        printf ("paramci: '%s' argument not supported yet.", varargin{1});

      otherwise
        error ("paramci: invalid NAME for optional argument.");

    endswitch
    varargin([1:2]) = [];
  endwhile

  ## Get confidence intervals for all parameters from selected distribution
  if (strcmpi (distname, "bino"))
    ntrials = pd.N;
    [~, ci] = mle (pd.InputData.data, "distribution", distname, ...
                   "alpha", alpha, "ntrials", ntrials, ...
                   "frequency", pd.InputData.freq);

  elseif (strcmpi (distname, "gp"))
    theta = pd.theta;
    [~, ci] = mle (pd.InputData.data, "distribution", distname, ...
                   "alpha", alpha, "theta", theta, ...
                   "frequency", pd.InputData.freq);

  elseif (strcmpi (distname, "hn"))
    mu = pd.mu;
    [~, ci] = mle (pd.InputData.data, "distribution", distname, ...
                   "alpha", alpha, "mu", mu, "frequency", pd.InputData.freq);

  elseif (! pd.CensoringAllowed)
    [~, ci] = mle (pd.InputData.data, "distribution", distname, ...
                   "alpha", alpha, "frequency", pd.InputData.freq);

  else
    [~, ci] = mle (pd.InputData.data, "distribution", distname, ...
                   "alpha", alpha, "censoring", pd.InputData.cens, ...
                   "frequency", pd.InputData.freq);

  endif

  ## Return ci only for requested parameters
  ci(:, !param) = [];
endfunction
