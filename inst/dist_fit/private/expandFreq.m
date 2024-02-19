## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn  {statistics} {[@var{xExp}, @var{cExp}] =} expandFreq (@var{freq}, @var{x}, @var{censor})
##
## Expand data and censoring vectors according to frequencies.
##
## @end deftypefn

function [xExp, cExp] = expandFreq (freq, x, c)
  ## Create censoring data if not given
  if (isempty (c))
    c = zeros (size (x));
  endif
  if (isempty (freq))
    xExp = x;
    cExp = c;
    return
  endif
  ## Check  X and FREQ vectors match (CENSOR vector is expected to match)
  if (! isequal (size (x), size (freq)))
    error ("expandFreq: X and FREQ vector mismatch.");
  endif
  ## Remove any NaNs
  is_nan = isnan (freq) | isnan (x) | isnan (c);
  freq(is_nan) = [];
  x(is_nan) = [];
  c(is_nan) = [];
  ## Check for valid values in FREQ
  if (any (freq != round (freq)) || any (freq < 0))
    error ("expandFreq: FREQ vector must contain non-negative integer values.");
  endif
  ## Compute expansion
  xExp = [];
  cExp = [];
  for i = 1:numel (freq)
    xExp = [xExp, repmat(x(i), 1, freq(i))];
    cExp = [cExp, repmat(c(i), 1, freq(i))];
  endfor
endfunction
