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
## FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
## more details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{y0} =} missingResponse (@var{Y}, @var{W})
##
## The response a regression model predicts for a row missing a predictor.
##
## @var{y0} is the weighted lower median of the training response @var{Y}
## under the observation weights @var{W}: the smallest response at which the
## weight of the responses up to it reaches half the total.  With equal
## weights and an even count it is the lower of the two middle values.
## Measured on MATLAB R2024a, whose linear, kernel, support vector, Gaussian
## process and neural network regressions all predict it.  A missing
## response takes no part.
##
## @end deftypefn

function y0 = missingResponse (Y, W)

  Y = Y(:);
  W = W(:);
  keep = ! (isnan (Y) | isnan (W));
  Y = Y(keep);
  W = W(keep);
  if (isempty (Y) || ! (sum (W) > 0))
    y0 = NaN;
    return;
  endif
  [Y, o] = sort (Y);
  cw = cumsum (W(o));
  ## Equal weights reach one half exactly, and the cumulative sum of fifty
  ## hundredths falls short of it in the last bits.
  i = find (cw >= cw(end) / 2 * (1 - 1e-12), 1);
  y0 = Y(i);

endfunction
