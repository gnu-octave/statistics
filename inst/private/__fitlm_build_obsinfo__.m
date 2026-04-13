## Copyright (C) 2026 Jayant Chauhan <0001jayant@gmail.com>
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
## @deftypefn {Private Function} {@var{ObservationInfo} =} __fitlm_build_obsinfo__ (@
##   @var{weights}, @var{excl_mask}, @var{raw_X}, @var{raw_y}, @var{obs_names})
##
## Build the ObservationInfo table for a LinearModel.
##
## Confirmed column names from MATLAB Block 7-P4:
##   {Weights, Excluded, Missing, Subset}
## Subset = true means the row was used in the fit (not excluded, not missing).
## @end deftypefn

function ObservationInfo = __fitlm_build_obsinfo__ (weights, excl_mask, raw_X, raw_y, obs_names)

  n = rows (raw_y);

  ## Missing: any NaN in X rows or y
  nan_in_X = any (isnan (raw_X), 2);
  nan_in_y = isnan (raw_y);
  missing  = logical (nan_in_X | nan_in_y);

  excluded = logical (excl_mask(:));
  subset   = (! excluded) & (! missing);

  ObservationInfo.Weights  = weights(:);
  ObservationInfo.Excluded = excluded;
  ObservationInfo.Missing  = missing;
  ObservationInfo.Subset   = subset;

endfunction
