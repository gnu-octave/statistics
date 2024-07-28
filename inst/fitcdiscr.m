## Copyright (C) 2024 Ruchika Sonagote <ruchikasonagote2003@gmail.com>
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

function obj = fitcdiscr (X, Y, varargin)

  ## Check input parameters
  if (nargin < 2)
    error ("fitcdiscr: too few arguments.");
  endif

  if (mod (nargin, 2) != 0)
    error ("fitcdiscr: name-value arguments must be in pairs.");
  endif

  ## Check predictor data and labels have equal rows
  if (rows (X) != rows (Y))
    error ("fitcdiscr: number of rows in X and Y must be equal.");
  endif

  ## Parse arguments to class def function
  obj = ClassificationDiscriminant (X, Y, varargin{:});

endfunction