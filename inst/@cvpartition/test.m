## Copyright (C) 2014 Nir Krakauer
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn{Function File}{@var{inds} =} test (@var{C}, [@var{i}])
## Return logical vector for testing-subset indices from a cvpartition object.
##
## @var{C} should be a cvpartition object. @var{i} is the fold index (default is 1).
##
## @seealso{cvpartition, @@cvpartition/training}
## @end deftypefn

## Author: Nir Krakauer

function inds = test (C, i = [])

  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  if nargin < 2 || isempty (i)
    i = 1;
  endif

  switch C.Type
    case {'kfold' 'given'}
      inds = C.inds == i;
    case 'holdout'
      inds = C.inds;
    case 'leaveout'
      inds = zeros(C.NumObservations, 1, "logical");
      inds(i) = true;
    case 'resubstitution'
      inds = ones(C.NumObservations, 1, "logical");
  endswitch
