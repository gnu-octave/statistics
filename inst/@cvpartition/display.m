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
## @deftypefn{Function File} display (@var{C})
## Display a cvpartition object.
##
## @seealso{cvpartition}
## @end deftypefn

## Author: Nir Krakauer

function display (C)

  if nargin != 1
    print_usage ();
  endif
  
  switch C.Type
    case 'kfold'
      str = 'K-fold';
    case 'given'
      str = 'Given';
    case 'holdout'
      str = 'HoldOut';
    case 'leaveout'
      str = 'Leave-One-Out';
    case 'resubstitution'
      str = 'Resubstitution';
    otherwise
      str = 'Unknown-type';
  endswitch 
  
disp([str ' cross validation partition'])
disp(['          N: ' num2str(C.NumObservations)])
disp(['NumTestSets: ' num2str(C.NumTestSets)])
disp(['  TrainSize: ' num2str(C.TrainSize')])
disp(['   TestSize: ' num2str(C.TestSize')])
