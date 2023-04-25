## Copyright (C) 2014 Nir Krakauer
##
## This file is part of the statistics package for GNU Octave.
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
## @deftypefn  {statistics} {} display (@var{C})
##
## Display a @samp{cvpartition} object, @var{C}.
##
## @seealso{@@cvpartition/cvpartition}
## @end deftypefn

function display (C)

  if (nargin != 1)
    print_usage ();
  endif

  switch C.Type
    case "kfold"
      str = "K-fold";
    case "given"
      str = "Given";
    case "holdout"
      str = "HoldOut";
    case "leaveout"
      str = "Leave-One-Out";
    case "resubstitution"
      str = "Resubstitution";
    otherwise
      str = "Unknown-type";
  endswitch

  disp([str " cross validation partition"])
  disp(["          N: " num2str(C.NumObservations)])
  disp(["NumTestSets: " num2str(C.NumTestSets)])
  disp(["  TrainSize: " num2str(C.TrainSize')])
  disp(["   TestSize: " num2str(C.TestSize')])

endfunction

%!test
%! C = cvpartition (ones (10, 1), "KFold", 5);
%! s = evalc ("display (C)");
%! sout = "K-fold cross validation partition";
%! assert (strcmpi (s(1:length (sout)), sout), true);
%!test
%! C = cvpartition (ones (10, 1), "HoldOut", 5);
%! s = evalc ("display (C)");
%! sout = "HoldOut cross validation partition";
%! assert (strcmpi (s(1:length (sout)), sout), true);
%!test
%! C = cvpartition (ones (10, 1), "LeaveOut", 5);
%! s = evalc ("display (C)");
%! sout = "Leave-One-Out cross validation partition";
%! assert (strcmpi (s(1:length (sout)), sout), true);
%!test
%! C = cvpartition (ones (10, 1), "resubstitution", 5);
%! s = evalc ("display (C)");
%! sout = "Resubstitution cross validation partition";
%! assert (strcmpi (s(1:length (sout)), sout), true);
%!test
%! C = cvpartition (ones (10, 1), "Given", 5);
%! s = evalc ("display (C)");
%! sout = "Given cross validation partition";
%! assert (strcmpi (s(1:length (sout)), sout), true);

%!error<Invalid call to display.  Correct usage is> display ()
