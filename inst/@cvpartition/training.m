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
## @deftypefn  {statistics} {@var{inds} =} training (@var{C}, [@var{i}])
##
## Return logical vector for training-subset indices from a @samp{cvpartition}
## object, @var{C}.  @var{i} is the fold index (default is 1).
##
## @seealso{@cvpartition/cvpartition, @cvpartition/test}
## @end deftypefn

function inds = training (C, i = [])

  if (nargin < 1 || nargin > 2)
    print_usage ();
  endif

  if (nargin < 2 || isempty (i))
    i = 1;
  endif

  switch (C.Type)
    case  {"kfold", "given"}
      inds = C.inds != i;
    case "holdout"
      inds = ! C.inds;
    case "leaveout"
      inds = ones (C.NumObservations, 1, "logical");
      inds(i) = false;
    case "resubstitution"
      inds = ones (C.NumObservations, 1, "logical");
  endswitch

endfunction


%!shared C
%! C = cvpartition (ones (10, 1), "KFold", 5);
%!assert (training (C, 1), logical ([0 0 1 1 1 1 1 1 1 1]'))
%!assert (training (C, 2), logical ([1 1 0 0 1 1 1 1 1 1]'))
%!assert (training (C, 3), logical ([1 1 1 1 0 0 1 1 1 1]'))
%!assert (training (C, 4), logical ([1 1 1 1 1 1 0 0 1 1]'))
%!assert (training (C, 5), logical ([1 1 1 1 1 1 1 1 0 0]'))
%!test
%! C = set (C, "inds", [1 2 2 2 3 4 3 4 5 5]');
%!assert (training (C), logical ([0 1 1 1 1 1 1 1 1 1]'))
%!assert (training (C, 2), logical ([1 0 0 0 1 1 1 1 1 1]'))
%!assert (training (C, 3), logical ([1 1 1 1 0 1 0 1 1 1]'))
