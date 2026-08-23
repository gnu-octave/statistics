## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Private Function} {@var{O} =} kfoldOpts (@var{args}, @var{validLoss}, @var{classname}, @var{caller}, @var{K})
##
## Parse the options the @code{kfold} methods share.
##
## @var{O} carries @qcode{Mode}, either @qcode{'average'} or
## @qcode{'individual'}; @qcode{Folds}, a row of fold indices defaulting to
## all @var{K} of them; and @qcode{LossFun}, empty unless one was given.
##
## @var{validLoss} lists the losses the caller accepts; pass an empty cell
## for a method that takes no @qcode{'LossFun'} at all, and the option is
## then refused like any other unknown name.
##
## @end deftypefn

function O = kfoldOpts (args, validLoss, classname, caller, K)

  O = struct ('Mode', 'average', 'Folds', 1:K, 'LossFun', '');

  if (mod (numel (args), 2) != 0)
    error (strcat ("%s.%s: optional arguments must be given in", ...
                   " Name-Value pairs."), classname, caller);
  endif

  while (numel (args) > 0)
    switch (lower (args{1}))

      case 'mode'
        O.Mode = args{2};
        if (! (ischar (O.Mode)
               && any (strcmpi (O.Mode, {'average', 'individual'}))))
          error (strcat ("%s.%s: 'Mode' must be either 'average' or", ...
                         " 'individual'."), classname, caller);
        endif
        O.Mode = lower (O.Mode);

      case 'folds'
        O.Folds = args{2};
        if (! (isnumeric (O.Folds) && isreal (O.Folds)
               && isvector (O.Folds) && ! isempty (O.Folds)
               && all (fix (O.Folds) == O.Folds)
               && all (O.Folds >= 1) && all (O.Folds <= K)))
          error (strcat ("%s.%s: 'Folds' must hold integers between 1", ...
                         " and %d."), classname, caller, K);
        endif
        O.Folds = O.Folds(:)';

      case 'lossfun'
        if (isempty (validLoss))
          error (strcat ("%s.%s: invalid parameter name in optional pair", ...
                         " arguments."), classname, caller);
        endif
        O.LossFun = args{2};
        if (! (ischar (O.LossFun) && any (strcmpi (O.LossFun, validLoss))))
          error ("%s.%s: 'LossFun' must be %s.", classname, caller, ...
                 listing (validLoss));
        endif
        O.LossFun = lower (O.LossFun);

      otherwise
        error (strcat ("%s.%s: invalid parameter name in optional pair", ...
                       " arguments."), classname, caller);

    endswitch
    args(1:2) = [];
  endwhile

endfunction

## The accepted values, quoted and listed the way an error message wants
## them, with 'or' before the last.
function s = listing (names)

  q = cellfun (@(n) sprintf ("'%s'", n), names, 'UniformOutput', false);
  if (numel (q) == 1)
    s = q{1};
  elseif (numel (q) == 2)
    s = sprintf ("%s or %s", q{1}, q{2});
  else
    s = sprintf ("%s, or %s", strjoin (q(1:end-1), ', '), q{end});
  endif

endfunction
