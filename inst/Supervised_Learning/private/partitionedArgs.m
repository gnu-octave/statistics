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
## @deftypefn {Private Function} {[@var{P}, @var{rest}] =} partitionedArgs (@var{args}, @var{classname})
##
## Split the options a cross-validated model owns from the ones its folds
## are fitted with.
##
## @var{P} is a structure with fields @qcode{ClassNames}, @qcode{Prior},
## @qcode{Cost}, @qcode{Weights}, @qcode{ScoreTransform} and
## @qcode{ResponseTransform}, each empty when it was not given.  @var{rest}
## is everything else, in the order it was given, including the option that
## says how to partition.
##
## Those six belong to the parent because they describe the data or the
## reported scores rather than the fit: the classes, the prior and the cost
## are resolved once over the whole data and handed down, the weights are
## sliced per fold, and a transform is applied once to the assembled scores
## with the folds left carrying none.
##
## @end deftypefn

function [P, rest] = partitionedArgs (args, classname)

  P = struct ('ClassNames', [], 'Prior', [], 'Cost', [], 'Weights', [], ...
              'ScoreTransform', [], 'ResponseTransform', []);
  rest = {};

  while (numel (args) > 0)
    if (numel (args) < 2)
      error (strcat ("%s: optional arguments must be given in Name-Value", ...
                     " pairs."), classname);
    endif
    switch (lower (args{1}))
      case 'classnames'
        P.ClassNames = args{2};
      case 'prior'
        P.Prior = args{2};
      case 'cost'
        P.Cost = args{2};
      case 'weights'
        P.Weights = args{2};
      case 'scoretransform'
        P.ScoreTransform = args{2};
      case 'responsetransform'
        P.ResponseTransform = args{2};
      otherwise
        rest = [rest, args(1:2)];
    endswitch
    args(1:2) = [];
  endwhile

endfunction
