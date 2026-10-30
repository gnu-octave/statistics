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
## @deftypefn {Private Function} {@var{tf} =} zeroPriorNames (@var{ClassNames}, @var{Prior}, @var{Cost})
##
## Which of the classes named in @qcode{'ClassNames'} the user gave a prior
## probability of 0, or a row of misclassification costs that is all zeros.
##
## Such a class is removed before fitting, so it may be absent from the
## response, where any other named class must be present.  @var{Prior} and
## @var{Cost} are the options as the user gave them, before they are
## resolved: a numeric vector or a structure for the prior, a matrix or a
## structure for the cost, or anything else, which marks no class.
## @var{tf} is a logical column, one element per name.
## @end deftypefn

function tf = zeroPriorNames (ClassNames, Prior, Cost)

  K = classCount (ClassNames);
  tf = false (K, 1);
  if (K == 0)
    return;
  endif

  if (isstruct (Prior) && isfield (Prior, 'ClassNames')
      && isfield (Prior, 'ClassProbs'))
    try
      Prior = priorFromStruct (Prior, ClassNames, '');
    catch
      Prior = [];
    end_try_catch
  endif
  if (isnumeric (Prior) && isreal (Prior) && numel (Prior) == K)
    tf |= (Prior(:) == 0);
  endif

  if (isstruct (Cost))
    [Cost, errmsg] = costMatrix (Cost, ClassNames);
    if (! isempty (errmsg))
      Cost = [];
    endif
  endif
  if (isnumeric (Cost) && isreal (Cost) && isequal (size (Cost), [K, K])
      && K > 1)
    tf |= all (Cost == 0, 2);
  endif

endfunction
