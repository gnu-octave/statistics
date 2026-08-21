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
## @deftypefn {Private Function} priorFromStruct (@var{S}, @var{ClassNames}, @var{classname})
##
## Resolve a structure form of @qcode{Prior} into a probability row vector.
##
## @var{S} carries a @qcode{ClassNames} field and a @qcode{ClassProbs} field.
## The probabilities are returned in the order of the model's @var{ClassNames}
## rather than the order they were given in, so naming the classes out of
## order still assigns each its own probability.
##
## @end deftypefn

function pr = priorFromStruct (S, ClassNames, classname)

  if (! (isfield (S, 'ClassNames') && isfield (S, 'ClassProbs')))
    error (strcat (classname, ": a structure 'Prior' must have", ...
                   " 'ClassNames' and 'ClassProbs' fields."));
  endif
  sn = S.ClassNames;
  sp = S.ClassProbs;
  if (numel (sn) != numel (sp))
    error (strcat (classname, ": 'ClassNames' and 'ClassProbs' must have", ...
                   " the same number of elements."));
  endif
  K = numel (ClassNames);
  pr = zeros (1, K);
  for i = 1:K
    if (iscellstr (ClassNames))
      j = find (strcmp (sn, ClassNames{i}));
    elseif (ischar (ClassNames))
      j = find (strcmp (sn, ClassNames(i,:)));
    else
      j = find (sn == ClassNames(i));
    endif
    if (isempty (j))
      error (strcat (classname, ": 'ClassNames' in the 'Prior' structure", ...
                     " must name every class of the model."));
    endif
    pr(i) = sp(j(1));
  endfor

endfunction
