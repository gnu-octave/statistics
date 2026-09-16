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
  if (classCount (sn) != numel (sp))
    error (strcat (classname, ": 'ClassNames' and 'ClassProbs' must have", ...
                   " the same number of elements."));
  endif
  ## Textual names are matched whole.  A character matrix holds one name per
  ## row padded with blanks, which cellstr removes, and a categorical or
  ## string array is compared by its text.
  textual = ! (isnumeric (ClassNames) || islogical (ClassNames));
  if (textual)
    names = cellstr (ClassNames);
    if (ischar (sn) || isa (sn, 'categorical') || isa (sn, 'string'))
      sn = cellstr (sn);
    endif
  endif
  K = classCount (ClassNames);
  pr = zeros (1, K);
  for i = 1:K
    if (textual)
      j = find (strcmp (sn, names{i}));
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
