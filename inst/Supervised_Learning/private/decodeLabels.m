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
## FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Private Function} {@var{data} =} decodeLabels (@var{data})
##
## Rebuild the categorical and string arrays of a loaded model.
##
## @var{data} is the structure @code{load} returns for a saved model.  A field
## written by @code{encodeLabels} is turned back into its categorical or string
## array.  So is a field written by version 1.9.2, which saved such an array as
## the structure @code{save} falls back to for a classdef object: fields
## @qcode{cats}, @qcode{code}, @qcode{isMissing}, @qcode{isOrdinal} and
## @qcode{isProtected} for a categorical array, and @qcode{isMissing} and
## @qcode{strs} for a string array.  A field holding a structure array is
## searched the same way, element by element, so class names a model keeps
## inside a structure are rebuilt too.  Every other field is left as it is.
##
## @seealso{encodeLabels}
## @end deftypefn

function data = decodeLabels (data)

  fn = fieldnames (data);
  for e = 1:numel (data)
    for j = 1:numel (fn)
      data(e).(fn{j}) = decodeValue (data(e).(fn{j}));
    endfor
  endfor

endfunction

## One saved value: an encoded array rebuilt, a structure searched in turn.
function v = decodeValue (v)

  if (! isstruct (v))
    return;
  endif
  fn = sort (fieldnames (v));
  if (isscalar (v) && isfield (v, 'LabelType')
      && strcmp (v.LabelType, 'categorical'))
    v = categoricalFrom (v.Codes, v.Categories, v.Ordinal, v.Protected, ...
                         v.Size);
  elseif (isscalar (v) && isfield (v, 'LabelType')
          && strcmp (v.LabelType, 'string'))
    v = stringFrom (v.Text, v.Missing, v.Size);
  elseif (isscalar (v) && isequal (fn, {'cats'; 'code'; 'isMissing'; ...
                                        'isOrdinal'; 'isProtected'}))
    codes = double (v.code);
    codes(v.isMissing) = NaN;
    v = categoricalFrom (codes, v.cats, v.isOrdinal, v.isProtected, ...
                         size (v.code));
  elseif (isscalar (v) && isequal (fn, {'isMissing'; 'strs'}))
    v = stringFrom (v.strs, v.isMissing, size (v.strs));
  else
    v = decodeLabels (v);
  endif

endfunction

## A categorical array from its codes, NaN marking an undefined element.
function C = categoricalFrom (codes, cats, ordinal, protected, sz)

  text = repmat ({''}, numel (codes), 1);
  ok = ! isnan (codes(:));
  text(ok) = cats(codes(ok));
  C = categorical (text, cats, 'Ordinal', ordinal, 'Protected', protected);
  C = reshape (C, sz);

endfunction

## A string array from its text and the mask of its missing elements.
function S = stringFrom (text, missing, sz)

  S = string (text(:));
  S(missing(:)) = string (NaN);
  S = reshape (S, sz);

endfunction
