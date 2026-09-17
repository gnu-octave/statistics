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
## @deftypefn {Private Function} {@var{v} =} encodeLabels (@var{v})
##
## A categorical or string array as a structure @code{save} can write.
##
## @code{save} cannot write a classdef object, so a model saving its class
## names or its response in either type writes this structure in their place,
## and @code{decodeLabels} rebuilds the array when the model is loaded.  A
## categorical array keeps its codes, categories, order and protection, and a
## string array its text and which of its elements are missing.  A structure
## array has every field of every element encoded in turn, which reaches class
## names a model keeps inside a structure.  Any other value is returned as it
## is.
##
## @seealso{decodeLabels}
## @end deftypefn

function v = encodeLabels (v)

  if (isa (v, 'categorical'))
    v = struct ('LabelType', 'categorical', 'Size', size (v), ...
                'Codes', double (v), 'Categories', {categories(v)}, ...
                'Ordinal', isordinal (v), 'Protected', isprotected (v));
  elseif (isa (v, 'string'))
    v = struct ('LabelType', 'string', 'Size', size (v), ...
                'Text', {cellstr(v)}, 'Missing', ismissing (v));
  elseif (isstruct (v))
    fn = fieldnames (v);
    for e = 1:numel (v)
      for j = 1:numel (fn)
        v(e).(fn{j}) = encodeLabels (v(e).(fn{j}));
      endfor
    endfor
  endif

endfunction
