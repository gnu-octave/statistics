## Copyright (C) 2025 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## this program; if not, see <http://www.gnu.org/licens

function [varargout] = pairedArgs (optNames, dfValues, args)
  foundNames = [];
  nargs = numel (args);
  ## Search through all input arguments
  for ii = nargs-1:-1:1
    if (ischar (args{ii}))
      if (ismember (tolower (args{ii}), tolower (optNames)))
        idx = find (strcmpi (args{ii}, optNames));
        varargout{idx} = args{ii+1};
        foundNames = [foundNames, idx];
        args(ii:ii+1) = [];
      endif
    endif
  endfor
  ## Find optNames that were not in args and add defaults
  allNames = 1:numel (optNames);
  notfound = ! ismember (allNames, foundNames);
  defNames = optNames(notfound);
  for ii = 1:numel (defNames)
    idx = find (strcmpi (defNames{ii}, optNames));
    varargout{idx} = dfValues{idx};
  endfor
  ## Append remaining input arguments in varargout
  idx = numel (optNames) + 1;
  varargout{idx} = args(:);
endfunction
