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
## @deftypefn {Private Function} {[@var{T}, @var{errmsg}] =} templateStruct (@var{method}, @var{args})
##
## Build a learner template from a list of name-value arguments.
##
## @var{method} names the learner and @var{args} is the caller's
## @code{varargin}.  @var{T} is a structure carrying @qcode{Method},
## @qcode{Type} and one field per name given, the value stored as it stands.
##
## The names themselves are not checked against the learner's own list.  The
## learner owns that list and validates it when the template is consumed, so
## restating it here would give two lists to keep in step and a template that
## silently accepts what its learner refuses once they drift apart.  What is
## checked here is only what must hold for a structure to be built at all.
##
## @var{errmsg} is empty when the template was built, and otherwise the body
## of a message the caller raises under its own name.
##
## @seealso{templateSVM, templateTree, templateKNN}
## @end deftypefn

function [T, errmsg] = templateStruct (method, args)

  T = struct ('Method', method, 'Type', 'classification');
  errmsg = '';

  if (mod (numel (args), 2) != 0)
    errmsg = "name-value arguments must be in pairs.";
    return;
  endif

  for i = 1:2:numel (args)
    name = args{i};
    if (! (ischar (name) && isrow (name)))
      errmsg = "parameter names must be character vectors.";
      return;
    endif
    ## A name that is not an identifier cannot be a field, and no learner
    ## takes one, so it is refused here rather than from inside setfield.
    if (! isvarname (name))
      errmsg = strcat ("'", name, "' is not a valid parameter name.");
      return;
    endif
    T.(name) = args{i+1};
  endfor

endfunction
