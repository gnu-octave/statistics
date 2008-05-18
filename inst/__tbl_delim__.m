## Copyright (C) 2008 Bill Denney
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {[@var{d}, @var{err}] = } tblread (@var{d})
## Return the delimiter for tblread or tblwrite.
##
## The delimeter, @var{d} may be any single character or
## @itemize
## @item "space" " " (default)
## @item "tab" "\t"
## @item "comma" ","
## @item "semi" ";"
## @item "bar" "|"
## @end itemize
##
## @var{err} will be empty if there is no error, and @var{d} will be NaN
## if there is an error.  You MUST check the value of @var{err}.
## @seealso{tblread, tblwrite}
## @end deftypefn

## Author: Bill Denney <bill@denney.ws>

function [d, err] = __tbl_delim__ (d)

  ## Check arguments
  if nargin != 1
    print_usage ();
  endif

  err = "";
  ## Format the delimiter
  if ischar (d)
    ## allow for escape characters
    d = sprintf (d);
    if numel (d) > 1
      ## allow the word forms
      s.space = " ";
      s.tab = "\t";
      s.comma = ",";
      s.semi = ";";
      s.bar = "|";
      if ! ismember (d, fieldnames (s))
        err = ["tblread: delimiter must be either a single " ...
               "character or one of\n" ...
               sprintf("%s, ", fieldnames (s))(1:end-2)];
        d = NaN;
      else
        d = s.(d);
      endif
    endif
  else
    err = "delimiter must be a character";
    d = NaN;
  endif
  if isempty (d)
    err = "the delimiter may not be empty";
    d = NaN;
  endif

endfunction

## Tests
## The defaults
%!test
%! [d err] = __tbl_delim__ (" ");
%! assert (d, " ");
%! assert (err, "");
## Named delimiters
%!test
%! [d err] = __tbl_delim__ ("space");
%! assert (d, " ");
%! assert (err, "");
%!test
%! [d err] = __tbl_delim__ ("tab");
%! assert (d, sprintf ("\t"));
%! assert (err, "");
%!test
%! [d err] = __tbl_delim__ ("comma");
%! assert (d, ",");
%! assert (err, "");
%!test
%! [d err] = __tbl_delim__ ("semi");
%! assert (d, ";");
%! assert (err, "");
%!test
%! [d err] = __tbl_delim__ ("bar");
%! assert (d, "|");
%! assert (err, "");
## An arbitrary character
%!test
%! [d err] = __tbl_delim__ ("x");
%! assert (d, "x");
%! assert (err, "");
## An arbitrary escape string
%!test
%! [d err] = __tbl_delim__ ('\r');
%! assert (d, sprintf ('\r'))
%! assert (err, "");
## Errors
%!test
%! [d err] = __tbl_delim__ ("bars");
%! assert (isnan (d));
%! assert (! isempty (err));
%!test
%! [d err] = __tbl_delim__ ("");
%! assert (isnan (d));
%! assert (! isempty (err));
%!test
%! [d err] = __tbl_delim__ (5);
%! assert (isnan (d));
%! assert (! isempty (err));
%!test
%! [d err] = __tbl_delim__ ({"."});
%! assert (isnan (d));
%! assert (! isempty (err));

