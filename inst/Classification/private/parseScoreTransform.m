## Copyright (C) 2024 Andreas Bertsatos <abertsatos@biol.uoa.gr>
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
## @deftypefn {Private Function} parseScoreTransform (@var{pd}, @var{classname})
##
## Parse ScoreTransform for Classification objects.
##
## @end deftypefn

function f = parseScoreTransform (ScoreTransform, classname)

  stList = {"doublelogit", "invlogit", "ismax", "logit", "none", ...
            "identity", "sign", "symmetric", "symmetricismax", ...
            "symmetriclogit"};
  if (! (ischar (ScoreTransform) ||
         strcmp (class (ScoreTransform), "function_handle")))
    error (strcat (["%s: 'ScoreTransform' must be a character"], ...
                   [" vector or a function handle."]), classname);
  endif
  if (! ismember (ScoreTransform, stList))
    error ("%s: unrecognized 'ScoreTransform' function.", classname);
  endif
  ## Handle ScoreTransform here
  if (is_function_handle (ScoreTransform))
    m = eye (5);
    if (! isequal (size (m), size (ScoreTransform (m))))
      error (strcat (["%s: function handle for 'ScoreTransform' must"], ...
                     [" return the same size as its input."]), classname);
    endif
    f = ScoreTransform;
  else
    if (strcmpi ("doublelogit", ScoreTransform))
      f = @(x) 1 ./ (1 + exp .^ (-2 * x));
    elseif (strcmpi ("invlogit", ScoreTransform))
      f = @(x) log (x ./ (1 - x));
    elseif (strcmpi ("ismax", ScoreTransform))
      f = eval (sprintf ("@(x) ismax (x)"));
    elseif (strcmpi ("logit", ScoreTransform))
      f = @(x) 1 ./ (1 + exp .^ (-x));
    elseif (strcmpi ("identity", ScoreTransform))
      f = 'none';
    elseif (strcmpi ("sign", ScoreTransform))
      f = @(x) sign (x);
    elseif (strcmpi ("symmetric", ScoreTransform))
      f = @(x) 2 * x - 1;
    elseif (strcmpi ("symmetricismax", ScoreTransform))
      f = eval (sprintf ("@(x) symmetricismax (x)"));
    elseif (strcmpi ("symmetriclogit", ScoreTransform))
      f = @(x) 2 ./ (1 + exp .^ (-x)) - 1;
    endif
  endif

endfunction

## Helper functions for ScoreTransform
function out = ismax (score)
  out = score;
  out(score == max (score)) = 1;
  out(score != max (score)) = 0;
endfunction

function out = symmetricismax (score)
  out = score;
  out(score == max (score)) = 1;
  out(score != max (score)) = -1;
endfunction
