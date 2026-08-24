## Copyright (C) 2026 Andreas Bertsatos <abertsatos@biol.uoa.gr>
## Copyright (C) 2026 Avanish Salunke <avanishsalunke16@gmail.com>
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
## @deftypefn {Private Function} {@var{t} =} robusttune (@var{wfun})
##
## Default tuning constant, giving 95% efficiency, for a robust weight function.
##
## @var{wfun} is one of the names @code{robustwfun} accepts, or a function
## handle, which carries no default and takes 1.
##
## @seealso{robustfit, nlinfit, robustwfun}
## @end deftypefn

function t = robusttune (wfun)

  if (is_function_handle (wfun))
    t = 1;
    return;
  endif
  switch (wfun)
    case "andrews";      t = 1.339;
    case "bisquare";     t = 4.685;
    case "cauchy";       t = 2.385;
    case "fair";         t = 1.400;
    case "huber";        t = 1.345;
    case "logistic";     t = 1.205;
    case "ols";          t = 1;
    case "talwar";       t = 2.795;
    case "welsch";       t = 2.985;
  endswitch

endfunction
