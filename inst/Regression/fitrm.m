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
## @deftypefn  {statistics} {@var{rm} =} fitrm (@var{t}, @var{modelspec})
## @deftypefnx {statistics} {@var{rm} =} fitrm (@var{t}, @var{modelspec}, @var{name}, @var{value})
##
## Fit a repeated measures model.
##
## @code{@var{rm} = fitrm (@var{t}, @var{modelspec})} fits the repeated
## measures held in the table @var{t}, one row per subject, on the
## between-subject model @var{modelspec}, and returns a
## @code{RepeatedMeasuresModel} object.  @var{modelspec} is a character vector
## or a string scalar of the form @qcode{'@var{responses} ~ @var{terms}'}: the
## responses a range of the table's variables such as @qcode{'y1-y6'}, a comma
## list such as @qcode{'y1,y2,y3'}, or both, and the terms a Wilkinson formula
## over the other variables, such as @qcode{'species'} or @qcode{'g*x'}, or
## @qcode{'1'} for an intercept alone.
##
## @code{@var{rm} = fitrm (@dots{}, @var{name}, @var{value})} takes the
## within-subject design and model as the name-value arguments
## @qcode{'WithinDesign'} and @qcode{'WithinModel'}; see
## @code{RepeatedMeasuresModel} for both.
##
## @code{ranova}, @code{epsilon} and @code{mauchly} test the fitted model.
##
## @seealso{RepeatedMeasuresModel, anova2, manova1}
## @end deftypefn

function rm = fitrm (t, modelspec, varargin)

  if (nargin < 2)
    error ("fitrm: too few input arguments.");
  endif
  rm = RepeatedMeasuresModel (t, modelspec, varargin{:});

endfunction

%!demo
%! ## Four measurements on each iris, and whether they differ by species
%! load fisheriris
%! t = table (species, meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!            'VariableNames', {'species', 'meas1', 'meas2', 'meas3', 'meas4'});
%! Meas = table ([1 2 3 4]', 'VariableNames', {'Measurements'});
%! rm = fitrm (t, 'meas1-meas4 ~ species', 'WithinDesign', Meas)
%! ranova (rm)
%! mauchly (rm)
%! epsilon (rm)

%!test
%! load fisheriris
%! t = table (species, meas(:,1), meas(:,2), meas(:,3), meas(:,4), ...
%!            'VariableNames', {'species', 'meas1', 'meas2', 'meas3', 'meas4'});
%! rm = fitrm (t, 'meas1-meas4 ~ species');
%! assert_equal (class (rm), 'RepeatedMeasuresModel');

%!error<fitrm: too few input arguments.> fitrm ()
%!error<fitrm: too few input arguments.> fitrm (1)
