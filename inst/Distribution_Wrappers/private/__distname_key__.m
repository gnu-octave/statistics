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
## @deftypefn {Private Function} {@var{key} =} __distname_key__ (@var{name})
##
## Fold a distribution name to the key its wrapper looks up.
##
## The key is @var{name} in lower case with spaces and hyphens removed, so
## that every spelling of a distribution reaches the same distribution in
## every wrapper.  @qcode{'Extreme Value'}, @qcode{'ExtremeValue'} and
## @qcode{'extreme-value'} all give @qcode{'extremevalue'}.  Both separators
## occur in names the package already ships -- @qcode{'Birnbaum-Saunders'},
## @qcode{'location-scale T'} -- which is why they are folded.
##
## The accepted set is deliberately a superset of MATLAB's, in two ways.
##
## MATLAB takes the spaced and the squashed spelling of a name but refuses
## the hyphenated one, so @qcode{'Birnbaum-Saunders'} and
## @qcode{'Log-Logistic'} are errors there.  Octave has accepted those since
## they are the display names its own tables carry, and withdrawing them to
## match MATLAB would break working code for no gain.  This is an extension,
## not a correction: MATLAB is not wrong to have a narrower set.
##
## MATLAB is, however, inconsistent with itself over
## @qcode{'tLocationScale'}: @code{makedist} accepts that name while
## @code{cdf} refuses it for the same distribution, which @code{cdf} calls
## @qcode{'location-scale T'}.  Measured against R2024a on 2026-08-07.  Every
## Octave wrapper accepts both names, so a name that builds a distribution
## also evaluates it.
##
## @end deftypefn

function key = __distname_key__ (name)

  key = tolower (name);
  key(key == " " | key == "-") = [];

endfunction
