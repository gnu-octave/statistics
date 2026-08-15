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
## @deftypefn  {statistics} {@var{value} =} statget (@var{options}, @var{name})
## @deftypefnx {statistics} {@var{value} =} statget (@var{options}, @var{name}, @var{default})
##
## Read one option out of a statistics options structure.
##
## @code{@var{value} = statget (@var{options}, @var{name})} returns the value
## the option @var{name} carries in @var{options}, or @code{[]} when that
## option is unset.  @var{options} is a structure as built by @code{statset},
## although any structure is accepted.
##
## @code{@var{value} = statget (@var{options}, @var{name}, @var{default})}
## returns @var{default} instead whenever the option is unset, which is the
## form a calling function uses to fall back on its own default.  Note that
## @var{default} is returned when the option is @emph{empty}, not only when it
## is absent, since an empty option is precisely how @code{statset} spells
## @qcode{"unset"}.
##
## @var{name} is matched case-insensitively, and may be abbreviated to any
## leading portion that singles out one option: @code{statget (@var{options},
## "MaxI")} reads @qcode{"MaxIter"}.  An abbreviation matching more than one
## option raises, rather than choosing between them; an exact match is taken
## as exact even where it is also a prefix of a longer name.
##
## @seealso{statset}
## @end deftypefn

function value = statget (options, name, default)

  if (nargin < 2)
    print_usage ();
  endif

  if (nargin < 3)
    default = [];
  endif

  if (! (isstruct (options) && isscalar (options)))
    error (strcat ("statget: OPTIONS must be a scalar structure, as", ...
                   " built by statset."));
  endif

  if (! (ischar (name) || (isa (name, 'string') && isscalar (name))))
    error ("statget: NAME must be a character vector or a string scalar.");
  endif
  name = char (name);

  fn = fieldnames (options);

  ## An exact match wins outright, so that a name which is also the prefix of
  ## a longer one is never treated as ambiguous.
  idx = find (strcmpi (name, fn));

  if (isempty (idx))
    idx = find (strncmpi (name, fn, numel (name)));
    if (numel (idx) > 1)
      error (strcat ("statget: '%s' matches more than one option name:", ...
                     " %s."), name, strjoin (fn(idx)', ", "));
    endif
  endif

  if (isempty (idx))
    error ("statget: unrecognized option name '%s'.", name);
  endif

  value = options.(fn{idx});

  if (isempty (value))
    value = default;
  endif

endfunction

%!demo
%! ## Read an option, falling back on a default when it is unset
%! options = statset ('nlinfit');
%! maxiter = statget (options, 'MaxIter')
%! tolbnd = statget (options, 'TolBnd', 1e-6)

## Reading a set option
%!test
%! assert_equal (statget (statset ('factoran'), 'TolX'), 1e-8);

%!test
%! assert_equal (statget (statset ('nlinfit'), 'MaxIter'), 200);

%!test
%! assert_equal (statget (statset ('nlinfit'), 'WgtFun'), 'bisquare');

## An unset option reads as empty
%!test
%! assert_equal (statget (statset ('factoran'), 'TolBnd'), []);

%!test
%! assert_equal (statget (statset (), 'MaxIter'), []);

## The default is returned only when the option is unset
%!test
%! assert_equal (statget (statset ('factoran'), 'TolBnd', 42), 42);

%!test
%! assert_equal (statget (statset ('factoran'), 'TolX', 42), 1e-8);

## The default may be of any type
%!test
%! assert_equal (statget (statset (), 'Display', 'final'), 'final');

%!test
%! assert_equal (statget (statset (), 'OutputFcn', {@sin}), {@sin});

## The name is matched case-insensitively
%!test
%! assert_equal (statget (statset ('factoran'), 'tolx'), 1e-8);

%!test
%! assert_equal (statget (statset ('factoran'), 'MAXITER'), 100);

## The name may be abbreviated when the abbreviation is unique
%!test
%! assert_equal (statget (statset ('factoran'), 'MaxI'), 100);

%!test
%! assert_equal (statget (statset ('factoran'), 'Displ'), 'off');

%!test
%! assert_equal (statget (statset ('nlinfit'), 'Deriv'), 6.0554544523933429e-06, 1e-20);

## An exact match is preferred over the longer names it prefixes
%!test
%! assert_equal (statget (statset ('nnmf'), 'TolX'), 1e-4);

%!test
%! assert_equal (statget (statset ('fitnlm'), 'Robust'), 'off');

## Any structure is accepted, not only one built by statset
%!test
%! assert_equal (statget (struct ('MaxIter', 7), 'MaxIter'), 7);

%!test
%! assert_equal (statget (struct ('MaxIter', 7), 'MaxI'), 7);

## Error conditions
%!error<Invalid call to statget> statget ()

%!error<Invalid call to statget> statget (statset ())

%!error<statget: OPTIONS must be a scalar structure, as built by statset.> ...
%! statget (1, 'MaxIter')

%!error<statget: OPTIONS must be a scalar structure, as built by statset.> ...
%! statget (struct ('MaxIter', {1, 2}), 'MaxIter')

%!error<statget: NAME must be a character vector or a string scalar.> ...
%! statget (statset (), 5)

%!error<statget: 'Tol' matches more than one option name: TolBnd, TolFun, TolTypeFun, TolX, TolTypeX.> ...
%! statget (statset (), 'Tol')

%!error<statget: 'TolT' matches more than one option name: TolTypeFun, TolTypeX.> ...
%! statget (statset (), 'TolT')

%!error<statget: unrecognized option name 'NoSuchOption'.> ...
%! statget (statset (), 'NoSuchOption')

%!error<statget: unrecognized option name 'TolX'.> ...
%! statget (struct ('MaxIter', 7), 'TolX')
