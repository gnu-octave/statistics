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
## @deftypefn {Private Function} {[@var{C}, @var{errmsg}] =} costMatrix (@var{val}, @var{ClassNames})
## Resolve a misclassification cost to a plain double matrix in the model's
## own class order.
##
## @var{val} is a numeric matrix, taken in the order the model holds its
## classes, or a struct with the
## fields @qcode{ClassNames} and @qcode{ClassificationCosts}, whose matrix is
## written in the order that struct names and is permuted into the order of
## @var{ClassNames}.  The struct form is how a caller supplies costs without
## having to know which order the model sorted its classes into.
##
## A cost must be floating point, not sparse, not complex, square with one row
## and column per class, non-negative, zero down the diagonal, and free of
## @qcode{NaN} and @qcode{Inf}.  A @code{single} is widened to @code{double}
## rather than refused.  All of it measured against MATLAB R2024a.
##
## @var{errmsg} is the body of the message the caller should raise, or empty
## when there is nothing wrong; the caller emits it under its own
## @code{class.method} name, as the package's shared validation helpers do.
## An empty @var{val} is returned unchanged, the caller owning what an empty
## cost resets to.
## @end deftypefn

function [C, errmsg] = costMatrix (val, ClassNames)

  C = val;
  errmsg = '';

  if (isempty (val))
    return;
  endif

  K = classCount (ClassNames);

  ## A struct names the order its matrix is written in.  Permuting it into the
  ## model's order is not a convenience: the two orders genuinely differ, and
  ## MATLAB permutes entry by entry, so this is measured behaviour and not a
  ## reading of the matrix as given.
  if (isstruct (val))
    if (! (isfield (val, 'ClassNames') && isfield (val, 'ClassificationCosts')))
      errmsg = strcat ("'Cost' given as a struct must have the fields", ...
                       " 'ClassNames' and 'ClassificationCosts'.");
      return;
    endif
    given = val.ClassNames;
    C = val.ClassificationCosts;
    n = classCount (given);
    if (! isequal (size (C), [n, n]))
      errmsg = strcat ("'Cost' given as a struct must have one row and", ...
                       " one column of 'ClassificationCosts' per name in", ...
                       " 'ClassNames'.");
      return;
    endif
    ## Matching is by name, so the struct must name the classes in the type
    ## the model holds them in.  A type it cannot compare reaches ismember
    ## and raises rather than coming back as a message, so it is caught here
    ## and named for what it is; labelIndices speaks of a response, which a
    ## cost is not.
    perm = [];
    try
      perm = labelIndices (given, ClassNames);
    catch
      perm = [];
    end_try_catch
    if (isempty (perm))
      errmsg = strcat ("'Cost' given as a struct must name the classes in", ...
                       " the same type as the model holds them.");
      return;
    endif
    ## A class of the model that the struct does not name has no cost to
    ## permute, which MATLAB refuses rather than defaulting.
    missing = find (perm == 0, 1);
    if (! isempty (missing))
      errmsg = sprintf (strcat ("'Cost' given as a struct must name every", ...
                                " class; there is no cost for class %d."), ...
                        missing);
      return;
    endif
    C = C(perm, perm);
  endif

  ## Char, cell and logical are none of them costs.  An integer type is a
  ## number but not a cost either: MATLAB refuses it in its own right, and
  ## says so differently, so the two are kept apart here as well.
  if (! isnumeric (C))
    errmsg = strcat ("'Cost' must be a numeric matrix, or a struct with", ...
                     " the fields 'ClassNames' and 'ClassificationCosts'.");
    return;
  endif
  if (isinteger (C))
    errmsg = "'Cost' must be a floating point matrix.";
    return;
  endif
  if (issparse (C))
    errmsg = "'Cost' must not be sparse.";
    return;
  endif
  ## Complex means a nonzero imaginary part, not the complex attribute:
  ## complex (C, 0) is accepted and stored as a double, measured on R2024a.
  ## Writing this as ! isreal (C) would refuse a value MATLAB takes.
  if (any (imag (C(:)) != 0))
    errmsg = "'Cost' must not be complex.";
    return;
  endif

  ## Single is widened rather than refused, and the imaginary part of an
  ## accepted complex value is dropped with it.
  C = double (real (C));

  if (! isequal (size (C), [K, K]))
    errmsg = strcat ("the number of rows and columns in 'Cost' must", ...
                     " correspond to selected classes in Y.");
    return;
  endif

  if (any (C(:) < 0))
    errmsg = "'Cost' must not contain negative values.";
  elseif (any (diag (C) != 0))
    errmsg = "'Cost' must have zeros on its diagonal.";
  elseif (any (! isfinite (C(:))))
    errmsg = "'Cost' must not contain NaN or Inf values.";
  endif

endfunction
