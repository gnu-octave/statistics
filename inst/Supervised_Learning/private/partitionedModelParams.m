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
## @deftypefn  {Private Function} {@var{MP} =} partitionedModelParams (@var{Mdl}, @var{kfold}, @var{method}, @var{type})
## @deftypefnx {Private Function} {@var{MP} =} partitionedModelParams (@dots{}, @var{learner})
##
## The @qcode{ModelParameters} of a cross-validated model.
##
## @var{Mdl} is the learner the folds were fitted from, @var{kfold} the number
## of folds, @var{method} the name the partitioned class reports and
## @var{type} either @qcode{'classification'} or @qcode{'regression'}.
## @var{learner} names the backing where the class publishes one.
##
## The learner's own parameters are carried through, so that a cross-validated
## model still says what its folds were fitted with.  Its @qcode{Version},
## @qcode{Method} and @qcode{Type} tags describe the learner rather than the
## partition and are reissued for the class doing the reporting, else a
## cross-validated SVM would answer @qcode{'SVM'} where it is a partitioned
## model.
##
## @end deftypefn

function MP = partitionedModelParams (Mdl, kfold, method, type, learner)

  MP = Mdl.ModelParameters;
  if (! isstruct (MP))
    MP = struct ();
  endif

  ## The tags belong to the reporting class, not to the learner.
  tags = {'Version', 'Method', 'Type'};
  for i = 1:numel (tags)
    if (isfield (MP, tags{i}))
      MP = rmfield (MP, tags{i});
    endif
  endfor

  if (nargin > 4 && ! isempty (learner))
    MP.LearnerTemplates = learner;
  endif
  MP.NLearn = kfold;
  MP.Version = 1;
  MP.Method = method;
  MP.Type = type;

endfunction
