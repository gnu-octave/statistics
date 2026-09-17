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

function post_install (desc)

  ## Octave before 11.2.0 rebuilds every doc-cache when the package is
  ## installed, one directory at a time and with only that directory on the
  ## path, so a class cannot resolve a superclass held in another directory
  ## and the install stops outright.  That walk takes back only the
  ## directories it puts on the path itself, so one added here stays on it
  ## for the whole of the walk.  This is what lets the model classes of
  ## Regression/ derive from PredictiveModel, which lives beside the learners
  ## it is mostly for.
  ##
  ## Measured on CI 2026-09-17: without this, both Octave 11.1.0 jobs failed
  ## on Regression/ and both 11.3.0 jobs passed, 11.2.0 having gained the
  ## check that keeps a valid shipped doc-cache rather than rebuilding it.
  ##
  ## The guard is deliberate: on 11.2.0 and later nothing is rebuilt, so the
  ## path is left exactly as the installation found it.  Where the version
  ## cannot be read, the path is added, which is the harmless way to be wrong.

  rebuilds = true;
  try
    rebuilds = compare_versions (version (), "11.2.0", "<");
  catch
    rebuilds = true;
  end_try_catch

  if (rebuilds)
    addpath (fullfile (desc.dir, "Supervised_Learning"));
  endif

endfunction
