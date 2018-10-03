## Copyright (C) 2018 Olaf Till <i7tiol@t-online.de>
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
## @deftypefn {Function File} {} conditional_installation (@var{srcdir}, @var{destdir}, @var{crosscompiling}, @var{indexfilename})
## Undocumented internal function.
## @end deftypefn

function conditional_installation (srcdir, destdir, crosscompiling, indexfilename)

  ## Cross-compilation: some version of Octave must be callable at the
  ## build system, and this file here has to be configured as follows:
  ##
  ## Set this to true if you have considered the next step.
  configured = false;
  ##
  ## If core Octave called at the build system has the same statistics
  ## functions as the core Octave you build for, leave the following
  ## unchanged. _Only_ if not, uncomment the following. If you
  ## uncomment the following, the sourced file has to be customized.
  ##
  ## source (fullfile (fileparts (mfilename ("fullpath")),
  ##                   "functions_to_install")); # defines variable
  ##                                             # 'install_functions'

  if exist ("isfolder") == 0
    isfolder = @(n) isdir(n);
  endif

  installed_functions = {};
  subdirs = {"base", "distributions", "models", "tests"};
  
  if (crosscompiling && ! configured)

    error ("for cross-compilation you have to customize 'conditional_installation.m'");

  endif

  ## assert destination subdirs and install private functions
  for tsubdir = subdirs

    subdir = tsubdir{:};

    assert_dir (fullfile (destdir, subdir));

    if (isfolder (private_dir = fullfile (srcdir, subdir, "private"))
        && ! ([status, msg] = ...
              copyfile (private_dir, fullfile (destdir, subdir))))

      error ("could not copy %s: %s", private_dir, msg);

    endif

  endfor

  if (crosscompiling && exist ("install_functions") == 1)

    for tfcn = install_functions

      fcn = tfcn{:};
  
      installed_functions{end+1} = fcn;

      if (! ([status, msg] = ...
             copyfile (fullfile (srcdir, fcn), fullfile (destdir, fcn))))

        error ("could not copy %s: %s", fcn, msg);

      endif

    endfor

  else

    if (crosscompiling)

      warning ("core Octave statistics functions are assumed to be the same at build- and host-systems; if this is not the case, customize the statistics package according to the comment in 'conditional_installation.m'");

    endif

    for tsubdir = subdirs

      subdir = tsubdir{:};

      files = glob (fullfile (srcdir, subdir, "*.m"));

      for tfile = files'

        file = tfile{:};

        [~, fcn] = fileparts (file);

        ## FIXME: the following needs linked code, similar to that
        ## in Octaves load_path::package_info::add_to_fcn_map
        if ((flag = exist (fcn)) == 1 || flag == 7)

          error ("can't check presence of '%s' function in core Octave since a variable or a directory in the search path exists with that name",
                 fcn);

        endif

        if (! flag)

          installed_functions{end+1} = fcn;

          if (! ([status, msg] = ...
                 copyfile (file, fullfile (destdir, subdir, [fcn, ".m"]))))

            error ("could not copy %s: %s", fcn, msg);

          endif

        endif

      endfor

    endfor

  endif

  update_index_file (indexfilename, fullfile (srcdir, "INDEX.in"), installed_functions);
endfunction

function assert_dir (directory)

  if exist ("isfolder") == 0
    isfolder = @(n) isdir(n);
  endif

  if (! isfolder (directory))

    if (! ([succ, msg] = mkdir (directory)))

      error ("Could not create '%s': %s", directory, msg);

    endif

  endif

endfunction

