#!/usr/bin/env bash
#
# Regenerate the per-directory `doc-cache` files shipped with the statistics
# package.  Octave builds these during `pkg install`; shipping up-to-date copies
# lets the installer skip that (slow) step and keeps `lookfor` current.
#
# The cache is generated per function directory (non-recursive), so there is one
# `doc-cache` per directory that holds functions.  Directories are discovered
# automatically (src, inst, and every inst/ subdirectory except private,
# datasets, demos, and tests), so new function directories are picked up without
# editing this script.
#
# Usage:
#   ./regen-doc-cache.sh                 # FULL: rebuild every directory (~minutes)
#   ./regen-doc-cache.sh <directory>     # rebuild one directory's cache only,
#                                        #   e.g.  ./regen-doc-cache.sh inst/dist_fun
#                                        #   (a bare name like 'dist_fun' is
#                                        #    resolved under inst/ as well)
#   ./regen-doc-cache.sh <name>          # update ONE function/class entry in its
#                                        #   directory's existing cache, e.g.
#                                        #   ./regen-doc-cache.sh anova1
#
# The single-name form is the fast path: it edits only the one column of the
# existing cache instead of rebuilding the whole directory, so refreshing a
# cache after a docstring change costs a second rather than several minutes.
# (If the name's directory has no cache yet, the whole directory is built.)
#
# Requirements:
#   * the oct-files must be built (run `make -C src` first) or src/doc-cache
#     will be incomplete;
#   * the `datatypes` dependency must be installed (it is loaded below).
#
#   OCTAVE=/path/to/octave ./regen-doc-cache.sh [arg]
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

OCTAVE="${OCTAVE:-octave}"

if ! ls src/*.oct >/dev/null 2>&1; then
  echo "warning: no built oct-files in src/ -- src/doc-cache will be incomplete." >&2
  echo "         run 'make -C src' first if the compiled functions changed." >&2
fi

DC_ARG="${1:-}" "$OCTAVE" --no-gui -q --eval "
  crash_dumps_octave_core (false);   ## no octave-workspace dump if interrupted
  pkg load datatypes;
  addpath (genpath (fullfile (pwd, 'inst')));

  arg = getenv ('DC_ARG');

  ## Decide the mode: no arg -> full; an existing folder -> that directory;
  ## anything else -> a single function/class name.
  mode = 'name';
  if (isempty (arg))
    mode = 'all';
  elseif (isfolder (arg))
    mode = 'dir';  d = arg;
  elseif (isfolder (fullfile ('inst', arg)))
    mode = 'dir';  d = fullfile ('inst', arg);
  endif

  if (strcmp (mode, 'all'))
    ## Discover the function directories that ship a doc-cache.
    skip = {'.', '..', 'private', 'datasets', 'demos', 'tests'};
    cand = {'src', 'inst'};
    s = dir ('inst');
    for k = 1:numel (s)
      if (s(k).isdir && ! any (strcmp (s(k).name, skip)))
        cand{end+1} = fullfile ('inst', s(k).name);
      endif
    endfor
    ## Keep only directories that actually contain functions.
    for k = 1:numel (cand)
      d = cand{k};
      if (! isfolder (d))
        continue;
      endif
      if (isempty (glob (fullfile (d, '*.m'))) && isempty (glob (fullfile (d, '*.oct'))))
        continue;
      endif
      t = tic;
      doc_cache_create (fullfile (d, 'doc-cache'), d);
      printf ('  regenerated %-26s (%.1fs)\n', fullfile (d, 'doc-cache'), toc (t));
      fflush (stdout);
    endfor
    printf ('doc-cache regeneration complete.\n');

  elseif (strcmp (mode, 'dir'))
    t = tic;
    doc_cache_create (fullfile (d, 'doc-cache'), d);
    printf ('  regenerated %-26s (%.1fs)\n', fullfile (d, 'doc-cache'), toc (t));

  else   ## single function/class name
    f = arg;
    fpath = which (f);
    if (isempty (fpath))
      error ('regen-doc-cache: function or class \"%s\" not found on the path.', f);
    endif
    d = fileparts (fpath);
    cachefile = fullfile (d, 'doc-cache');

    if (exist (cachefile, 'file') != 2)
      ## No cache in this directory yet -- build the whole thing.
      printf ('  no existing cache in %s; regenerating the directory ...\n', d);
      doc_cache_create (cachefile, d);
      printf ('  regenerated %s\n', cachefile);
    else
      ## Build the single cache entry exactly as doc_cache_create does per
      ## function: name / plain-text help / first help sentence.
      entry = {};
      if (! strncmp (f, '__', 2))
        [text, format] = get_help_text (f);
        status = 1;
        switch (lower (format))
          case 'plain text'
            status = 0;
          case 'texinfo'
            [text, status] = __makeinfo__ (text, 'plain text');
        endswitch
        if (status == 0 && ! isempty (text))
          entry = {f; text; get_first_help_sentence(f)};
        else
          warning ('regen-doc-cache: unusable help text in %s', fpath);
        endif
      endif

      ## Splice the entry into the existing cache (native Octave -text data).
      tmp = load (cachefile);
      cache = tmp.cache;
      idx = find (strcmp (cache(1,:), f));
      if (isempty (entry))
        ## Internal (__name__) or unusable -> drop any stale column.
        if (! isempty (idx))
          cache(:, idx) = [];
          printf ('  removed stale entry \"%s\" from %s\n', f, cachefile);
        else
          printf ('  \"%s\" is internal/unusable; nothing to cache.\n', f);
        endif
      elseif (isempty (idx))
        cache(:, end+1) = entry;
        printf ('  added \"%s\" to %s\n', f, cachefile);
      else
        cache(:, idx) = entry;
        printf ('  updated \"%s\" in %s\n', f, cachefile);
      endif

      ## Save with the same header doc_cache_create writes.
      save_header_format_string (['# doc-cache created by Octave ' OCTAVE_VERSION]);
      save ('-text', cachefile, 'cache');
    endif
  endif
"
