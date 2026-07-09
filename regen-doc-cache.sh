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
# Requirements:
#   * the oct-files must be built (run `make -C src` first) or src/doc-cache
#     will be incomplete;
#   * the `datatypes` dependency must be installed (it is loaded below).
#
# Usage:  ./regen-doc-cache.sh          (run from anywhere)
#         OCTAVE=/path/to/octave ./regen-doc-cache.sh
#
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

OCTAVE="${OCTAVE:-octave}"

if ! ls src/*.oct >/dev/null 2>&1; then
  echo "warning: no built oct-files in src/ -- src/doc-cache will be incomplete." >&2
  echo "         run 'make -C src' first if the compiled functions changed." >&2
fi

"$OCTAVE" --no-gui -q --eval "
  pkg load datatypes;
  addpath (genpath (fullfile (pwd, 'inst')));

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
"
