## Copyright 2015-2016 Carnë Draug
## Copyright 2015-2016 Oliver Heimlich
##
## Copying and distribution of this file, with or without modification,
## are permitted in any medium without royalty provided the copyright
## notice and this notice are preserved.  This file is offered as-is,
## without any warranty.

PACKAGE = $(shell grep "^Name: " DESCRIPTION | cut -f2 -d" ")
VERSION = $(shell grep "^Version: " DESCRIPTION | cut -f2 -d" ")

RELEASE_DIR     = $(PACKAGE)-$(VERSION)
RELEASE_TARBALL = $(PACKAGE)-$(VERSION).tar.gz
HTML_DIR        = $(PACKAGE)-html
HTML_TARBALL    = $(PACKAGE)-html.tar.gz

M_SOURCES   = $(wildcard inst/*.m)
PKG_ADD     = $(shell grep -Pho '(?<=// PKG_ADD: ).*' $(M_SOURCES))

OCTAVE ?= octave

.PHONY: help dist html release install check run clean

help:
	@echo "Targets:"
	@echo "   dist    - Create $(RELEASE_TARBALL) for release"
	@echo "   html    - Create $(HTML_TARBALL) for release"
	@echo "   release - Create both of the above and show md5sums"
	@echo
	@echo "   install - Install the package in GNU Octave"
	@echo "   check   - Execute package tests (w/o install)"
	@echo "   run     - Run Octave with development in PATH (no install)"
	@echo
	@echo "   clean   - Remove releases, html documentation, and oct files"

$(RELEASE_DIR): .hg/dirstate
	@echo "Creating package version $(VERSION) release ..."
	-rm -rf "$@"
	hg archive --exclude ".hg*" --exclude "Makefile" --type files "$@"
	chmod -R a+rX,u+w,go-w "$@"

$(RELEASE_TARBALL): $(RELEASE_DIR)
	tar cf - --posix "$<" | gzip -9n > "$@"

$(HTML_DIR): install
	@echo "Generating HTML documentation. This may take a while ..."
	-rm -rf "$@"
	$(OCTAVE) --silent \
	  --eval "pkg load generate_html; " \
	  --eval "pkg load $(PACKAGE);" \
	  --eval 'generate_package_html ("${PACKAGE}", "$@", "octave-forge");'
	chmod -R a+rX,u+w,go-w "$@"

$(HTML_TARBALL): $(HTML_DIR)
	tar cf - --posix "$<" | gzip -9n > "$@"

dist: $(RELEASE_TARBALL)
html: $(HTML_TARBALL)

release: dist html
	md5sum $(RELEASE_TARBALL) $(HTML_TARBALL)
	@echo "Upload @ https://sourceforge.net/p/octave/package-releases/new/"
	@echo 'Execute: hg tag "release-${VERSION}"'

install: $(RELEASE_TARBALL)
	@echo "Installing package locally ..."
	$(OCTAVE) --silent --eval 'pkg ("install", "${RELEASE_TARBALL}")'

check:
	$(OCTAVE) --silent \
	  --eval 'addpath (fullfile ([pwd filesep "inst"]));' \
	  --eval '${PKG_ADD}' \
	  --eval 'runtests ("inst");'

run:
	$(OCTAVE) --no-gui --silent --persist \
	  --eval 'addpath (fullfile ([pwd filesep "inst"]));' \
	  --eval '${PKG_ADD}'

clean:
	rm -rf $(RELEASE_DIR) $(RELEASE_TARBALL) $(HTML_TARBALL) $(HTML_DIR)

