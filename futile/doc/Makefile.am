EXTRA_DIST = \
	Futile_Logo.png \
	extra.css

CLEANFILES = *.aux *.toc *.log *.xml

#doc_DATA=$(doxygen-doc)

if BUILD_DEVDOC
insdevdocs = install-devdoc
uninsdevdocs = uninstall-devdoc
else
insdevdocs =
uninsdevdocs =
endif

install-data-hook: $(insdevdocs)

uninstall-local: $(uninsdevdocs)

install-devdoc: Doxygen/html
	test -z "$(docdir)" || /bin/mkdir -p "$(docdir)"
	cp -rp Doxygen/html '$(docdir)'

uninstall-devdoc:
	rm -rf '$(docdir)/html'

doxygen-doc: Doxygen/html

Doxygen/html: Doxygen.stamp web.stamp

Doxygen.stamp: $(builddir)/Doxyfile
	\rm -rf Doxygen;\
	doxygen $(builddir)/Doxyfile; \
	touch Doxygen.stamp

web.stamp: Doxygen.stamp
	touch web.stamp

clean-local:
	rm -rf Doxygen lcov* *.aux *.log  *.toc Doxygen *.stamp doxygen.warnings

PHONY = install-devdoc
