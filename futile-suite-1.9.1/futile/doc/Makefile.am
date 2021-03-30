SPHINX_RST = conf.py \
	_static \
	f_prec.rst \
	index.rst \
	inputfile.rst \
	pyfutile.rst \
	yaml_parse.rst \
	yaml_str.rst

EXTRA_DIST = \
	Futile_Logo.png \
	extra.css \
	$(SPHINX_RST)

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

install-devdoc: sphinx
	test -z "$(docdir)" || /bin/mkdir -p "$(docdir)"
	cp -rp sphinx '$(docdir)'

devdoc: sphinx-doc

sphinx-doc: sphinx

sphinx: $(SPHINX_RST)
	sphinx-build $(srcdir) "sphinx"

uninstall-devdoc:
	rm -rf '$(docdir)/html'

doxygen-doc: Doxygen/html

Doxygen/html: Doxygen.stamp web.stamp

Doxygen.stamp: $(top_builddir)/config/Doxyfile
	\rm -rf Doxygen;\
	doxygen $(top_builddir)/config/Doxyfile; \
	touch Doxygen.stamp

web.stamp: Doxygen.stamp
	touch web.stamp

ford-doc: FORD

FORD: FORD.stamp

FORD.stamp: $(top_builddir)/config/ford.md
	\rm -rf FORD
	ford --src_dir $(top_srcdir)/flib \
		--include $(top_srcdir)/flib --include $(top_builddir)/flib \
		--output_dir FORD $(top_srcdir)/config/ford.md

clean-local:
	rm -rf FORD Doxygen lcov* *.aux *.log  *.toc Doxygen *.stamp doxygen.warnings

PHONY = install-devdoc devdoc
