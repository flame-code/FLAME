## Process this file with automake to produce Makefile.in

EXTRA_DIST = fldiff_yaml.py.in report.py.in README $(pytest_DATA)

SUBDIRS = \
	flib \
	fft \
	wrappers \
	c-bindings

#Installation of some python file (use $pythondir to access into).
python_PYTHON = yaml_hl.py fldiff.py f_regtest.py UniParse.py
nodist_python_PYTHON = fldiff_yaml.py report.py
dist_python_DATA = yaml_hl.cfg yaml_hl_remarks.cfg

#Used by jupyter notebooks
#pytest_DATA = test_input_variables_definition.yaml test-input.yaml
pytestdir = .

#remove some files in order to create it with the right path (possible to use nodist option)
dist-hook:
	rm -f report.py fldiff.py

if USE_OCL
oclrun_message = oclrun
else
oclrun_message =
endif

if USE_MPI
mpirun_message = mpirun
else
mpirun_message =
endif

BUILT_SOURCES = fldiff_yaml.py report.py

fldiff_yaml.py: fldiff_yaml.py.in
	@$(SED)							\
		-e s!\@prefix\@!$(prefix)!			\
		-e s!\@exec_prefix\@!$(exec_prefix)!		\
		-e s!\@pyexecdir\@!$(pyexecdir)!		\
		-e s!\@AX_PYYAML_PATH\@!$(AX_PYYAML_PATH)!	\
		-e s!\@abs_srcdir\@!$(abs_srcdir)!		\
		< $< > $@

report.py: report.py.in
	@$(SED)							\
		-e s!\@prefix\@!$(prefix)!			\
		-e s!\@exec_prefix\@!$(exec_prefix)!		\
		-e s!\@pyexecdir\@!$(pyexecdir)!		\
		-e s!\@AX_PYYAML_PATH\@!$(AX_PYYAML_PATH)!	\
		-e s!\@abs_srcdir\@!$(abs_srcdir)!		\
		< $< > $@

head_message: all-recursive
	@echo "========================================================"
	@echo "To run the tests, type 'make check' (might take a while)"

foot_message: $(mpirun_message) $(oclrun_message) head_message
	@echo "========================================================"

mpirun: head_message
	@echo "     Use the environment variable run_parallel"
	@echo "         ex: export run_parallel='mpirun -np 2'  "

oclrun: head_message $(mpirun_message)
	@echo "     Use the environment variable run_ocl"
	@echo "         ex: export run_ocl='on' to use OpenCL acceleration"
	@echo "     Use run_ocl='CPU' or 'ACC' to force use of hardware"
	@echo "         different than GPU"
	@echo "     and the environment variable ocl_platform"
	@echo "         ex: export ocl_platform='NVIDIA'"
	@echo "     and the environment variable ocl_devices"
	@echo "         ex: export ocl_devices='K20'"

CLEANFILES = report.py fldiff_yaml.py

all: head_message foot_message

check: report pytest

report: report.py
	@if test $(MAKELEVEL) = 0 ; then python ./report.py ; fi

pytest: $(pytest_DATA)
	py.test --nbval --disable-warnings $(srcdir)/InputVariablesTest.ipynb

