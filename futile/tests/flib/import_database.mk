#makefile portion to transform a set of yaml files into a 
#set of c source files including such files as a string buffer
#together with corresponding routine to be called from fortran
#to copy such buffer into a character array
# the list of the files have to be defined in the including makefile
# as a value of the YAML_DATABASE variable
# and the variable DATABASE_SOURCES have to be used in the corresponding
# places in the including makefile.
# Also the EXTRA_DIST and the CLEANFILES variable have to be defined

EXTRA_DIST += $(YAML_DATABASE)
DATABASE_SOURCES=$(YAML_DATABASE:.yaml=.c)
CLEANFILES += $(DATABASE_SOURCES)
$(DATABASE_SOURCES): $(YAML_DATABASE)
.yaml.c:
	@base=$(shell basename $< .yaml);\
	routinelc="get_"$$base;\
	routineuc=`echo $$routinelc | tr a-z A-Z`;\
	array=$$base"_arr";\
	file=$@;\
	echo "#include <config.h>" > $$file &&\
	echo "#define _GNU_SOURCE" >> $$file &&\
	echo "#include <fcntl.h>" >> $$file &&\
	echo "#include <stdio.h>" >> $$file &&\
	echo "#include <stdlib.h>" >> $$file &&\
	echo "#include <string.h>" >> $$file &&\
	echo "static const char $$array[] =" >> $$file &&\
	$(SED) -e "s/^/\"/;s/$$/\\\n\"/" $< >> $$file &&\
	echo "  ;" >> $$file &&\
	echo "void FC_FUNC_($$routinelc, $$routineuc)(char* db_ptr,int* db_len)" >> $$file &&\
	echo "{" >> $$file &&\
	echo "  if (*db_len==0)" >> $$file &&\
	echo "    {" >> $$file &&\
	echo "    *db_len=strlen($$array);" >> $$file &&\
	echo "    return;" >> $$file &&\
	echo "    }" >> $$file &&\
	echo "  else" >> $$file &&\
	echo "    {" >> $$file &&\
	echo "      memcpy(db_ptr,$$array, sizeof(char) * (*db_len));" >> $$file &&\
	echo "    }" >> $$file &&\
	echo "}" >> $$file
