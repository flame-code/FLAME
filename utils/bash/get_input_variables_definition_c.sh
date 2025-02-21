#!/bin/bash

FILE_YAML=src/input_variables_definition.yaml
FILE_C=src/input_variables_definition.c
SED=sed

top_srcdir=$PWD

base=$(basename $FILE_YAML .yaml);\
datafilesdir=`echo "$top_srcdir/datafiles" | $SED -e 's!\/!\\\/!g'`;\
routinelc=`echo get_$base | tr A-Z a-z`;\
routineuc=`echo $routinelc | tr a-z A-Z`;\
array=$base"_arr";\
file=$@;\
echo "#include <config.h>" > $FILE_C &&\
echo "#define _GNU_SOURCE" >> $FILE_C &&\
echo "#include <fcntl.h>" >> $FILE_C &&\
echo "#include <stdio.h>" >> $FILE_C &&\
echo "#include <stdlib.h>" >> $FILE_C &&\
echo "#include <string.h>" >> $FILE_C &&\
echo "static const char input_variables_definition_arr[] =" >> $FILE_C &&\
$SED -e "s/^/\"/;s/$/\\\n\"/" $FILE_YAML >> $FILE_C &&\
$SED -i "s/DATAFILESDIR/$datafilesdir/" $FILE_C &&\
echo "  ;" >> $FILE_C &&\
echo "void FC_FUNC_($routinelc, $routineuc)(char* db_ptr,int* db_len)" >> $FILE_C &&\
echo "{" >> $FILE_C &&\
echo "  if (*db_len==0)" >> $FILE_C &&\
echo "    {" >> $FILE_C &&\
echo "    *db_len=strlen($array);" >> $FILE_C &&\
echo "    return;" >> $FILE_C &&\
echo "    }" >> $FILE_C &&\
echo "  else" >> $FILE_C &&\
echo "    {" >> $FILE_C &&\
echo "      memcpy(db_ptr,$array, sizeof(char) * (*db_len));" >> $FILE_C &&\
echo "    }" >> $FILE_C &&\
echo "}" >> $FILE_C
