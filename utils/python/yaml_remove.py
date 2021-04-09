#!/usr/bin/env python
import yaml
import argparse

try :
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader

str1 = "This script reads input.yaml and deletes the specific structures from input file"
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_del', action='store' ,type=str, help="Name of the delete list in yaml format")
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input structures file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output structures file in yaml format")
args=parser.parse_args()

stream = open(args.fn_del,'r')
stream_2=open(args.fn_inp)
output_stream=open(args.fn_out,"w")

del_list = yaml.load(stream)
for i,d in enumerate(del_list):
    del_list[i]=int(d)-1
    print((del_list[i]), end=' ')
    print(('\n'), end=' ')

dict_list = list(yaml.load_all(stream_2, Loader=Loader))
conf_list=list(range(len(dict_list)))

print(('\n'), end=' ')
print(('\n'), end=' ')
print(('\n'), end=' ')

conf_keep_list=list(set(conf_list) - set(del_list))
#print(conf_keep_list)

for i in conf_keep_list:
    print((i), end=' ')
    print(('\n'), end=' ')

for conf in conf_keep_list:
    output_stream.write('---\n')
    yaml.dump(dict_list[conf],output_stream)
