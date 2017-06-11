import yaml
    
def kw_pop(*args,**kwargs):
    """ Treatment of kwargs. Eliminate from kwargs the tuple in args"""
    arg=kwargs.copy()
    key,default=args
    if key in arg:
        return arg,arg.pop(key)
    else:
        return arg,default

#function which removes from a set of lines the yaml_fields contained in the to_remove list
def clean_logfile(logfile_lines,to_remove):
  line_rev=logfile_lines #list of the lines of the logfile
  #loop in the reversed from (such as to parse by blocks)
  extra_lines=20 #internal variable to be customized
  line_rev.reverse()
  #clean the log
  cleaned_logfile=[]
  removed=[]
  #for line in line_rev: #line_iter:
  while len(line_rev) >0:
    line=line_rev.pop()
    to_print=line
    #check if the line contains interesting information
    for remove_it in to_remove :
      stream_list=[]
      #line without comments
      valid_line=line.split('#')[0]
      spaces='nospace'
      #control that the string between the key and the semicolon is only spaces
      if remove_it in valid_line and ":" in valid_line:
        #print "here",remove_it,remove_it in valid_line and ":" in valid_line,valid_line
        starting_point=valid_line.find(remove_it)
        tmp_buf=valid_line[:starting_point]
        #find the closest comma to the staring point, if exists
        tmp_buf=tmp_buf[::-1]
        starting_comma=tmp_buf.find(',')
        if starting_comma <0: st=0
        tmp_buf=tmp_buf[st:]
        tmp_buf=tmp_buf[::-1]
        tmp_buf=tmp_buf.strip(' ')
        #print "there",tmp_buf,'starting',starting_point,len(tmp_buf)
        valid_line= valid_line[starting_point+len(remove_it):]
        spaces= valid_line[1:valid_line.find(':')]
        #if remove_it+':' in line.split('#')[0]:
      if len(spaces.strip(' ')) == 0 and len(tmp_buf)==0: #this means that the key has been found
         #creates a new Yaml document starting from the line
         #treat the rest of the line following the key to be removed
         header=''.join(line.split(':')[1:])
         header=header.rstrip()+'\n'
         #eliminate the anchor
         header=header.lstrip(' ')
         header=header.lstrip('*') 
         if len(header) > 0 :
            stream_list.append(header)
         #part to be printed, updated
         to_print = line.split(':')[0] + ": <folded> \n"

         #then check when the mapping will end:
         while True:
            #create a stream with extra_lines block           
            for i in range(0,min(extra_lines,len(line_rev))):
               stream_list.append(line_rev.pop())
            #create a stream to be parsed
            stream=''.join(stream_list)
            #then parse the stream until the last valid position has been found
            try:
              for i in yaml.parse(stream,Loader=yaml.CLoader):
                endpos=i.end_mark.index
            except Exception, e:
              #  print 'error',str(e),stream
              #convert back the valid stream into a list
              #if needed the stream can be loaded into a document
              item_list=stream[:endpos].split('\n')
              #if lengths are different there is no need to add lines
              if len(item_list) != len(stream_list):
                #last line might be shorter, therefore treat it separately
                last_line=item_list.pop()
                #purge the stream
                for item in item_list:
                  stream_list.remove(item+'\n')
                #extract the remaining line which should be compared with the last one
                strip_size=len(last_line.rstrip())
                if strip_size > 0:
                  first_line=stream_list.pop(0)[strip_size:]
                  if '*' in first_line or '&' in first_line:
                    first_line='' #eliminate anchors
                else:
                  first_line=''
                #then put the rest in the line to be treated
                to_print.rstrip('\n')
                to_print += first_line+'\n'
                # the item has been found
                break
         stream_list.reverse()
         #put back the unused part in the document
         line_rev.extend(stream_list)
         # mark that the key has been removed
         if (remove_it not in removed):
           removed.append(remove_it)
           print 'removed: ',remove_it
  # then print out the line 
    cleaned_logfile.append(to_print)

  # check that everything has been removed, at least once
  if (set(removed) != set(to_remove)):
    print 'WARNING, not all the requested items have been removed!'
    print 'To_remove : ',to_remove
    print 'removed   : ',removed
    print 'Difference: ',list(set(to_remove) - set(removed) )
  return cleaned_logfile

def load(file=None,stream=None):
    strm=stream if stream else open(file,'r')
    try:
        ld=yaml.load(strm,Loader=yaml.CLoader)
    except:
        ld=yaml.load_all(strm,Loader=yaml.CLoader)
    return ld

def dump(data,filename=None,raw=False,tar=None):
    todump=str(data) if raw else yaml.dump(data)
    if filename:
        if tar:
            import tarfile
            from cStringIO import StringIO
            import time
            f=tarfile.TarInfo(filename)
            f.size=len(str(todump))
            f.mtime=time.time()
            tar.addfile(f,StringIO(str(todump)))
        else:
            f=open(filename,'w')
            f.write(todump)
    else:
        import sys
        sys.stdout.write(todump)

class YamlDB(dict):
    """Yaml Database, read from a file or a stream
    Attributes:

    file: take the database from a file
    stream: associate the stream as the value of the dictionary
    singledoc: guarantees that the provided stream can only contain one document
    ignore: A list of keys that will not be considered in the loading.
            Useful to parse log logfiles in less time
    YamlDB(file) -> Database from a file
    YamlDB(stream) -> Database from a stream
    YamlDB() -> new empty dictionary
    YamlDB(mapping) -> new dictionary initialized from a mapping object's
              (key, value) pairs
    YamlDB(key, stream=stream) -> new YamlDB with unparsed stream written as the value of the key.
    YamlDB(iterable) -> new dictionary initialized as if via:
           d = {}
           for k, v in iterable:
               d[k] = v
    YamlDB(**kwargs) -> new dictionary initialized with the name=value pairs
      in the keyword argument list.  For example:  YamlDB(one=1, two=2)
    """
    def __init__(self,*args,**kwargs): #stream=None,singledoc=False):
        """Extract the database information"""
        newkw,self.singledoc=kw_pop('singledoc',False,**kwargs)
        newkw,self.file=kw_pop('file',None,**newkw)
        newkw,self.stream=kw_pop('stream',None,**newkw)
        newkw,self.ignore=kw_pop('ignore',None,**newkw)
        listargs=list(args)
        if len(listargs)>0:
            if type(listargs[0])==type('string'):
                import os
                tmp=listargs.pop(0)
                if self.file is None and os.path.isfile(tmp):
                    self.file=tmp
                elif self.stream is None:
                    self.stream=tmp
        #suppose that the stream is a file
        if self.file is not None:
            sl=open(self.file,'r')
            prestream=''.join(sl.readlines())
            sl.close()
        elif self.stream is not None:
            prestream=self.stream
        else:
            super(YamlDB,self).__init__(*listargs,**newkw)
            return
        if self.singledoc:
            self.stream=prestream
        else:
            #process stream, only one document
            startpos=0
            i=0
            for endpos in self._doc_finder(prestream):
                if i==1: self.docs=[YamlDB(stream=doctmp,singledoc=True,ignore=self.ignore)]
                doctmp=prestream[startpos:endpos]
                if i>=1: self.docs.append(YamlDB(stream=doctmp,singledoc=True,ignore=self.ignore))
                startpos=endpos
                i+=1
            if i==1: self.stream=doctmp
        if self.stream is not None: 
            dd=self._load()
            if dd: super(YamlDB,self).update(dd)
    def __len__(self):
        if hasattr(self,'docs'):
            return len(self.docs)
        else:
            return 0
    def __getitem__(self,key):
        if isinstance(key,int):
            if hasattr(self,'docs'):
                return self.docs[key]
        else:
            return self.get(key)
    def documents(self):
        "Generator over the whole documents"
        if hasattr(self,'docs'):
            for d in self.docs:
                yield dict(d)
        else:
            yield dict(self)
    def _doc_finder(self,stream):
        #first find the stream of the document start
        startpos=0
        import sys
        print 'Loading...'
        while startpos < len(stream):
            try:
                startpos+=self._event_finder(stream[startpos:],yaml.events.DocumentEndEvent)
            except:
                startpos=len(stream)
            print int((100.0*startpos)/len(stream)),'%'
            yield startpos
    def _load(self):
        #clean the logfile if needed
        if self.ignore is not None:
            listlog=[ a+'\n' for a in self.stream.split('\n')]
            cleanedlog=clean_logfile(listlog,self.ignore)
            stream='\n'.join(cleanedlog)
        else:
            stream=self.stream
        try:
            return yaml.load(stream,Loader=yaml.CLoader)
        except Exception, e:
            #here we might put a bigger 
            return None
    def _event_finder(self,present,event,end=False):
        for i in yaml.parse(present,Loader=yaml.CLoader):
            if isinstance(i,event):
                if end:
                    key=i.value
                    startpos=i.start_mark.index
                    endpos=startpos+self._find_endblock(present[startpos:])
                    substream=present[i.end_mark.index:endpos]
                    return key,substream,endpos
                else:
                    startpos=i.end_mark.index
                    return startpos
        #in the case of absent event
        if end:
            print 'here',present,'end'
            return 'nothing','',len(present)
        else:
            return len(present)  
    def _find_endblock(self,stream):
        """Find the end of the block which is yaml-compliant"""
        endpos=0
        try:
            for i in yaml.parse(stream,Loader=yaml.CLoader):
                endpos=i.end_mark.index
        except yaml.YAMLError,e:
            #stop at the last carriage return
            endpos=e.problem_mark.index
            endpos=stream.rfind('\n',0,endpos)
        return endpos
