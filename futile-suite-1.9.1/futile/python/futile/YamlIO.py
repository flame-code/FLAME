import yaml
from futile.Utils import write

#function which removes from a set of lines the yaml_fields contained in the to_remove list
def clean_logfile(logfile_lines,to_remove):
  """Remove yaml fields from a list of lines.

  Removes from a set of lines the yaml_fields contained in the to_remove list.

  Arguments:
      logfile_lines (list): list of the lines of the logfile. Generated from a file by e.g. :py:meth:`~io.IOBase.readlines`.
      to_remove (list): list of keys to remove from logfile_lines

  Returns:
      list of lines where the removed keys have as values the `"<folded>"` string
  """
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
            except Exception(e):
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
           write('removed: ',remove_it)
  # then print out the line
    cleaned_logfile.append(to_print)

  # check that everything has been removed, at least once
  if (set(removed) != set(to_remove)):
    write('WARNING, not all the requested items have been removed!')
    write('To_remove : ',to_remove)
    write('removed   : ',removed)
    write('Difference: ',list(set(to_remove) - set(removed) ))
  return cleaned_logfile


def load_from_archive(arch, member=None):
    import tarfile
    tar = tarfile.open(arch)
    members = [tar.getmember(member)] if member is not None else tar.getmembers()
    # print members
    for memb in members:
        f = tar.extractfile(memb)
        dicts += load(stream=f.read())
        # Add the label (name of the file)
        # dicts[-1]['label'] = memb.name


def load(file=None, stream=None, doc_lists=True, safe_mode=False,
         archive=None, member=None):
    """Encapsulate the loading of yaml documents.

    Provides a dictionary, or a list of dictionaries, which
    represents the structure of the stream to be loaded.
    It also wraps the yaml loader to perform a optimized parsing when the
    `minloader` of PyYaml 3.13 is available.
    This wrapper ensures to extract from the stream the maximum possible
    information by choosing the best loader available.

    Arguments:
        file (str): path of the yaml-compliant file containing the stream
             to be loaded

        stream (str): the stream to load, overrides the ``file`` argument
             if present

        archive (str): path of the archive to be used for the retrieval
             of the stream

        member (str): name of the file to be extracted from the archive.
            the entire archive is parsed if absent

        doc_lists (bool): if True, ensures that the results is always in a form
           of lists of documents, even in the case of a single doc
           When False, the return type is either a dictionary or a generator
           according to the specifications of yaml.load and
           yaml.load_all respectively.

        safe_mode (bool): When true, in the case of multiple documents
           in the stream, it loads the document one after another.
           This is useful to avoid losing of all the document list
           in the case when one of the document is
           not yaml compliant, like in the case of a broken logfile.
           It may works only when the separation of the
           documents is indicated by the usual syntax ``"---\\n"``
           (i.e. no yaml tags between documents)

    Returns:
        * a list of dictionaries, if ``doc_lists`` is set to ``True``;
        * a dictionary, if the stream or the file, or the archive contains a
             single yaml document;
        * a generator if the parsed stream is made of multiple
             documents *and* ``safe_mode`` = ``False``;
        * a list of dictionaries if the stream is made of multiple documents
             and ``safe_mode`` is ``True``.
    """
    # choose the loader
    try:
        ldr = yaml.MinLoader
    except Exception as e:
        try:
            ldr = yaml.CLoader
        except Exception as f:
            ldr = yaml.Loader

    # load the documents
    ld = []

    # verify if we are in the archive case
    if archive is not None:
        import tarfile
        tar = tarfile.open(archive)
        if member is not None:
            members = [tar.getmember(member)]
        else:
            members = tar.getmembers()
        ld = []
        for memb in members:
            f = tar.extractfile(memb)
            if safe_mode:
                try:
                    ll = yaml.load(f.read(), Loader=ldr)
                except Exception as f:
                    write('Document', member, 'of archive NOT loaded, error:',
                          f)
            else:
                ll = yaml.load(f.read(), Loader=ldr)
            ld += ll
            if len(members) == 1 and not doc_lists:
                return ll
        return ld

    # Detect None otherwise a doc == '' gives an error
    strm = stream if stream is not None else open(file, 'r').read()
    try:
        ld = yaml.load(strm, Loader=ldr)
        if doc_lists:
            ld = [ld]
    except Exception as e:
        if safe_mode:
            ld = []
            documents = [v for v in strm.split('---\n') if len(v) > 0]
            for i, raw_doc in enumerate(documents):
                try:
                    ld.append(yaml.load(raw_doc, Loader=ldr))
                except Exception as f:
                    write('Document', i, 'of stream NOT loaded, error:', f)
        else:
            ld = yaml.load_all(strm, Loader=ldr)
            if doc_lists:
                ld = [l for l in ld]
    return ld


def dump(data, filename=None, raw=False, tar=False):
    """Encapsulate the dumping of dictionaries.

    This function is useful to dump a dictionary in yaml or json-compliant form.
    This may be used as an alternative to the usual
    :py:meth:`yaml.dump <https://pyyaml.org/wiki/PyYAMLDocumentation>` method,
    especially when the dictionary to be dump'ed is heavy.
    No particular attention is paid in human readability of the output.
    The dumped information can then be parsed either from json or yaml
    interpreter.

    Arguments:
       data (dict,list): the information to be dumped
       filename (str): path of the file in which the information will
          be stored.
          If absent, the information is written on :py:func:`sys.stdout`.
       raw (bool): if ``True`` the output is in json-style, otherwise it is
          pre-processed by :py:meth:yaml.dump,
          but ``None`` is passed to ``default_flow_style``.
       tar (bool): if ``True`` the filename is assumed to be a compressed
          tarfile. The :py:mod:`tarfile` module is used to create and
          append information.
    """
    todump = str(data) if raw else yaml.safe_dump(data,
                                                  default_flow_style=None)
    if filename:
        if tar:
            import tarfile
            from cStringIO import StringIO
            import time
            f = tarfile.TarInfo(filename)
            f.size = len(str(todump))
            f.mtime = time.time()
            tar.addfile(f, StringIO(str(todump)))
        else:
            f = open(filename, 'w')
            f.write(todump)
    else:
        import sys
        sys.stdout.write(todump)


class YamlDB(dict):
    """
    Yaml Database, read from a file or a stream

    :param str file: take the database from a file
    :param str stream: associate the stream as the value of the dictionary
    :param bool singledoc: guarantees that the provided stream can only contain one document
    :param list ignore: A list of keys that will not be considered in the loading.
            Useful to parse log logfiles in less time

    * YamlDB(file) -> Database from a file
    * YamlDB(stream) -> Database from a stream
    * YamlDB() -> new empty dictionary
    * YamlDB(mapping) -> new dictionary initialized from a mapping object's (key, value) pairs
    * YamlDB(key, stream=stream) -> new YamlDB with unparsed stream written as the value of the key.
    * YamlDB(iterable) -> new dictionary initialized as if via:
           d = {}
           for k, v in iterable:
               d[k] = v
    * YamlDB(**kwargs) -> new dictionary initialized with the name=value pairs
      in the keyword argument list.  For example:  YamlDB(one=1, two=2)
    """
    def __init__(self,*args,**kwargs): #stream=None,singledoc=False):
        """Extract the database information"""
        from futile.Utils import kw_pop
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
        write('Loading...')
        while startpos < len(stream):
            try:
                startpos+=self._event_finder(stream[startpos:],yaml.events.DocumentEndEvent)
            except:
                startpos=len(stream)
            write(int((100.0*startpos)/len(stream)),'%')
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
        except Exception(e):
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
            write( 'here',present,'end')
            return 'nothing','',len(present)
        else:
            return len(present)
    def _find_endblock(self,stream):
        """Find the end of the block which is yaml-compliant"""
        endpos=0
        try:
            for i in yaml.parse(stream,Loader=yaml.CLoader):
                endpos=i.end_mark.index
        except yaml.YAMLError(e):
            #stop at the last carriage return
            endpos=e.problem_mark.index
            endpos=stream.rfind('\n',0,endpos)
        return endpos
