"""
This file contains some low-level useful functions
"""

from __future__ import print_function


def write(*args, **kwargs):
    """
    Wrapper for print function or print to ensure compatibility with python 2
    The arguments are used similarly as the print_function
    They can also be generalized to python 2 cases
    """
    return print(*args, **kwargs)


def push_path(inp, *keys):
    """
    Follow in the dictionary inp the path indicated by the keys.
    If this path does not exists creates it.

    Args:
       inp (dict): dictionary
       keys (str): keys of the path to follow

    Returns:
       (``branch``,``key``) tuple, where

       * ``branch`` (dict): the dictionary of the second-last item of the path
       * ``key`` (str): the last item of the path

    Example:

       >>> inp={}
       >>> d,key=push_path(inp,'dft','nspin','mpol')
       >>> print (d,key)
       >>> print (inp)
       {},'mpol'
       {'dft': {'nspin': {}}}

       >>> inp={'dft': {'nspin': {'mpol': 2}}}
       >>> d,key=push_path(inp,'dft','nspin','mpol')
       >>> print (d,key)
       >>> print (inp)
       {'mpol': 2},'mpol'
       {'dft': {'nspin': {'mpol': 2}}}

    """
    tmp = inp
    for i, key in enumerate(keys):
        k = key
        if i == len(keys)-1:
            break
        tmp.setdefault(key, {})
        tmp = tmp[key]
    return tmp, k


def dict_set(inp, *subfields):
    """Ensure the provided fields and set the value

    Provide a entry point to the dictionary.
    Useful to define a key in a dictionary that may not have the
    previous keys already defined.

    Arguments:
       inp (dict): the top-level dictionary
       subfields (str,object): keys, ordered by level, that have to be
          retrieved from topmost level of ``inp``.
          The last item correspond to the value to be set .

    Example:

       >>> inp={}
       >>> dict_set(inp,'dft','nspin','mpol',2)
       >>> print (inp)
       {'dft': {'nspin': {'mpol': 2}}}

    """
    if len(subfields) <= 1:
        raise ValueError('invalid subfields, the sequence should be longer than one item as the last one is the value to be given')
    keys = subfields[:-1]
    tmp, key = push_path(inp, *keys)
    tmp[key] = subfields[-1]


def dict_get(inp, *subfields):
    """Find the value of the provided sequence of keys in the dictionary,
    if available.

    Retrieve the value of the dictionary in a sequence of keys if it is
    available. Otherwise it provides as default value the last item of the
    sequence ``subfields``.

    Args:
       inp (dict): the top-level dictionary. Unchanged on exit.
       subfields (str,object): keys, ordered by level, that have to be
           retrieved from topmost level of ``inp``. The last item correspond
           to the value to be set.

    Returns:
       The value provided by the sequence of subfields if available,
       otherwise the default value given as the last item of the ``subfields``
       sequence.

    """
    if len(subfields) <= 1:
        raise ValueError('invalid subfields, the sequence should be longer than one item as the last one is the value to be given')
    tmp = inp
    keys = subfields[:-1]
    val = subfields[-1]
    for key in keys:
        tmp = tmp.get(key)
        if tmp is None:
            return val
    return tmp


def sort_lists(sort_by, ascending, *lists):
    """
    Sort lists altogether following the lists indicated by the ``sort_by`` index.

    Args:

       sort_by (int):  the index of the list which has to be taken as reference for sorting
       ascending (bool): Sort is performed in ascending order if True

       *lists: sequence of lists to be mutually sorted. They have to be of the same length.

    Returns:
       tuple of sorted lists

    Example:
    >>> l1=[5,3,4]
    >>> l2=['c','t','q']
    >>> l3=[6,3,7]
    >>> print (sort_lists(0,True,l1,l2,l3))
    >>> print (sort_lists(2,True,l1,l2,l3))
    [(3, 4, 5), ('t', 'q', 'c'), (3, 7, 6)]
    [(3, 5, 4), ('t', 'c', 'q'), (3, 6, 7)]
    """
    import operator
    return zip(*sorted(zip(*lists), reverse=not ascending,
               key=operator.itemgetter(sort_by)))


def file_time(filename):
    """
    Determine the time of the last modification of a file.

    Args:
       filename (str): path of the file to inspect.

    Returns:
       float: time of the modified file. Returns 0 if the file does not exist.
    """
    import os
    if os.path.isfile(filename):
        return os.path.getmtime(filename)
    else:
        return 0


def non_null_size(filename):
    """
    Control if the file has nonzero size and exists
    """
    from os.path import getsize, isfile
    return isfile(filename) and (getsize(filename) > 0)


def more_recent_than_parent(filename, parent):
    """
    Filename should be more recent than parent
    """
    t = non_null_size(filename) and (file_time(parent) <= file_time(filename))
    return t


def fill_dictionary_in_parallel(nthreads, keys, func, **kwargs):
    """
    Fill a dictionary of a given set of keys with the return value
    of a function which accept this key as a first argument

    Args:
        nthreads(int):  the number of threads of the pool
        keys(list): the arguments of the function. Will be the list of the
            dictionary
        func(func): the python function that has the key as a last argument
        **kwargs: further arguments of the function, if needed

    Returns:
        dict: the key-> obj dictionary with the obj the return value of func
    """
    from multiprocessing import Pool
    import time
    from functools import partial
    refunc = partial(func, **kwargs)
    p = Pool(nthreads)
    start = time.time()
    res = p.map(refunc, keys)
    p.close()
    p.join()
    end = time.time()
    write(end - start)
    return {a: b for a, b in zip(keys, res)}


class Node():
    """
    An object that is associated to a queue.
    Has a requirement and a validity function as well as a generation function
    that is triggered in case the node is not valid.

    Args:
       obj (object): a generic object
    """
    def __init__(self, obj, valid_if=None):
        self.obj = obj
        if valid_if is not None:
            self.valid_if(func=valid_if)

    def valid_if(self, func):
        """
        Set the function as callback to check the validity of the node.

        Args:
           func(func): function that has the node object as a the first
              argument and returns a boolean which assess the validity
              of the function.
        """
        self.validity_func = func

    @property
    def valid(self):
        if hasattr(self, 'parent'):
            return self.validity_func(self.obj, self.parent.obj)
        else:
            return self.validity_func(self.obj)

    def requires(self, node, generator):
        """
        Set the dependency between two nodes.

        Args:
           node (Node): the node from which this one depends.
              If this is a valid node, and the present is not, it should employ
              the generator function for creating the node.
           generator (func): function which should be called to make the
              node valid. Should have as arguments the two objects associated
              the node (the current one is the first argument)
        """
        self.parent = node
        self.dependency_func = generator

    def validate(self):
        """
        This function makes the node valid by validating its dependency
        and by calling the generator function if this is not the case.

        Returns:
           bool: the value of self.valid. Should be true.
        """
        if hasattr(self, 'parent'):
            assert self.parent.validate()
        if not self.valid:
            self.dependency_func(self.obj, self.parent.obj)
        assert self.valid
        return self.valid


def dict_merge(dest, src):
    """ Recursive dict merge. Inspired by :meth:`dict.update`, instead of
    updating only top-level keys, dict_merge recurses down into dicts nested
    to an arbitrary depth, updating keys. The ``src`` is merged into
    ``dest``.  From :ref:`angstwad/dict-merge.py
    <https://gist.github.com/angstwad/bf22d1822c38a92ec0a9>`

    Arguments:
       dest (dict): dict onto which the merge is executed
       src (dict): dict merged into dest

    """
    import collections
    for k, v in src.items():
        if (k in dest and isinstance(dest[k], dict)
                and isinstance(src[k], collections.Mapping)):
            dict_merge(dest[k], src[k])
        else:
            dest[k] = src[k]


def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


def unique_filename(prefix):
    """
    Provides a filename with a unique id appended

    Args:
        prefix (str): the prefix of the file

    Returns:
        str: filename
    """
    from uuid import uuid1
    unique_id = str(uuid1().hex)
    return prefix + unique_id


def file_list(directory, suffix=None, prefix=None, exclude=None,
              include_directory_path=False):
    """
    Return the list of the files inside a given directory

    Args:
       directory (str): path of the directory to search into
       suffix (str): the suffix that the files should have
       prefix (str): the prefix that the files should have
       exclude (str): exclude the files which matches this string from the list
       include_directory_path (bool): if True includes the path of
           the directory in the list.

    Returns:
       list: list of the files that matches the requirements.
    """
    files = []
    from os import listdir
    from os.path import join
    for filename in listdir(directory):
        ok = exclude not in filename if exclude is not None else True
        if ok and suffix is not None:
            ok = filename.endswith(suffix)
        if ok and prefix is not None:
            ok = filename.startswith(prefix)
        if ok:
            if include_directory_path:
                to_append = join(directory, filename)
            else:
                to_append = filename
            files.append(to_append)
    return files


def make_dict(inp):
    """
    Transform the instance ``inp`` into a python dictionary.
    If inp is already a dictionary, it perfroms a copy.

    Args:
       inp (dict): a instance of a Class which inherits from dict

    Returns:
       dict: the copy of the class, converted as a dictionary
    """
    import copy
    local_tmp = copy.deepcopy(inp)
    local_input = {}
    local_input.update(local_tmp)
    return local_input


def execute_code_if(condition, code, glob=None, loc=None):
    """
    Execute code if condition is true

    Args:
        condition (bool): if true the code is executed
        code_if_obj_not_found (str): the code to be evaluated if the object has
             not been found.
        globals (dict): the global variables dictionary
        locals (dict): the local variables dictionary

    Returns:
         the object returned from the code executed, None otherwise
    """
    if not condition:
        return None
    if glob is None:
        glob = globals()
    if loc is None:
        loc = locals()
    return eval(code, glob, loc)


def ensure_object(filename, code_if_obj_not_found=None, glob=None, loc=None):
    """
    Identify a pickle file to save a given object on it.
    In case this file is present, read the object from it.
    Otherwise, assume that the object is ready to be dumped and write
    it in the file.

    Args:
       filename (str): the path of the file in which the object is saved/loaded
       code_if_obj_not_found (str): the code to be evaluated if the object has
       not been found.
       globals (dict): the global variables dictionary
       locals (dict): the local variables dictionary
    Returns:
       object: The result from the pickle file or from the code to be executed
    """
    import pickle
    try:
        with open(filename, "rb") as ifile:
            obj = pickle.load(ifile)
    except Exception as e:
        obj = execute_code_if(True, code_if_obj_not_found, glob=glob, loc=loc)
        with open(filename, "wb") as ofile:
                pickle.dump(obj, ofile)
    return obj


def property_attribute(self, attribute, code_if_not_found):
    obj = execute_code_if(not hasattr(self, attribute), code_if_not_found)
    if obj is None:
        obj = getattr(self, attribute)
    return obj


def function_signature_regenerator(target_kwargs_function, fun_name='',
                                   fun_docstring='', **kwargs):
    '''
    Generate the function of the name provided by `fun_name`,
    with signature provided by the kwargs dictionary.

    Args:
       target_kwargs_function (func): keyword arguments function that will be
           used for the generated function.
       fun_name (str): name of the regenerated function. If empty it will be
           the ``target_kwargs_functon.__name__`` prefixed by ``regenerated``,
           which will be copied in the docstring of the regenerated function.
       fun_docstring (str): docstring of the generated function, if empty it
           will take the docstring from ``target_kwargs_function``.
       **kwargs: keyword arguments which will represent the signature of the
           generated function.

    Example:
        >>> def write_kwargs(**kwargs):
        >>>     """
        >>>     Convert keyword arguments into a string
        >>>     """
        >>>     return str(kwargs)
        >>> write_opts=function_signature_regenerator(write_kwargs,
        >>>                                           fun_name='write_opts',
        >>>                                           opt1='default1',
        >>>                                           opt2='default2')
        >>> help(write_opts)
        >>> print (write_opts())
        Help on function write_opts:

        write_opts(opt1='default1', opt2='default2')
              Convert keyword arguments into a string

        {'opt1': 'default1', 'opt2': 'default2'}

    '''
    signature = option_line_generator(',', **kwargs).lstrip(',')
    docstring = target_kwargs_function.__doc__ if not fun_docstring else fun_docstring
    if docstring is None:
        docstring = "Automatically generated function from the target function '" + target_kwargs_function.__name__ + "'"
    docstring = '   """\n'+docstring+'\n   """'
    fname = "regenerated_" + target_kwargs_function.__name__ if not fun_name else fun_name
    function = "def %s(%s):\n%s\n   return target_function(**locals())" % (fname, signature, docstring)
    gen_locals = {}
    gen_object = compile(function, 'generated_fun', 'exec')
    eval(gen_object, {'target_function': target_kwargs_function}, gen_locals)
    return gen_locals[fname]


def option_line_generator(separator='--', **kwargs):
    """
    Associate to each of the keyword arguments a command line argument.

    Args:
       separator (str): The string needed to separate the options.
       Might be '--' for command-line arguments, but also ','
       for function signatures.

    Warning:
        The separator comes **before** the first argument therefore pay
        attention to lstrip it in case you want to use it as a function
        signature string.

    Example:
        >>> option_line_generator(arg1='val1',arg2='val2')
        '--arg1=val1 --arg2=val2'
    """
    command = ''
    for option, value in kwargs.items():
        command += separator + option + '="' + str(value) + '" '
    return command


def kw_pop(*args, **kwargs):
    """
    Treatment of kwargs. Eliminate from kwargs the tuple in args.
    """
    arg = kwargs.copy()
    key, default = args
    if key in arg:
        return arg, arg.pop(key)
    else:
        return arg, default


def find_files(regexp, archive=None):
    """
    Returns a list of the paths to the files that follow the regular expression
    regexp. They are searched from the current working directory or from an
    archive given as optional argument.


    :param regexp: A regular expression
    :type regexp: string
    :param archive: an opened tarfile archive (optional)
    :type archive:
    :returns: a list of all the paths that agree with the regexp
    :rtype: list of strings
    :raises: ValueError if the regexp does not find a single path.


    Example::

        #Find all python files in the current working directory
        find_files('*py')

        #An exmple outside of the current working directory
        find_files('*/log-*.yaml')

        #Example using a tarfile
        import tarfile
        my_archive = tarfile.open('archive.tar.gz')
        find_files('*/*/log-*.yaml', archive=my_archive)
    """
    import os

    # Get a list of all paths to files satisfying the regexp
    if archive is not None:
        paths = _find_files_from_archive(regexp, archive)
    else:
        paths = os.popen('ls '+regexp).read().splitlines()

    # Test that the regexp found files
    if paths == []:
        raise ValueError('The regexp "{}" leads to no file. '\
                         'Consider using another one.'.format(regexp))
    else:
        return paths


def _find_files_from_archive(re, archive):
    """
    This function retrieves the list of Logfiles instances
    from the file archived satisfying a regular expression.
    #function to identify an archive out of its regexp,
    #solves the bug in re for '*' (solved in Python 2.7.6)
    """
    import tarfile

    # Open the archive
    with tarfile.open(archive, 'r') as arch:
    # Return paths to logfiles satisfying the regexp
        return [f for f in arch.getnames()
                if all(pattern in f for pattern in re.split('*'))]


def ensure_copy(src, dest):
    """Copy src into dest.

    Guarantees that the file indicated by ``dest`` is a copy of the file ``src``

    Args:
      src (str): path of the source file. Should be valid.
      dest (src): path of the destination file

    Returns:
      bool: ``True`` if the file needed to be copied, ``False`` if ``src``
         and ``dest`` are identical
    """
    import shutil, os
    copied = False
    if (os.path.isfile(dest) and os.stat(dest) != os.stat(src)) or not os.path.isfile(dest):
        shutil.copy2(src, os.path.dirname(dest))
        copied = True
    return copied


def ensure_dir(file_path):
    """
    Guarantees the existance on the directory given by the (relative) file_path

    Args:
       file_path (str): path of the directory to be created

    Returns:
       bool: True if the directory needed to be created, False if it existed already
    """
    import os
    directory = file_path
    created=False
    if not os.path.exists(directory):
        os.makedirs(directory)
        created=True
    return created


if __name__ == '__main__':
    import os

    #Tests of the find_files function
    #
    print("Test finding all python files in this directory")
    print(find_files("*py"))
    print()

    #
    print("Test finding the Utils.py file in this directory")
    print(find_files("Utils.py"))
    print()

    #
    print("Test raising a ValueError because the regexp leads to no files")
    try:
        find_files('*html')
    except ValueError as e:
        print('This raised the following ValueError:')
        print(e)
    print()

    #
    print("Test raising an exception because there is no such archive.")
    fname = 'file.tar.gz'
    if fname in os.popen('ls'): os.system('rm '+fname)
    #os.system('rm '+fname)
    try:
        find_files('*py', archive=fname)
    except Exception as e:
        #print(dir(e))
        print('This raised the following Exception:')
        print(e)
    print()

    #
    print("Test without error using an archive")
    os.system('find * -name "*py" | tar -zcvf '+fname+' -T -')
    os.system('ls '+fname)
    find_files('*py', archive=fname)
    os.system('rm '+fname)


def create_tarball(filename, files, objects={}):
    """
    Assemble files and objects in a tarball

    Args:
        filename (str): the name of the archive. Determine the tarball
           compression method from its extension.
        files(dict,set): file paths that have to be included in the tarball.
           If it is a dict, it should be in the form "{arcname : file}",
           where `file` is the path of the file to be put, and `arcname`
           is the name of thefile that would be used in the archive.
           If it is a set, the file will preserve its name in the archive
        objects(dict): dictionary '{arcname: buffer}' of the buffers that
           will have to be serialized in the `arcname` position
           the buffers are given as `class::StringIO.StringIO` instance,
           following specification of the `func:serialize_objects` function.
    """
    import tarfile
    from os.path import splitext
    extension = splitext(filename)[-1].lstrip('.')
    arch = tarfile.open(filename, mode='w:'+extension)
    isdict = isinstance(files, dict)
    for arcname in files:
        name = files[arcname] if isdict else arcname
        arch.add(name=name, arcname=arcname)
    for arcname, string in objects.items():
        tarinfo = tarfile.TarInfo(arcname)
        tarinfo.size = len(string.getvalue())
        arch.addfile(tarinfo=tarinfo, fileobj=string)
    arch.close()


def unpack_tarball(archive, tmpdir_prefix='tmp_'):
    """
    Open an archive in a temporary directory

    Args:
         archive (str): the path of th archive to open
         tmpdir_prefix (str): prefix of the temporary directory to untar the
             archive to.
    Returns:
         tuple: tmpdir, files path of the temporary directory and names of
              the files extracted by the tarfile
    """
    from os.path import basename
    import tempfile
    import tarfile
    # creates the directory
    tmpdir = tempfile.mkdtemp(prefix=tmpdir_prefix + basename(archive) + '_')
    # extract the archive
    arch = tarfile.open(archive)
    arch.extractall(path=tmpdir)
    files = arch.getnames()
    arch.close()
    return tmpdir, files


def serialize_objects(objects, extra_encoder_functions=[]):
    """
    Convert a dictionary of objects into buffers.
    Employs json serialization into `StringIO` instances

    Args:
       objects(dict): dictionary of key/value pair of objects to be serialized
       extra_encoder_functions (dict)): dictionary of the format
           {'cls': Class, 'func': function} which is employed in the
           serialization

    Returns:
       dict: dictionary of key/buffer pairs
    """
    from io import BytesIO
    import json

    class CustomEncoder(json.JSONEncoder):
        """ Special json encoder for numpy types and System"""
        def default(self, obj):
            import numpy as np
            if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                                np.int16, np.int32, np.int64, np.uint8,
                                np.uint16, np.uint32, np.uint64)):
                return int(obj)
            elif isinstance(obj, (np.float_, np.float16, np.float32,
                                  np.float64)):
                return float(obj)
            elif isinstance(obj, (np.ndarray,)):
                return obj.tolist()
            elif isinstance(obj, (set,)):
                return list(obj)
            else:
                for spec in extra_encoder_functions:
                    if isinstance(obj, (spec['cls'],)):
                        return spec['func'](obj)
            return json.JSONEncoder.default(self, obj)

    return {key: BytesIO((json.dumps(obj, cls=CustomEncoder)).encode("utf-8"))
            for key, obj in objects.items()}
