"""Handle and automatize the parsing of input arguments

This module uses the same convention of the :f:mod:`yaml_parse` fortran module to define the command line arguments.
Such module can be used to define command line arguments of python scripts that follows the same conventions
or to generate python functions that have the same signature than the provided command arguments.

"""

def get_python_function(target_kwargs_function,func_name,func_spec):
    """Convert a argparse spec into a python function

    This function provides a python function with a signature indicated by the ``fun_spec`` dictionary
    With the conventions of the :f:mod:`yaml_argparse` modules.
    The :py:func:`futile.Utils.function_signature_regenerator` function is used for the conversion

    Args:
       target_kwargs_function (func): the keyword arguments function we want to give the signature to.
       func_name (str): Name of the function, usually the key of the dictionary whose ``func_spec`` is the value
       func_spec (dict) : dictionary of the function specifications to be provided to the
          :py:func:`futile.Utils.function_signature_regenerator` function.

    Returns:
       func: the genreated function with signature given by the arguments of ``func_spec``
         defaulting to their default value.

    Todo:
       Create the docstring of the generated function by also including the docstring of the arguments
    """
    from copy import deepcopy
    from futile.Utils import function_signature_regenerator as fsr
    fspec=deepcopy(func_spec)
    docstring=fspec.pop("help")
    if "shorthelp" in fspec: fspec.pop("shorthelp")

    key_args = {key: val["default"] for (key, val) in fspec["args"].items()}

    return fsr(target_kwargs_function, fun_name=func_name,
               fun_docstring=docstring,**key_args)
