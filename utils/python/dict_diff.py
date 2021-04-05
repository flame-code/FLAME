def dict_diff(first, second):
#    """ Return a dict of keys that differ with another config object.  If a value is
#        not found in one fo the configs, it will be represented by KEYNOTFOUND.
#        @param first:   Fist dictionary to diff.
#        @param second:  Second dicationary to diff.
#        @return diff:   Dict of Key => (first.val, second.val)
#    """
    KEYNOTFOUND = '<KEYNOTFOUND>'       # KeyNotFound for dictDiff
    diff = {}
    # Check all keys in first dict
    for key in list(first.keys()):
        if (key not in second):
            diff[key] = (first[key], KEYNOTFOUND)
        elif (first[key] != second[key]):
            diff[key] = (first[key], second[key])
    # Check all keys in second dict to find missing
    for key in list(second.keys()):
        if (key not in first):
            diff[key] = (KEYNOTFOUND, second[key])
    return diff
