#!/usr/bin/env python

from f2py import crackfortran
import sys

def toCType(tt):
  if (tt['typespec'] == 'logical'):
    return 'bool'
  elif (tt['typespec'] == 'integer'):
    kind = tt.get('kindselector', {'kind': 4})['kind']
    if kind in ['f_long', '8']:
      return 'long'
    else:
      return 'int'
  elif (tt['typespec'] == 'real'):
    kind = tt.get('kindselector', {'kind': 4})['kind']
    if kind in ['dp', 'gp', 'wp', 'f_double', '8']:
      return 'double'
    else:
      return 'float'
  elif (tt['typespec'] == 'character'):
    return 'char'
  elif (tt['typespec'] == 'type'):
    if 'pointer' in tt.get('attrspec', []):
      return 'f90_' + tt['typename'] + '_pointer'
    else:
      return 'f90_' + tt['typename']
  raise ValueError('Unknown Fortran type %s.' % tt['typespec'])

def toCReturnType(tt):
  if tt['typespec'] == 'type' and 'pointer' not in tt.get('attrspec', []):
    return toCType(tt) + '_pointer'
  else:
    return toCType(tt)

def toFortranType(tt, lbl):
  if tt['typespec'] == 'type':
    out = 'type(%s)' % tt['typename']
  elif tt['typespec'] == 'character':
    out = 'character(len = %s)' % (tt['charselector'].get('*', tt['charselector']['len'] if 'len' in tt['charselector'] else '') if isFixedSizeCharacter(tt) else '%s_len' % lbl)
  else:
    out = tt['typespec']
    if 'kindselector' in tt:
      out += '(kind = %s)' % tt['kindselector']['kind']
  if 'dimension' in tt:
    out += ', dimension('
    out += tt['dimension'][0] if tt['dimension'][0] != ':' else ('%s_dim_0' % lbl)
    for i, d in enumerate(tt['dimension'][1:]):
      out += ', '
      out += d if d != ':' else ('%s_dim_%d' % (lbl, i + 1))
    out += ')'
  return out

def isFixedSizeArray(tt):
  return (':' not in tt['dimension'] and '*' not in tt['dimension'])

def isFixedSizeCharacter(tt):
  return tt['charselector'].get('*', ('(' + tt['charselector']['len'] + ')') if 'len' in tt['charselector'] else '') != '(*)'

def isFixedSizeTypeArray(tt):
  return 'dimension' in tt and isFixedSizeArray(tt) and tt['typespec'] == 'type' and 'intent' in tt

def isCTypeAsPointer(tt):
  return ('dimension' in tt and not isFixedSizeArray(tt)) \
    or tt['typespec'] == 'type' \
    or (tt['typespec'] == 'character' and not isFixedSizeCharacter(tt))

class BindArg(object):
  def __init__(self, name, aDef):
    self.block = aDef
    self.name = name
    
    if aDef.get('typespec', '') == 'type' and aDef['typename'] == 'f_enumerator' and \
      ('out' in aDef.get('intent', []) or 'inout' in aDef.get('intent', [])):
      raise NotImplementedError('f_enumerator as return type.')

    self.ctype = aDef['ctype'] if 'ctype' in aDef else toCType(aDef)
    if 'in' in aDef.get('intent', ['inout', ]):
      self.ctype = 'const ' + self.ctype

    self.byValue = aDef.get('byValue', False)
      
    isOptional = 'optional' in aDef.get('attrspec', [])
    isFixedSize = 'dimension' in aDef and isFixedSizeArray(aDef)
    isPointer = 'typespec' in aDef and isCTypeAsPointer(aDef)
    if isOptional and 'dimension' in aDef:
      first = '[0]' * len(aDef['dimension'])
      self.call = '%s ? &(*%s)%s : NULL' % (name, name, first)
    elif isOptional and aDef['typespec'] == 'character':
      self.call = name
    elif isPointer and isFixedSize:
      self.call = '*%s' % name
    elif isFixedSize:
      first = '[0]' * len(aDef['dimension'])
      self.call = '&%s%s' % (name, first)
    elif self.byValue or isOptional or isPointer:
      self.call = name
    else:
      self.call = '&%s' % name
    # Special case where we look at a specific index from an array
    self.index = aDef['index'] if 'index' in aDef else None
    if self.index is not None:
      self.call += '[%d]' % self.index

    self.ftype = toFortranType(aDef, name) if 'typespec' in aDef else None
    
  def toDecl(self):
    return self.ctype + ('*' if not self.byValue else '')

  def toCall(self):
    return self.call

  def toFDecl(self):
    if self.index is not None:
      return '%s_%d' % (self.name, self.index)
    else:
      return self.name

  def toFDefinition(self):
    if self.ftype is None:
      return ''
    out = self.ftype
    for a in self.block.get('attrspec', []):
      out += ', %s' % a
    if 'intent' in self.block:
      out += ', intent(%s)' % self.block['intent'][0]
    if self.index is not None:
      out += ' :: %s_%d\n' % (self.name, self.index)
    else:
      out += ' :: %s\n' % self.name
    return out

class BindReturnArg(BindArg):
  def __init__(self, name, aDef):
    tt = aDef.copy()
    if tt['typespec'] == 'type' and 'pointer' not in tt.get('attrspec', []):
      # Cannot return a type, so work on pointers instead
      if 'attrspec' in tt:
        tt['attrspec'].append('pointer')
      else:
        tt['attrspec'] = ['pointer', ]
    if 'pointer' not in tt.get('attrspec', []):
      tt['intent'] = ['out', ]
    if tt['typespec'] == 'type' and tt['typename'] == 'f_enumerator' and \
      not tt.get('constructor', False):
      raise NotImplementedError('f_enumerator as return type.')

    super(BindReturnArg, self).__init__(name, tt)

  def toType(self, ctype):
    if self.ctype == ctype:
      return self.name
    if ctype.endswith('*') and self.ctype == ctype[:-1] + '_pointer':
      return 'F_TYPE(%s)' % self.name
    raise ValueError('Cannot convert from %s to %s' % (self.ctype, ctype))

  def toInit(self):
    out = ''
    if 'pointer' in self.block.get('attrspec', []):
      out += '  nullify(%s)\n' % self.name
      if self.block['typespec'] == 'type' and self.block['typename'] != 'f_enumerator':
        out += '  allocate(%s)\n' % self.name
    return out

  def toCall(self):
    out = super(BindReturnArg, self).toCall()
    if out[0] == '&':
      return out
    else:
      return '&' + out

  def toValue(self):
    if self.block['typespec'] == 'type' and self.block['typename'] == 'f_enumerator':
      return '%s =>' % self.name
    else:
      return '%s =' % self.name
    
class CArg(object):
  def __init__(self, name, aDef):
    if 'external' in aDef.get('attrspec', []):
      raise NotImplementedError('external argument')
    isPointer = isCTypeAsPointer(aDef)
    isFixedSize = ('dimension' in aDef and isFixedSizeArray(aDef))

    self.block = aDef
    
    self.name = name
    self.optional = 'optional' in aDef.get('attrspec', [])
    self.isString = aDef['typespec'] == 'character' and not isFixedSizeCharacter(aDef)
    self.isAssumeSize = 'dimension' in aDef and not isFixedSizeArray(aDef)

    self.ctype = aDef['ctype'] if 'ctype' in aDef else toCType(aDef)
    if 'in' in aDef.get('intent', ['inout', ]) and (isPointer or isFixedSize):
      self.ctype = 'const ' + self.ctype
    if isPointer and not(self.optional and self.isString):
      self.ctype += '*'

    self.dims = ''
    if isFixedSize:
      for dim in aDef['dimension']:
        self.dims += '[%s]' % dim
      if aDef['typespec'] == 'character' and isFixedSizeCharacter(aDef):
        ln = aDef['charselector'].get('*', aDef['charselector']['len'])
        if ln != 1:
          self.dims += '[%s]' % ln

  def toDecl(self):
    lbl = ('(*' + self.name + ')') if self.optional else self.name
    return self.ctype + ' ' + lbl + self.dims

class CReturnArg(CArg):
  def __init__(self, name, aDef):
    if 'dimension' in aDef and not isFixedSizeArray(aDef):
      raise NotImplementedError('dynamic array as a return value.')

    super(CReturnArg, self).__init__(name, aDef)

    self.ctype = toCType(aDef)
    if aDef['typespec'] == 'type' and 'pointer' not in aDef.get('attrspec', []):
      if aDef['typename'] == 'f_enumerator':
        self.ctype += '*'
      else:
        self.ctype += '_pointer'

    self.isArray = 'dimension' in aDef and isFixedSizeArray(aDef)

class BindFunc(object):
  def __init__(self, block, generic = None, fbody = None):
    self.block = block
    self.generic = generic if generic is not None else block['name']
    self.fbody = fbody
    self.args = []
    self.retArg = None
    self.strArgs = []
    self.isValid = False
    try:
      self.__analyseBlock()
      self.isValid = True
    except NotImplementedError as e:
      print '%s error:' % block['name'], e

  def __analyseBlock(self):
    # List of C / Fortran arguments
    for a in self.block['args']:
      arg = CArg(a, self.block['vars'][a])
      # Prefix a size argument of every assume size array
      if arg.isAssumeSize:
        for i, s in enumerate(arg.block['dimension']):
          self.args.append(CArg('%s_dim_%d' % (a, i), {'typespec': 'integer', 'intent': ['in', ], 'ctype': 'size_t'}))
      self.args.append(arg)
    # Return argument if any
    retName = self.block.get('result', self.block['name']) if self.block['block'] == 'function' else None
    if retName is not None:
      self.retArg = CReturnArg(retName, self.block['vars'][retName])
      if self.retArg.isArray:
        # Convert returned array as out argument.
        tt = self.block['vars'][retName].copy()
        tt['intent'] = ['out', ]
        self.args.append(CReturnArg('arg_%s' % retName, tt))
        self.retArg = None
    # List of argument names that are Fortran character
    self.strArgs = [a for a in self.block['args'] if self.block['vars'][a]['typespec'] == 'character']
    # List of C / Fortran arguments for the binding routine
    self.bindArgs = []
    self.expandedArgs = []
    if self.retArg is not None:
      # Transform retArg into an out argument of the bind subroutine
      self.bindArgs.append(BindReturnArg('out_%s' % self.retArg.name, self.retArg.block))
    for a in self.args:
      if isFixedSizeTypeArray(a.block):
        self.expandedArgs.append(a)
        tt = a.block.copy()
        tt.pop('dimension')
        # Expanded fixed-size type array into references on type
        for i in range(int(a.block['dimension'][0])):
          tt['index'] = i
          self.bindArgs.append(BindArg(a.name, tt))
      else:
        if a.isString:
          # Prefix string length
          self.bindArgs.append(BindArg('%s_len' % a.name, {'typespec': 'integer', 'intent': ['in', ], 'ctype': 'size_t'}))
        self.bindArgs.append(BindArg(a.name, a.block))
    for a in self.strArgs:
      # GFortran is adding a by-value size for bounds-check of strings
      self.bindArgs.append(BindArg('%s_len' % a, {'ctype': 'size_t', 'byValue': True}))

  def __toDecl(self, out):
    if self.retArg is not None:
      out.write('%s %s' % (self.retArg.ctype, self.block['name']))
    else:
      out.write('void %s' % self.block['name'])
    out.write('(')
    if len(self.args) > 0:
      out.write(self.args[0].toDecl())
    else:
      out.write('void')
    for a in self.args[1:]:
      out.write(",\n  ")
      out.write(a.toDecl())
    out.write(')')

  def __toFCFunc(self):
    return 'FC_FUNC%s(bind_%s, BIND_%s)' % ("_" if "_" in self.block['name'] else "", self.block['name'].lower(), self.block['name'].upper())

  def __toBindCDecl(self, out):
    out.write('(')
    if len(self.bindArgs) > 0:
      out.write(self.bindArgs[0].toDecl())
    else:
      out.write('void')
    for a in self.bindArgs[1:]:
      out.write(',\n  ')
      out.write(a.toDecl())
    out.write(')')

  def __callBindCFunc(self, out):
    out.write('  %s\n' % self.__toFCFunc())
    out.write('    (')
    if len(self.bindArgs) > 0:
      out.write(self.bindArgs[0].toCall())
    for a in self.bindArgs[1:]:
      out.write(', ')
      out.write(a.toCall())
    out.write(');\n')
    
  def toHeader(self, out):
    if not self.isValid:
      out.write('/* ')
    self.__toDecl(out)
    out.write(';')
    if not self.isValid:
      out.write(' */')
    out.write('\n\n')

  def toImpl(self, out):
    if not self.isValid:
      return
    out.write('void %s' % self.__toFCFunc())
    self.__toBindCDecl(out)
    out.write(';\n')
    
    self.__toDecl(out)
    out.write('\n{\n')

    # Add some local variables, if needed
    for a in self.strArgs:
      if isFixedSizeCharacter(self.block['vars'][a]):
        out.write('  size_t %s_len = %s;\n' %(a, self.block['vars'][a]['charselector'].get('*', self.block['vars'][a]['charselector']['len'])))
      else:
        out.write('  size_t %s_len = %s ? strlen(%s) : 0;\n' %(a, a, a))
    if len(self.bindArgs) > 0 and isinstance(self.bindArgs[0], BindReturnArg):
      out.write('  %s %s;\n' % (self.bindArgs[0].ctype, self.bindArgs[0].name))
    
    self.__callBindCFunc(out)
    if len(self.bindArgs) > 0 and isinstance(self.bindArgs[0], BindReturnArg):
      out.write('  return %s;\n' % self.bindArgs[0].toType(self.retArg.ctype))
    out.write('}\n\n')

  def toBind(self, out, m, uses):
    if not self.isValid:
      return
    out.write('subroutine bind_%s( &\n    ' % self.block['name'])
    if len(self.bindArgs) > 0:
      out.write(self.bindArgs[0].toFDecl())
    for a in self.bindArgs[1:]:
      if a.ftype is not None:
        out.write(', &\n    ')
        out.write(a.toFDecl())
    out.write(')\n')
    uses_ = self.block.get('use', {}).copy()
    uses_.update(uses)
    for (u, a) in uses_.items():
      out.write('  use %s' % u)
      if 'map' in a:
        out.write(', only: ')
        it = a['map'].items()
        if len(it) > 0:
          out.write('%s => %s' % (it[0][0], it[0][1]))
        for fr, to in it[1:]:
          out.write(', &\n    %s => %s' % (fr, to))
      out.write('\n')
    out.write('  use %s\n' % m)
    out.write('  implicit none\n')
    for a in [a for a in self.bindArgs if a.ftype is not None]:
      out.write('  ' + a.toFDefinition())
    for a in self.expandedArgs:
      out.write('  %s :: %s\n' % (toFortranType(a.block, ''), a.name))
    out.write('\n')
    for a in [a for a in self.bindArgs if a.index is not None \
                and ('in' in a.block.get('intent', ['inout', ]) or \
                       'inout' in a.block.get('intent', ['inout', ]))]:
      out.write('  %s(%d) = %s_%d\n' % (a.name, a.index + 1, a.name, a.index))
    if len(self.bindArgs) > 0 and isinstance(self.bindArgs[0], BindReturnArg) and 'pointer' in self.bindArgs[0].block.get('attrspec', []):
      out.write(self.bindArgs[0].toInit())
    if self.fbody is None:
      end = -1
      if len(self.bindArgs) > 0 and isinstance(self.bindArgs[0], BindReturnArg):
        out.write('  %s %s( &\n' % (self.bindArgs[0].toValue(), self.generic))
      elif len(self.args) > 0 and isinstance(self.args[-1], CReturnArg):
        out.write('  %s = %s( &\n' % (self.args[-1].name, self.generic))
        end = -2
      else:
        out.write('  call %s( &\n' % self.generic)
      for a in self.block['args'][:end]:
        out.write('    %s, &\n' % a)
      if len(self.block['args']) > 0:
        out.write('    %s' % self.block['args'][end])
      out.write(')\n')
    else:
      out.write(self.fbody)
    for a in [a for a in self.bindArgs if a.index is not None \
                and ('out' in a.block.get('intent', ['inout', ]) or \
                       'inout' in a.block.get('intent', ['inout', ]))]:
      out.write('  %s_%d = %s(%d)\n' % (a.name, a.index, a.name, a.index + 1))
    out.write('end subroutine bind_%s\n\n' % self.block['name'])
    
class BindType(object):
  def __init__(self, aDef):
    self.name = aDef['name']
    tt = {'block': 'function', 'name': 'f90_%s_new' % self.name, \
            'args': [], \
            'result': 'self', \
            'vars': {'self': \
                       {'typespec': 'type', 'typename': self.name}} \
            }
    self.newFunc = BindFunc(tt, fbody = '')
    tt = {'block': 'subroutine', 'name': 'f90_%s_free' % self.name, \
            'args': ['self', ], \
            'vars': {'self': \
                       {'typespec': 'type', 'typename': self.name, 'attrspec': ['pointer', ]}} \
            }
    self.freeFunc = BindFunc(tt, fbody = '  deallocate(self)\n  nullify(self)\n')

  def toDecl(self, out):
    out.write('F_DEFINE_TYPE(%s);\n' % self.name)
    self.newFunc.toHeader(out)
    self.freeFunc.toHeader(out)

  def toImpl(self, out):
    self.newFunc.toImpl(out)
    self.freeFunc.toImpl(out)

  def toBind(self, out, m):
    self.newFunc.toBind(out, m, {})
    self.freeFunc.toBind(out, m, {})

tt = crackfortran.crackfortran(sys.argv[1])[0]

if tt['block'] != 'module':
  raise TypeError('Not a module file')
modname = tt['name']

def publicBlocks(tt, kind):
  interfaces = []
  if kind != 'interface':
    # Some funny people like to name some routines like the interface itself
    interfaces = [i['name'] for i in publicBlocks(tt, 'interface')]
  return [b for b in tt['body'] if b['block'] == kind and b['name'] in tt['vars'] and 'public' in tt['vars'][b['name']]['attrspec'] and b['name'] not in interfaces]

def publicInterfaces(tt):
  return [(b['name'], i) for b in publicBlocks(tt, 'interface') if not b['name'].startswith('operator') for i in tt['body'] if i['name'] in b['implementedby'] and (i['block'] == 'function' or i['block'] == 'subroutine')]

def publicEnums(tt):
  return [(k, b) for k, b in tt['vars'].items() if 'typename' in b and b['typename'] == 'f_enumerator' and 'public' in b['attrspec']]

def publicUsedTypes(tt):
  local = [b['name'] for b in publicBlocks(tt, 'type')]
  return {v['typename'] for t in publicBlocks(tt, 'type') for v in t['vars'].values() if v['typespec'] == 'type' and v['typename'] not in local} \
    | {v['typename'] for f in publicBlocks(tt, 'function') for a in f['args'] for k, v in f['vars'].items() if k == a and v['typespec'] == 'type' and v['typename'] not in local} \
    | {v['typename'] for f in publicBlocks(tt, 'subroutine') for a in f['args'] for k, v in f['vars'].items() if k == a and v.get('typespec', []) == 'type' and v['typename'] not in local} \
    | {v['typename'] for i, f in publicInterfaces(tt) for a in f['args'] for k, v in f['vars'].items() if k == a and v.get('typespec', []) == 'type' and v['typename'] not in local}

types = [BindType(b) for b in publicBlocks(tt, 'type')]

funcs = [BindFunc(b) for b in publicBlocks(tt, 'function')] \
  + [BindFunc(b) for b in publicBlocks(tt, 'subroutine')] \
  + [BindFunc(b, i) for (i, b) in publicInterfaces(tt)]
for k, b in publicEnums(tt):
  aDef = {'block': 'function', 'name': k, \
            'args': [], \
            'result': 'self', \
            'vars': {'self': b.copy()}}
  aDef['vars']['self'].pop('attrspec')
  aDef['vars']['self']['constructor'] = True
  funcs.append(BindFunc(aDef, fbody = '  out_self => %s\n' % k))

with open(tt['name'] + '.h', 'w') as header:
  header.write('#ifndef %s_H\n' % tt['name'].upper())
  header.write('#define %s_H\n' % tt['name'].upper())
  header.write('\n')

  header.write('#include <futile/futile_cst.h>\n')
  for t in publicUsedTypes(tt):
    if t == 'dictionary':
      header.write('#include <futile/dict.h>\n')
    elif t == 'f_enumerator':
      header.write('#include <futile/enum.h>\n')
    else:
      header.write('#include "%s.h"\n' % t)
  header.write('\n')

  for b in types:
    b.toDecl(header)
  header.write('\n')

  for b in funcs:
    b.toHeader(header)

  header.write('#endif\n')

with open(tt['name'] + '_c.c', 'w') as cimpl:
  cimpl.write('#include "%s.h"\n' % tt['name'])
  cimpl.write('#include <config.h>\n')
  cimpl.write('#include <string.h>\n\n')

  for b in types:
    b.toImpl(cimpl)

  for b in funcs:
    b.toImpl(cimpl)

with open(tt['name'] + '_bind.f90', 'w') as fimpl:
  for b in types:
    b.toBind(fimpl, tt['name'])

  for b in funcs:
    b.toBind(fimpl, tt['name'], tt.get('use', {}))
