# -*- coding: utf-8 -*-

###############################################################################
#  Sage Robotics: Robotics Toolbox for Sage Mathematical Software
#
#      Copyright (C) 2010, 2011 Cristóvão Sousa <crisjss@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################



import sympy
from copy import deepcopy, copy


def sagematrix_to_code( auxvars, matrix ):
    code = [[],[]]
    for a in auxvars:
        code[0].append( ( sympy.sympify(a[0]._sympy_()) , sympy.sympify(a[1]._sympy_()) ) )
    for t in matrix.list():
        code[1].append( sympy.sympify(t._sympy_()) )
    return code



def code_cse( code, auxvarname = 'ccse' ):
    
    code_in = code

    import sympy.simplify.cse_main
    codecse = sympy.cse( sympy.sympify( [i[1] for i in code_in[0]] + code_in[1] ), \
                         sympy.cse_main.numbered_symbols(auxvarname) )
    
    auxv1_num = len(code_in[0])
    
    if auxv1_num == 0:
        return codecse
    
    A1 = zip( zip(*code_in[0])[0] , codecse[1][:auxv1_num] ) 
    Acse = codecse[0]
    
    def getdeps(expr,assignments):
        out = []
        search_assignments = copy(assignments)
        for assign in search_assignments:
            if sympy.sympify(expr).has(assign[0]) and assign in assignments:
                assignments.remove(assign)
                out += getdeps(assign[1],assignments) + [assign]
        return out
    
    codemerge = []
    for a1 in A1:
        codemerge += getdeps(a1[1],Acse)
        codemerge.append( a1 )
    for acse in Acse:
        codemerge.append( acse )
    
    retcode = [ codemerge , codecse[1][auxv1_num:] ]
    
    return retcode


def code_remove_not_compound( code ):
    
    debug = False
    removed=0
    
    retcode = deepcopy(code)
    toremove = []
    
    for i in range(len(retcode[0])):
        v = retcode[0][i][0]
        e = retcode[0][i][1]
        if e.is_Atom or (-e).is_Atom:
            if debug: print i,v,e,'is atom'
            for j in range(i+1,len(retcode[0])):
                retcode[0][j] = ( retcode[0][j][0] , retcode[0][j][1].subs({v:e}) )
            for j in range(len(retcode[1])):
                retcode[1][j] = retcode[1][j].subs({v:e})
            toremove.append(i)
            removed += 1
            if debug: print '  poped'

    for r in reversed(toremove):
        retcode[0].pop(r)
    
    if debug: print 'removed',removed
    return retcode


def code_remove_not_or_once_used( code ):
    
    debug = False
    removed=0
    
    retcode = deepcopy(code)
    
    for i in range(len(retcode[0])-1,-1,-1):
        v = retcode[0][i][0]
        e = retcode[0][i][1]
        uses = 0
        usedin_ai = None
        usedin_oi = None
        if debug: print i,v
        for ai in range(i+1,len(retcode[0])):
            if sympy.sympify(retcode[0][ai][1]).has(v):
                uses += 1
                usedin_ai = ai
                # if debug: print '  used in',ai
        for oi in range(len(retcode[1])):
            if sympy.sympify(retcode[1][oi]).has(v):
                uses += 1
                usedin_oi = oi
                # if debug: print '  used in out',oi
        
        if uses == 0:
            retcode[0].pop(i)
            removed += 1
            if debug: print '  not used: removed'
        elif uses == 1:
            if usedin_ai != None:
                retcode[0][usedin_ai] = ( retcode[0][usedin_ai][0] , retcode[0][usedin_ai][1].subs({v:e}) )
                if debug: print '  used once: removed and substituted in',retcode[0][usedin_ai][0]
            elif usedin_oi != None:
                retcode[1][usedin_oi] = retcode[1][usedin_oi].subs({v:e})
                if debug: print '  used once: removed and substituted in out',usedin_oi
            retcode[0].pop(i)
            removed += 1
    
    if debug: print 'removed',removed
    return retcode


def code_rename_auxvars(code, auxvnames='auxA', auxv_sagelatexnames=None ):
    
    retcode = deepcopy(code)
    
    for i in range(len(retcode[0])):
        
        Dsubs = { retcode[0][i][0] : sympy.var(auxvnames+str(i),real=True) }
        
        if auxv_sagelatexnames:
            from sage.all import var as sage_var, RR as sage_RR
            sage_var(auxvnames+str(i),latex_name=auxv_sagelatexnames+'_{'+str(i)+'}',domain=sage_RR)
        
        retcode[0][i] = ( retcode[0][i][0].subs( Dsubs ) , retcode[0][i][1] )
        
        for j in range(i+1,len(retcode[0])):
            retcode[0][j] = ( retcode[0][j][0], retcode[0][j][1].subs( Dsubs ) )
        
        for j in range(len(retcode[1])):
            retcode[1][j] = sympy.sympify(retcode[1][j]).subs( Dsubs )
    
    return retcode


def optimize_code(code) :
    retcode = deepcopy(code)
    retcode = code_remove_not_compound(retcode)
    retcode = code_remove_not_or_once_used(retcode)
    retcode = code_cse(retcode,'cse')
    retcode = code_remove_not_compound(retcode)
    retcode = code_remove_not_or_once_used(retcode)
    retcode = code_rename_auxvars(retcode, auxvnames='aux', auxv_sagelatexnames=r'\mathbf{A}')
    return retcode


def code_subs( code, subs_dict ):
    retcode = deepcopy(code)
    for i in range(len(retcode[0])):
        retcode[0][i] = ( retcode[0][i][0] , retcode[0][i][1].subs(subs_dict) )
    for i in range(len(retcode[1])):
        retcode[1][i] = retcode[1][i].subs(subs_dict)
    return retcode


def code_back_to_expressions(code):
    
    avars = deepcopy(code[0])
    exps = deepcopy(code[1])
    
    for i in range(len(avars)):
        
        Dsubs = { avars[i][0] : avars[i][1] }
        
        for j in range(i+1,len(avars)):
                if avars[j][1].has(avars[i][0]):
                    avars[j] = ( avars[j][0], avars[j][1].subs( Dsubs ) )
        
        for j in range(len(exps)):
                if sympy.sympify(exps[j]).has(avars[i][0]):
                    exps[j] = exps[j].subs( Dsubs )
    
    return exps


def code_to_string( code, outvarname='out', indent='', realtype='', line_end='' ):
    
    codestr = ''
    
    if realtype: realtype += ' '
    
    import sympy.printing
    
    for i in range( len(code[0]) ) :
        codestr += indent + realtype + sympy.ccode( code[0][i][0] ) + ' = ' + sympy.ccode( code[0][i][1] ) + line_end + '\n'
    
    codestr += '\n'
    for i in range( len(code[1]) ) :
        codestr += indent + outvarname + '['+str(i)+'] = ' + sympy.ccode( code[1][i] ) + line_end + '\n'
    
    return codestr


def codestring_count( codestring, resume=False ):
  ops = []
  ops += [( '=' , int(codestring.count('=')) )]
  ops += [( '+' , int(codestring.count('+')) )]
  ops += [( '-' , int(codestring.count('-')) )]
  ops += [( '*' , int(codestring.count('*')) )]
  ops += [( '/' , int(codestring.count('/')) )]
  ops += [( 'pow' , int(codestring.count('pow')) )]
  ops += [( 'sin' , int(codestring.count('sin')) )]
  ops += [( 'cos' , int(codestring.count('cos')) )]
  if not resume:
    return ops
  else:
    adds = int(codestring.count('+'))+int(codestring.count('-'))
    muls = int(codestring.count('*'))+int(codestring.count('/'))+int(codestring.count('pow'))
    return ops, {'add':adds, 'mul':muls, 'total':adds+muls }



def gen_py_func( code, func_parms, subs_pairs, funcname='func', outvarname='out' ):
    
    indent = 4*' '
    
    pycode = 'def ' + funcname + '('
    if func_parms:
        pycode += ' ' + func_parms[0]
        for parm in func_parms[1:] :
            pycode += ', ' + parm
        pycode += ' '
    pycode += ') :\n\n' 
    
    pycode += indent + outvarname + ' = [0]*' + str( len(code[1]) ) + '\n\n'
    
    mainpycode = code_to_string( code, outvarname, indent )
    for (old,new) in subs_pairs: mainpycode = mainpycode.replace(old,new)
    pycode += mainpycode
    
    pycode += '\n' + indent + 'return ' + outvarname + '\n'  
    return pycode


def gen_c_func( code, func_parms, subs_pairs, func_name='func', outvar_name='out' ):
    
    indent = 2*' ' 
    
    ccode = 'void ' + func_name + '( double* ' + outvar_name
    for parm in func_parms :
        ccode += ', const double* ' + parm
    ccode += ' )\n{\n'
    
    mainccode = code_to_string( code, outvar_name, indent, 'double', ';' )
    for (old,new) in subs_pairs: mainccode = mainccode.replace(old,new)
    
    ccode += mainccode + '}\n'
    return ccode



def sagematrix_to_func( lang, sageauxvars, sagematrix, func_name, func_parms, subs_pairs ):
    lang = lang.lower()
    if lang in ['python','py'] : gen_func = gen_py_func
    elif lang in ['c','c++'] : gen_func = gen_c_func
    else: raise Exception('chosen language not supported.')
    code = optimize_code( sagematrix_to_code( sageauxvars, sagematrix ) )
    return gen_func( code, func_parms, subs_pairs, func_name, func_name+'_out' )


