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


from code import  sagematrix_to_func

def _gen_q_dq_ddq_subs(dof):
    subs = []
    for i in reversed(range(dof)):
        subs.append( ( 'ddq'+str(i+1)  ,'ddq['+str(i)+']' ) )
        subs.append( ( 'dq'+str(i+1)  ,'dq['+str(i)+']' ) )
        subs.append( ( 'q'+str(i+1)  ,'q['+str(i)+']' ) )
    return subs

def _gen_parms_subs( parms_symbols, name = 'parms' ):
    subs = []
    for i in reversed(range(len(parms_symbols))):
        subs.append( ( str(parms_symbols[i]), name+'['+str(i)+']' ) )
    return subs

 
def dyn_matrix_to_func( lang, sageauxvars, sagematrix, funcname, qderivlevel, dof, dynparam_symbols=[] ):
    func_parms = []
    subs_pairs = []
    if dynparam_symbols:
        func_parms.append('parms')
        subs_pairs += _gen_parms_subs(dynparam_symbols,'parms')
    if qderivlevel >= 0:
        for i in range(qderivlevel+1):
            func_parms.append('d'*i+'q')
        subs_pairs += _gen_q_dq_ddq_subs(dof)
    return sagematrix_to_func(  lang, sageauxvars, sagematrix, funcname, func_parms, subs_pairs )


