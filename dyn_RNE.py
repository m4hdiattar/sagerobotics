# -*- coding: utf-8 -*-

###############################################################################
#  Sage Robotics: Robotics Toolbox for Sage Mathematical Software
#
#      Copyright (C) 2010 Cristóvão Sousa <crisjss@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################




from utils import *
from sage.all import ( var, restore, I, e, SR, RR, Integer, exp, cos, sin,
                       matrix, block_matrix, identity_matrix, zero_matrix )


def _gen_rbt_Si( rbt ):
    
    var('a theta alpha d',domain=RR)
    D_exp2trig = { e**(SR.wild(0)) : e**SR.wild(0).real_part() * ( cos(SR.wild(0).imag_part()) + I*sin(SR.wild(0).imag_part()) ) }
    M = matrix(SR,[[1,0,0,a],[0,cos(alpha),-sin(alpha),0],[0,sin(alpha),cos(alpha),d],[0,0,0,1]])
    P = matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])
    ePtM = (exp(P*theta)*M).expand().subs(D_exp2trig).simplify_rational()
    restore('a theta alpha d')
    if bool( ePtM != rbt.T_dh ):
        raise Exception('_gen_rbt_Si: interlink transformation does not follows implemented DH formulation')
    S = ( inverse_T(M) * P * M ).expand().simplify_trig()
    
    dof = rbt.dof
    
    Si = range(dof+1)
    for l in range(dof): Si[l+1] = se3unskew( S.subs(LoP_to_D(zip(rbt.params_dh , rbt.links_dh[l]))) )
    
    return Si


def _forward_RNE( rbt, Si, varify_func=None ):
  
    from copy import copy
    
    dof = rbt.dof
    
    Vi = range(0,dof+1)
    dVi = range(0,dof+1)
    
    Vi[0] = copy(zero_matrix(SR,6,1))
    dVi[0] = - zero_matrix(SR,3,1).stack( rbt.grav )
    
    # Forward
    for i in range(1,dof+1):
        Vi[i] =  Adj( rbt.Tdhi_inv[i], Vi[i-1] )  +  Si[i] * rbt.dq[i-1,0]
        if varify_func: Vi[i] = varify_func( Vi[i] , 'V_'+str(i) )
        
        dVi[i] =  Si[i] * rbt.ddq[i-1,0]  +  Adj( rbt.Tdhi_inv[i], dVi[i-1] )  +  adj(  Adj( rbt.Tdhi_inv[i], Vi[i-1] ),  Si[i] * rbt.dq[i-1,0] )
        if varify_func: dVi[i] = varify_func( dVi[i] , 'dV_'+str(i) )

    return Vi, dVi


def _backward_RNE( rbt, IIi, Si, Vi, dVi, usemotordyn=True, Imzi=None, varify_func=None ) :
  
    from copy import copy
    
    dof = rbt.dof
    Tdhi_inv = copy(rbt.Tdhi_inv)
    Tdhi_inv.append(identity_matrix(4))
    
    Fi = range(0,dof+2)
    tau = matrix(SR,dof,1)
    
    Fi[dof+1] = copy(zero_matrix(SR,6,1))

    if usemotordyn:
      nfi_motor = [zero_matrix(3,1) for i in xrange(dof+1)]
      for mi,m in enumerate(rbt.motors):
        if Imzi == None: Im = m[0]
        else: Im = Imzi[mi]
        qm = m[1]; l = m[2]; zm = m[3]
        dqm = qm.subs(rbt.D_q_v2f).derivative(rbt.t).subs(rbt.D_q_f2v)
        ddqm = dqm.subs(rbt.D_q_v2f).derivative(rbt.t).subs(rbt.D_q_f2v)
        nfi_motor[l] += ddqm * Im * zm + dqm * Im * cross_product( Vi[l][:3,:] , zm )
    
    # Backward
    for i in range(dof,0,-1):
        Fi[i] =  Adjdual( Tdhi_inv[i+1], Fi[i+1] )  +  IIi[i] * dVi[i]  -  adjdual( Vi[i],  IIi[i] * Vi[i] )

        if usemotordyn:
          Fi[i][:3,:] += nfi_motor[i]
        
        if varify_func: Fi[i] = varify_func( Fi[i] , 'F_'+str(i) )
        
        tau[i-1,0] =  ( Si[i].transpose() *  Fi[i] )[0,0]

        if usemotordyn:
          for mi,m in enumerate(rbt.motors):
            if Imzi == None: Im = m[0]
            else: Im = Imzi[mi]
            qm = m[1]; l = m[2]; zm = m[3]
            dqm = qm.subs(rbt.D_q_v2f).derivative(rbt.t).subs(rbt.D_q_f2v)
            ddqm = dqm.subs(rbt.D_q_v2f).derivative(rbt.t).subs(rbt.D_q_f2v)
            km = qm.coeff(rbt.q[i-1,0])
            if km != 0.0:
              tau[i-1,0] += km * Im * ( ( dVi[l][:3,:] + ddqm * zm + dqm * cross_product( Vi[l][:3,:], zm ) ).transpose() * zm )[0,0]
    
    return tau



def gen_tau_RNE( do_varify, rbt, usemotordyn=True, varify_trig = True ):
    
    Si = _gen_rbt_Si(rbt)
    
    IIi = range(0,rbt.dof+1)
    ### Linear dynamic parameters:
    for i in range(1,rbt.dof+1): IIi[i] = block_matrix( [ rbt.Ifi[i] , skew(rbt.mli[i]) , -skew(rbt.mli[i]) , rbt.mi[i]*identity_matrix(3) ] ,subdivide=False)
    ### Not linear dynamic parameters:
    #for i in range(1,dof+1): IIi[i] = block_matrix( [ rbt.Ifi_from_Ici[i] , skew(rbt.mli_e[i]) , -skew(rbt.mli_e[i]) , rbt.mi[i]*identity_matrix(3) ] ,subdivide=False)
    
    if do_varify:
        auxvars = []
        def custom_simplify( expression ):
        #    return trig_reduce(expression.expand()).simplify_rational()
            return expression.simplify_rational()
        
        def varify_func(exp,varrepr):
            if varify_trig : exp = exp.subs( LoP_to_D(rbt.LoP_trig_f2v) )
            exp = custom_simplify( exp )
            return m_varify(exp,auxvars,condition_func=is_compound)
            #return v6R_varify(exp,auxvars,varrepr)
        
    else:
        def varify_func(exp,varrepr):
            return exp
    
    Vi, dVi = _forward_RNE( rbt, Si, varify_func )
    tau = _backward_RNE( rbt, IIi, Si, Vi, dVi, usemotordyn, None, varify_func = varify_func )
    
    if do_varify:
        if varify_trig : auxvars = rbt.LoP_trig_v2f + auxvars
        return auxvars,tau
    else:
        return tau




def gen_regressor_RNE( do_varify, rbt, usemotordyn = True, usefricdyn = True, varify_trig = True ):
  
    from copy import copy
    
    Si = _gen_rbt_Si(rbt)

    if usefricdyn:
      fric = gen_Dyn_fricterm(rbt)
    
    def custom_simplify( expression ):
        #return trig_reduce(expression.expand()).simplify_rational()
        return expression.simplify_rational()
    
    varify_func = None
    
    if do_varify:
        auxvars = []
        def varify_func(exp,varrepr):
            if varify_trig : exp = exp.subs( LoP_to_D(rbt.LoP_trig_f2v) )
            exp = custom_simplify( exp )
            return m_varify(exp, auxvars, poolrepr='auxY', condition_func=is_compound)
    
    
    Vi, dVi = _forward_RNE( rbt, Si, varify_func )
    
    P = rbt.Parms(usemotordyn,usefricdyn)
    
    Y = matrix(SR,rbt.dof,P.nrows())
    
    for p in range(P.nrows()) :
        
        select =  subsm( P, zero_matrix(P.nrows(),1) )
        select.update( {P[p,0]:Integer(1)} )
        
        IIi = range(0,rbt.dof+1)
        for i in range(1,rbt.dof+1):
            IIi[i] = block_matrix( [ rbt.Ifi[i].subs(select) , skew(rbt.mli[i].subs(select)) , -skew(rbt.mli[i].subs(select)) , rbt.mi[i].subs(select)*identity_matrix(3) ] ,subdivide=False)
        
        if usemotordyn:
          Imzi = deepcopy(rbt.Imzi)
          for i,Im in enumerate(Imzi):
            Imzi[i] = Im.subs(select)
        else:
          Imzi = None
      
        Y[:,p] = _backward_RNE( rbt, IIi, Si, Vi, dVi, usemotordyn, Imzi, varify_func )
        
        if usefricdyn:
          Y[:,p] += fric.subs(select)
    
    if do_varify:
        if varify_trig : auxvars = rbt.LoP_trig_v2f + auxvars
        return auxvars , Y
    else:
        return Y



def gen_gravterm_RNE( do_varify, rbt, varify_trig=True ):
    from copy import deepcopy
    rbttmp = deepcopy(rbt)
    rbttmp.dq = zero_matrix(rbttmp.dof,1)
    rbttmp.ddq = zero_matrix(rbttmp.dof,1)
    rbttmp.gen_geom()
    return gen_tau_RNE( do_varify, rbttmp, varify_trig=varify_trig )


def gen_ccfterm_RNE( do_varify, rbt, usemotordyn=True, varify_trig=True ):
    from copy import deepcopy
    rbttmp = deepcopy(rbt)
    rbttmp.grav = zero_matrix(3,1)
    rbttmp.ddq = zero_matrix(rbttmp.dof,1)
    rbttmp.gen_geom()
    return gen_tau_RNE( do_varify, rbttmp, usemotordyn=usemotordyn, varify_trig=varify_trig )


def gen_massmatrix_RNE( do_varify, rbt, usemotordyn = True, varify_trig = True ):
  
    from copy import deepcopy
    
    Si = _gen_rbt_Si(rbt)
    
    IIi = range(0,rbt.dof+1)
    ### Linear dynamic parameters:
    for i in range(1,rbt.dof+1): IIi[i] = block_matrix( [ rbt.Ifi[i] , skew(rbt.mli[i]) , -skew(rbt.mli[i]) , rbt.mi[i]*identity_matrix(3) ] ,subdivide=False)
    ### Not linear dynamic parameters:
    #for i in range(1,dof+1): IIi[i] = block_matrix( [ rbt.Ifi_from_Ici[i] , skew(rbt.mli_e[i]) , -skew(rbt.mli_e[i]) , rbt.mi[i]*identity_matrix(3) ] ,subdivide=False)
    
    def custom_simplify( expression ):
        #return trig_reduce(expression.expand()).simplify_rational()
        return expression.simplify_rational()
    
    varify_func = None
    
    if do_varify:
        auxvars = []
        def varify_func(exp,varrepr):
            if varify_trig : exp = exp.subs( LoP_to_D(rbt.LoP_trig_f2v) )
            exp = custom_simplify( exp )
            return m_varify(exp, auxvars, poolrepr='auxM', condition_func=is_compound)
    
    M = matrix(SR,rbt.dof,rbt.dof)
    
    rbttmp = deepcopy(rbt)
    rbttmp.grav = zero_matrix(3,1)
    rbttmp.dq = zero_matrix(rbttmp.dof,1)
    
    for i in range(M.nrows()):
        rbttmp.ddq = zero_matrix(rbttmp.dof,1)
        rbttmp.ddq[i,0] = 1.0
        rbttmp.gen_geom()
        
        Vi, dVi = _forward_RNE( rbttmp, Si, varify_func )
        Mcoli = _backward_RNE( rbttmp, IIi, Si, Vi, dVi, usemotordyn, None, varify_func = varify_func )
        
        # Do like this since M is symmetric:
        M[:,i] = matrix(SR,i,1,M[i,:i].list()).stack( Mcoli[i:,:] )
    
    if do_varify:
        if varify_trig : auxvars = rbt.LoP_trig_v2f + auxvars
        return auxvars , M
    else:
        return M





def _gen_tau_atlRNE( do_varify, rbt, usemotordyn = True, varify_trig = True ):
  
    from copy import copy
    
    dof = rbt.dof
    
    Rdhi = copy(rbt.Rdhi)
    Rdhi.append( identity_matrix(SR,3) )
    
    wfi = range(0,dof+1)
    dwfi = range(0,dof+1)
    ddpfi = range(0,dof+1)
    ddpci = range(0,dof+1)
    ffi = range(0,dof+2)
    nfi = range(0,dof+2)
    tau = matrix(SR,dof,1)
    
    wfi[0] = copy(zero_matrix(SR,3,1))
    dwfi[0] = copy(zero_matrix(SR,3,1))
    ddpfi[0] = copy(zero_matrix(SR,3,1))
    ddpci[0] = copy(zero_matrix(SR,3,1))
    ffi[dof+1] = copy(zero_matrix(SR,3,1))
    nfi[dof+1] = copy(zero_matrix(SR,3,1))
    
    if do_varify:
        
        def custom_simplify( expression ):
        #    return trig_reduce(expression.expand()).simplify_rational()
            return expression.simplify_rational()
        
        def varify_func(exp,varpool,varrepr):
            if varify_trig : exp = exp.subs( LoP_to_D(rbt.LoP_trig_f2v) )
            exp = custom_simplify( exp )
            return m_varify(exp,varpool,condition_func=is_compound)
            #return v3R_varify(exp,varpool,varrepr)
        
        auxvars = []
    
    else:
        def varify_func(exp,varpool,varrepr):
            return exp
        auxvars = ()
    
    # Forward
    for i in range(1,dof+1):
        
        wfi[i] = Rdhi[i].transpose() * ( wfi[i-1] + rbt.dq[i-1,0] * rbt.zi[0] )
        wfi[i] = varify_func( wfi[i] , auxvars , 'wf_'+str(i) )
        
        dwfi[i] = Rdhi[i].transpose() * ( dwfi[i-1] + rbt.ddq[i-1,0] * rbt.zi[0] \
        + cross_product( rbt.dq[i-1,0] * wfi[i-1] , rbt.zi[0] ) )
        dwfi[i] = varify_func( dwfi[i] , auxvars , 'dwf_'+str(i) )
        
        ddpfi[i] = Rdhi[i].transpose() * ddpfi[i-1] + cross_product( dwfi[i] , rbt.pdhfi[i] ) \
        + cross_product( wfi[i] , cross_product( wfi[i] , rbt.pdhfi[i] ) )
        ddpfi[i] = varify_func( ddpfi[i] , auxvars , 'ddpf_'+str(i) )
        
        ddpci[i] = ddpfi[i] + cross_product( dwfi[i] , rbt.li[i] ) \
        + cross_product( wfi[i] , cross_product( wfi[i] , rbt.li[i] ) )
        ddpci[i] = varify_func( ddpci[i] , auxvars , 'ddpc_'+str(i) )

    
    # Backward

    if usemotordyn:
      nfi_motor = [zero_matrix(3,1) for i in xrange(dof+1)]
      for m in rbt.motors:
        Im = m[0]; qm = m[1]; l = m[2]; zm = m[3]
        dqm = qm.subs(rbt.D_q_v2f).derivative(rbt.t).subs(rbt.D_q_f2v)
        ddqm = dqm.subs(rbt.D_q_v2f).derivative(rbt.t).subs(rbt.D_q_f2v)
        nfi_motor[l] += ddqm * Im * zm + dqm * Im * cross_product( wfi[l] , zm )
    
    for i in range(dof,0,-1):
        
        ffi[i] = Rdhi[i+1] * ffi[i+1] + rbt.mi[i] * ( ddpci[i] - rbt.Ri[i].transpose()*rbt.grav )
        ffi[i] = varify_func( ffi[i] , auxvars , 'ff_'+str(i) )
        
        nfi[i] = cross_product( - ffi[i] , rbt.pdhfi[i] + rbt.li[i] ) + Rdhi[i+1] * nfi[i+1] + \
        cross_product( Rdhi[i+1] * ffi[i+1] , rbt.li[i] ) + rbt.Ici[i] * dwfi[i] + \
        cross_product( wfi[i] , rbt.Ici[i] * wfi[i] )
        if usemotordyn:
          nfi[i] += nfi_motor[i]
        nfi[i] = varify_func( nfi[i] , auxvars , 'nf_'+str(i) )
    
    for i in range(1,dof+1):
      
        tau[i-1,0] = ( nfi[i].transpose() * Rdhi[i].transpose() * rbt.zi[0] )[0,0]
        
        if usemotordyn:
          for m in rbt.motors:
            Im = m[0]; qm = m[1]; l = m[2]; zm = m[3]
            dqm = qm.subs(rbt.D_q_v2f).derivative(rbt.t).subs(rbt.D_q_f2v)
            ddqm = dqm.subs(rbt.D_q_v2f).derivative(rbt.t).subs(rbt.D_q_f2v)
            km = qm.coeff(rbt.q[i-1,0])
            if km != 0.0:
              tau[i-1,0] += km * Im * ( ( dwfi[l] + ddqm * zm + dqm * cross_product( wfi[l] , zm ) ).transpose() * zm )[0,0]
    
    if do_varify:
        if varify_trig : auxvars = rbt.LoP_trig_v2f + auxvars
        return auxvars,tau
    else:
        return tau



def _gen_gravterm_altRNE(do_varify, rbt, varify_trig=True):
    from copy import deepcopy
    rbttmp = deepcopy(rbt)
    rbttmp.dq = zero_matrix(rbttmp.dof,1)
    rbttmp.ddq = zero_matrix(rbttmp.dof,1)
    rbttmp.gen_geom()
    return _gen_tau_altRNE(do_varify, rbttmp, varify_trig)



