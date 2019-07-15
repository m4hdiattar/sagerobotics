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



import utils


def gen_fricterm( rbt ):
  from sage.all import zero_matrix, SR, sgn
  fric = zero_matrix(SR, rbt.dof, 1 )
  for i in xrange(1,rbt.dof+1):
    fric[i-1,0] = rbt.fvi[i] * rbt.dq[i-1,0] + rbt.fci[i] * sgn(rbt.dq[i-1,0])
  return fric

def tau_2_g(robot,tau):
    g = tau( dict( utils.subsm(robot.dq,zero_matrix(robot.dof,1 )).items() + utils.subsm(robot.ddq,zero_matrix(robot.dof,1 )).items()   )  )
    return g

def tau_2_v(robot,tau,g=None):
    Dg = robot._D_grav_2_zero()
    if Dg:
        v = tau.subs( dict( utils.subsm(robot.ddq,zero_matrix(robot.dof,1 )).items() + Dg.items() ) )
    else:
        if not g:
            raise Exception('tau_2_v: needs the g parameter')
        v = tau.subs( dict( utils.subsm(robot.ddq,zero_matrix(robot.dof,1 )) ) ) - g
    return v

def tau_2_M(robot,tau,g=None):
    from sage.all import zero_matrix, SR
    from copy import copy
    M = copy(zero_matrix(SR,robot.dof,robot.dof))
    Dg = robot._D_grav_2_zero()
    if Dg:
        for i in range(0 ,robot.dof) :
            M[:,i] = tau( dict( utils.subsm( robot.ddq , identity_matrix(robot.dof)[:,i] ).items() + utils.subsm( robot.dq , zero_matrix(robot.dof,1 ) ).items() + Dg.items()  ) )
    else:
        if not g:
            raise Exception('tau_2_M: needs the g parameter')
        for i in range(0 ,robot.dof) :
            M[:,i] = tau( dict( utils.subsm( robot.ddq , identity_matrix(robot.dof)[:,i] ).items() + utils.subsm( robot.dq , zero_matrix(robot.dof,1 ) ).items()  ) ) - g
    return M


def tau_to_regressor(tau,P,verify=False):
    '''
    given symbolic vectors tau and P,
    returns matrix Y
    which verifies
    Y * P = tau
    '''
    from sage.all import zero_matrix, SR

    m = tau.nrows()
    n = P.nrows()

    Y = copy( zero_matrix(SR,m,n) )

    for i in range(0 ,m):
        for j in range(0 ,P.nrows()):
            if verify:
                coeffs = ( tau[i,0 ].coeffs(P[j,0 ]) )
                for coeff in coeffs:
                    if coeff[1] > 1 or coeff[1] < 0 :
                        raise Exception('gen_dyn_lineqs: '+P[j,0 ].__repr__()+' is not linear')
                    if coeff[1] == 1 :
                        Y[i,j] = coeff[0]
            else:
                Y[i,j] = tau[i,0].coeff( P[j,0] )

    if verify and bool( ( Y * P - tau ).expand() != zero_matrix(m,1) ):
        raise Exception('gen_dyn_lineqs: linearity not verified')

    return Y



def parm_identbase(Y,P):
    from sage.all import zero_matrix
    
    m = Y.nrows()
    n = P.nrows()
    
    nonzerol = range(0,n)
    for i in range(0,n) :
    #if Y[:,i].is_zero() :
        if bool(Y[:,i] == zero_matrix(m,1)) :
            nonzerol.remove(i)
    
    PB = P[nonzerol,:]
    YB = Y[:,nonzerol]
    
    return YB,PB


def parm_all2base( Pvals, Psyms, PBsyms ):
    return PBsyms.subs( utils.LoP_to_D( zip( Psyms.list(), Pvals.list() ) ) )

def parm_base2all( PBvals, Psyms, PBsyms ):
    return Psyms.subs( utils.LoP_to_D( zip( PBsyms.list(), PBvals.list() ) ) )


def linear_dependencies( m , err_dec, round_dec ):

    import numpy as np
    from sage.all import matrix
    
    W = m.numpy()
    (nr,nc) = W.shape
    
    d = np.identity(nc)
    
    if not W[:,0].any():
        d[0,0] = 0.0

    for i in range(1,nc):
        if not W[:,i].any():
            d[i,i] = 0.0
        else:
            A = W[:,:i]
            b = W[:,i:i+1]
            Api = np.linalg.pinv(A)
            x = np.dot(Api,b)
            err = np.linalg.norm(b - np.dot(A,x))
            if err < 10.0**(-err_dec):
                d[:,i:i+1] = np.vstack(( x.reshape(i,1),np.zeros((nc-i,1)) ))
                W[:,i:i+1] = np.zeros((nr,1))
    
    return matrix(d).apply_map(lambda i: round(i,round_dec))


def parm_deps_numeric( dof, Pn, Y_func, samples=100, err_dec = 5, round_dec = 6 ):

    import numpy
    from sage.all import matrix

    Wnp = numpy.zeros( ( dof*samples, Pn ) )
    for i in range(samples):
        q = [ float( random()*2.0*pi - pi ) for j in range(dof) ]
        dq = [ float( random()*2.0*pi - pi ) for j in range(dof) ]
        ddq = [ float( random()*2.0*pi - pi ) for j in range(dof) ]
        Wnp[ i*dof : i*dof+dof , : ] = numpy.array( Y_func( q, dq, ddq ) ).reshape( dof, Pn )

    W = matrix(Wnp)

    deps = linear_dependencies( W, err_dec, round_dec )
    
    return deps


def parm_f2c( P, dof, usefricdyn = False ):
  
  from copy import copy
  from sage.all import matrix

  Pc = copy(P)
  P = P.list()
  
  for i in range(dof):
    
    o = i * (10 + 2*usefricdyn)
    oI = o+4
    
    m = P[o+0]
    ml = matrix([[P[o+1]],[P[o+2]],[P[o+3]]])
    If = matrix([ [ P[oI+0], P[oI+1], P[oI+2] ],
                  [ P[oI+1], P[oI+3], P[oI+4] ],
                  [ P[oI+2], P[oI+4], P[oI+5] ] ])
    
    l = ml / m
    Ic = If - m * utils.skew(l).transpose() * utils.skew(l)
    
    Pc[o+0,0] = m
    Pc[o+1,0] = l[0,0]
    Pc[o+2,0] = l[1,0]
    Pc[o+3,0] = l[2,0]
    Pc[oI+0,0] = Ic[0,0]
    Pc[oI+1,0] = Ic[0,1]
    Pc[oI+2,0] = Ic[0,2]
    Pc[oI+3,0] = Ic[1,1]
    Pc[oI+4,0] = Ic[1,2]
    Pc[oI+5,0] = Ic[2,2]
  
  return Pc


def parm_c2f( P, dof, usefricdyn = False ):

  from copy import copy
  from sage.all import matrix
  
  Pf = copy(P)
  P = P.list()
  
  for i in range(dof):
    
    o = i * (10 + 2*usefricdyn)
    print i * (10 + 2*usefricdyn)
    oI = o+4
    
    m = P[o+0]
    l = matrix([[P[o+1]],[P[o+2]],[P[o+3]]])
    Ic = matrix([ [ P[oI+0], P[oI+1], P[oI+2] ],
                  [ P[oI+1], P[oI+3], P[oI+4] ],
                  [ P[oI+2], P[oI+4], P[oI+5] ] ])
    
    ml = l * m
    If = Ic + m * utils.skew(l).transpose() * utils.skew(l)
    
    Pf[o+0,0] = m
    Pf[o+1,0] = ml[0,0]
    Pf[o+2,0] = ml[1,0]
    Pf[o+3,0] = ml[2,0]
    Pf[oI+0,0] = If[0,0]
    Pf[oI+1,0] = If[0,1]
    Pf[oI+2,0] = If[0,2]
    Pf[oI+3,0] = If[1,1]
    Pf[oI+4,0] = If[1,2]
    Pf[oI+5,0] = If[2,2]
  
  return Pf
