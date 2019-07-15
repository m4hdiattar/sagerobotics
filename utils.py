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


class ProgressTask:
  def __init__(self,tasksize):
    import sage
    self.tasksize = tasksize
    self.iter = 0
    self.t0 = sage.all.walltime()
  def increment(self,quantity=1):
    import sage
    self.tnow = sage.all.walltime(self.t0)
    self.iter += quantity
    return self
  def status(self):
    done = 100.0 *float( self.iter )/float( self.tasksize )
    eta = (self.tnow/done)*(100.0-done)
    return "%.1f%% - %.2fs / %.2fs - eta: %.2fs" % (done,self.tnow,self.tnow+eta,eta)


def LoP_to_D(LoP):
    D = {}
    for (a,b) in LoP: D[a] = b
    return D



def symround(expr, ndigits=0):
   if hasattr(expr,'operator') and expr.operator():
      opds = map( lambda y : symround(y,ndigits=ndigits) , expr.operands() )
      if len(opds) == 1:
        return expr.operator()(opds[0])
      else:
        r = opds[0]
        for i in xrange(1,len(opds)):
          r = expr.operator()(r,opds[i])
        return r
   try:
      from sage.all import round
      r = round(expr,ndigits=ndigits)
      if r == expr: return expr
      else: return r
   except TypeError:
      return expr


def trig_reduce(m):
  return m.apply_map(lambda i: i.trig_reduce())

def subsm(m,n):
    d = {}
    for i in range(0 ,m.nrows()):
        for j in range(0 ,m.ncols()):
            d[m[i,j]] = n[i,j]
    return d


def inverse_T(T):
  from sage.all import block_matrix, zero_matrix, Integer
  return block_matrix( [ [ T[0:3,0:3].transpose() , - T[0:3,0:3].transpose() * T[0:3,3] ] , [ zero_matrix(1,3) , Integer(1) ] ], subdivide=False )

def sagetosympy(m, print_progress = True) :
    import sympy
    sm = sympy.zeros( (m.nrows(), m.ncols()) )
    if print_progress : progress = ProgressTask( m.nrows()* m.ncols() )
    for i in range(m.nrows()) :
        for j in range(m.ncols()) :
            sm[i,j] = m[i,j]._sympy_()
            if print_progress : print progress.increment().status()
    return sm


def sympytosage(sm) :
    from copy import copy
    from sage.all import zero_matrix,SR
    m = copy(zero_matrix(SR,sm.shape[0 ],sm.shape[1 ]))
    for i in range(m.nrows()) :
        for j in range(m.ncols()) :
            m[i,j] = sm[i,j]._sage_()
    return m



def linear_dependencies( m , err_dec = 20, round_dec = 10 ):

    import numpy as np
    from sage.all import matrix, round
    
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



def cross_product(v1,v2):
  from sage.all import matrix
  return matrix( v1.transpose()[0].cross_product( v2.transpose()[0]  ) ).transpose()

liebracket = cross_product

def skew(v):
  from sage.all import matrix
  return matrix( [ [0 ,-v[2 ,0 ],v[1 ,0 ]] , [v[2 ,0 ],0 ,-v[0 ,0 ]] , [-v[1 ,0 ],v[0 ,0 ],0 ] ] )

def unskew(m):
    from sage.all import matrix
    if not m.is_skew_symmetric() :
        raise Exception('unskew: input matrix is not skew symmetric')
    return matrix( [ m[2,1], m[0,2], m[1,0] ] ).transpose()

def se3skew(g):
    from sage.all import block_matrix, zero_matrix, Integer
    w = g[0:3,0]
    v = g[3:6,0]
    return block_matrix( [ [ skew(w) , v ] , [ zero_matrix(1,3) , Integer(0) ] ] ,subdivide=False)

def se3unskew(g):
    w = unskew(g[0:3,0:3])
    v = g[0:3,3]
    return w.stack(v)


def Adj(G,g) :
    from sage.all import block_matrix, zero_matrix
    R = G[0:3,0:3]
    p = G[0:3,3]
    return block_matrix( [ [ R , zero_matrix(3,3) ] , [ skew(p)*R , R ] ] ,subdivide=False) * g

def Adjdual(G,g) :
    from sage.all import block_matrix, zero_matrix
    R = G[0:3,0:3]
    p = G[0:3,3]
    return block_matrix( [ [ R , zero_matrix(3,3) ] , [ skew(p)*R , R ] ] ,subdivide=False).transpose() * g

def adj(g,h) :
    from sage.all import block_matrix, zero_matrix
    wg = g[0:3,0]
    vg = g[3:6,0]
    return block_matrix( [ [ skew(wg) , zero_matrix(3,3) ] , [ skew(vg) , skew(wg) ] ] ,subdivide=False) * h

def adjdual(g,h) :
    from sage.all import block_matrix, zero_matrix
    wg = g[0:3,0]
    vg = g[3:6,0]
    return block_matrix( [ [ skew(wg) , zero_matrix(3,3) ] , [ skew(vg) , skew(wg) ] ] ,subdivide=False).transpose() * h



def is_compound(expression):
    import sympy
    exp = sympy.sympify(expression._sympy_())
    if exp.is_Atom or (-exp).is_Atom :
        return False
    else :
        return True


def varify( expression, varpool, varrepr = '', poolrepr = '', condition_func=None ):
    from sage.all import var, RR
    if condition_func and not condition_func(expression):
        return expression
    else:
        latexname = None
        if not varrepr:
            varrepr = str(len(varpool))
            if not poolrepr:
                poolrepr = 'auxA'
                latexname = r'\mathbf{a}_{'+str(len(varpool))+'}'
        name = poolrepr+varrepr
        from sage.all import var
        newvar = var(poolrepr+varrepr,latex_name=latexname,domain=RR)
        varpool.append( (newvar, expression) )
        return newvar


def m_varify( m_exp, varpool, varrepr = '', poolrepr = '', condition_func = None ):
    from copy import copy
    m_exp_out = copy(m_exp)
    repridx = None
    for e in [(i,j) for i in range(m_exp.nrows()) for j in range(m_exp.ncols())] :
        if varrepr : repridx = varrepr+'_'+str(e[0])+'_'+str(e[1])
        m_exp_out[e] = varify( m_exp[e] , varpool, repridx, poolrepr, condition_func )
    return m_exp_out

def v3R_varify( m_exp, varpool, varrepr = '', poolrepr = '', condition_func = None ):
    from copy import copy
    m_exp_out = copy(m_exp)
    repridx = None
    for i in range(3) :
        if varrepr : repridx = varrepr+'_'+('xyz'[i])
        m_exp_out[i] = varify( m_exp[i,0] , varpool, repridx, poolrepr, condition_func )
    return m_exp_out

def v6R_varify( m_exp, varpool, varrepr = '', poolrepr = '', condition_func = None ):
    from copy import copy
    m_exp_out = copy(m_exp)
    repridx = None
    for i in range(6) :
        if varrepr : repridx = varrepr+'_'+('uvwxyz'[i])
        m_exp_out[i] = varify( m_exp[i,0] , varpool, repridx, poolrepr, condition_func )
    return m_exp_out



