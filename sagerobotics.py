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



from sage.all import *   # import sage library

'''
TODO

prismatic joints

add motor rotor inertia model

add friction model

'''

import srbtx_utils as utils


class Robot:
  
  def __init__(self,dof,name,shortname=None):
    
    self.dof = int(dof)
    self.name = name
    if shortname == None :
      self.shortname = name.replace(' ','_').replace('.','_')
    else :
      self.shortname = shortname.replace(' ','_').replace('.','_')
    
    for i in range(1 ,1 +dof):
      var('q'+str(i),domain=RR)
    
    
    (alpha , a , d , theta) = var('alpha,a,d,theta',domain=RR)
    default_params_dh = [ alpha , a , d , theta  ]
    default_T_dh = matrix( \
    [[ cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta) ], \
    [ sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta) ], \
    [ 0.0, sin(alpha), cos(alpha), d ]
    ,[ 0.0, 0.0, 0.0, 1.0] ] )
    restore('alpha,a,d,theta')
    
    default_grav = matrix([[0.0],[0.0],[-9.81]])
    
    
    self.grav = default_grav
    
    self.params_dh = default_params_dh
    self.T_dh = default_T_dh
    
    self.links_dh = []
  
    self._gensymbols()


  def __deepcopy__(self, memo):
      import copy
      dpcpy = self.__class__(self.dof,self.name,self.shortname)
      memo[id(self)] = dpcpy
      for attr in dir(self):
          try:
              setattr(dpcpy, attr, copy.deepcopy(getattr(self, attr),memo))
          except:
              setattr(dpcpy, attr, copy.copy(getattr(self, attr)))
      return dpcpy
  
  
  def __repr__(self) :
    return 'Robot instance: ' + self.name  
  
  
  def _gensymbols(self) :
    
    self.t = var('t',domain=RR)

    self.q = matrix(SR,self.dof,1 )
    self.dq = matrix(SR,self.dof,1 )
    self.ddq = matrix(SR,self.dof,1 )
    for i in range(1 ,self.dof+1 ):
      self.q[i-1 ] = var('q'+str(i),domain=RR)
      self.dq[i-1 ] = var('dq'+str(i),latex_name=r'\dot{q}_'+str(i),domain=RR)
      self.ddq[i-1 ] = var('ddq'+str(i),latex_name=r'\ddot{q}_'+str(i),domain=RR)

    self.qt = matrix(SR,self.dof,1 )
    for i in range(1 ,self.dof+1 ):
      self.qt[i-1 ,0 ] = function('q'+str(i)+'t',t,latex_name=r'q_'+str(i))
      
      self.LoP_trig_f2v = []
      for i in range(1 ,self.dof+1 ):
          self.LoP_trig_f2v.append( ( cos(self.q[i-1 ,0 ]) , var('c'+str(i),domain=RR) ) )
          self.LoP_trig_f2v.append( ( sin(self.q[i-1 ,0 ]) , var('s'+str(i),domain=RR) ) )
      self.LoP_trig_v2f = zip(zip(*self.LoP_trig_f2v)[1],zip(*self.LoP_trig_f2v)[0])
    
    self.D_q_v2f={} ; self.D_q_f2v={}
    for i in range(0 ,self.dof):
      self.D_q_v2f.update( { self.q[i,0 ]:self.qt[i,0 ], self.dq[i,0 ]:self.qt[i,0 ].derivative(t), self.ddq[i,0 ]:self.qt[i,0 ].derivative(t,2 ) } )
      self.D_q_f2v.update( { self.qt[i,0 ]:self.q[i,0 ], self.qt[i,0 ].derivative(t):self.dq[i,0 ], self.qt[i,0 ].derivative(t,2 ):self.ddq[i,0 ] } )
    
    
    self.mi = range(0 ,self.dof+1 )
    self.li = range(0 ,self.dof+1 )
    self.Ici = range(0 ,self.dof+1 )
    self.Ifi = range(0 ,self.dof+1 )
    self.Ici_from_Ifi = range(0 ,self.dof+1 )
    self.Ifi_from_Ici = range(0 ,self.dof+1 )
    self.LoP_Ii_c2f = []
    self.LoP_Ii_f2c = []
    self.mli = range(0 ,self.dof+1 )
    self.mli_e = range(0 ,self.dof+1 )
    self.LoP_mli_v2e = []

    for i in range(1 ,self.dof+1 ):
      
      self.mi[i] = var('m'+str(i),domain=RR)
      
      aux1 = var('l_'+str(i)+'x',latex_name=r'l_{'+str(i)+',x}',domain=RR)
      aux2 = var('l_'+str(i)+'y',latex_name=r'l_{'+str(i)+',y}',domain=RR)
      aux3 = var('l_'+str(i)+'z',latex_name=r'l_{'+str(i)+',z}',domain=RR)
      self.li[i] = matrix([[aux1],[aux2],[aux3]])
      
      aux1 = var('ml_'+str(i)+'x',latex_name=r'\widehat{m_'+str(i)+'l_{'+str(i)+',x}}',domain=RR)
      aux2 = var('ml_'+str(i)+'y',latex_name=r'\widehat{m_'+str(i)+'l_{'+str(i)+',y}}',domain=RR)
      aux3 = var('ml_'+str(i)+'z',latex_name=r'\widehat{m_'+str(i)+'l_{'+str(i)+',z}}',domain=RR)
      self.mli[i] = matrix([[aux1],[aux2],[aux3]])
      self.mli_e[i] = self.mi[i] * self.li[i]
      self.LoP_mli_v2e.append( ( self.mli[i][0 ,0 ] , self.mli_e[i][0 ,0 ] ) )
      self.LoP_mli_v2e.append( ( self.mli[i][1 ,0 ] , self.mli_e[i][1 ,0 ] ) )
      self.LoP_mli_v2e.append( ( self.mli[i][2 ,0 ] , self.mli_e[i][2 ,0 ] ) )
    
      auxIcxx = var('I_'+str(i)+'xx',latex_name=r'I_{'+str(i)+',xx}',domain=RR)
      auxIcyy = var('I_'+str(i)+'yy',latex_name=r'I_{'+str(i)+',yy}',domain=RR)
      auxIczz = var('I_'+str(i)+'zz',latex_name=r'I_{'+str(i)+',zz}',domain=RR)
      auxIcxy = var('I_'+str(i)+'xy',latex_name=r'I_{'+str(i)+',xy}',domain=RR)
      auxIcxz = var('I_'+str(i)+'xz',latex_name=r'I_{'+str(i)+',xz}',domain=RR)
      auxIcyz = var('I_'+str(i)+'yz',latex_name=r'I_{'+str(i)+',yz}',domain=RR)
      self.Ici[i] = matrix([[auxIcxx,auxIcxy,auxIcxz],[auxIcxy,auxIcyy,auxIcyz],[auxIcxz,auxIcyz,auxIczz]])
      
      auxIfxx = var('If_'+str(i)+'xx',latex_name=r'\hat{I}_{'+str(i)+',xx}',domain=RR)
      auxIfyy = var('If_'+str(i)+'yy',latex_name=r'\hat{I}_{'+str(i)+',yy}',domain=RR)
      auxIfzz = var('If_'+str(i)+'zz',latex_name=r'\hat{I}_{'+str(i)+',zz}',domain=RR)
      auxIfxy = var('If_'+str(i)+'xy',latex_name=r'\hat{I}_{'+str(i)+',xy}',domain=RR)
      auxIfxz = var('If_'+str(i)+'xz',latex_name=r'\hat{I}_{'+str(i)+',xz}',domain=RR)
      auxIfyz = var('If_'+str(i)+'yz',latex_name=r'\hat{I}_{'+str(i)+',yz}',domain=RR)
      self.Ifi[i] = matrix([[auxIfxx,auxIfxy,auxIfxz],[auxIfxy,auxIfyy,auxIfyz],[auxIfxz,auxIfyz,auxIfzz]])
      
      self.Ifi_from_Ici[i] = self.Ici[i] + self.mi[i] * utils.skew(self.li[i]).transpose() * utils.skew(self.li[i])
      self.Ici_from_Ifi[i] = self.Ifi[i] - self.mi[i] * utils.skew(self.li[i]).transpose() * utils.skew(self.li[i])
      for e in [(a,b) for a in range(3) for b in range (3)]:
          self.LoP_Ii_f2c.append( ( self.Ifi[i][e] , self.Ifi_from_Ici[i][e] ) )
          self.LoP_Ii_c2f.append( ( self.Ici[i][e] , self.Ici_from_Ifi[i][e] ) )
    
    Pi = [ matrix([ self.mi[i], self.mli[i][0,0], self.mli[i][1,0], self.mli[i][2,0], self.Ifi[i][0,0], self.Ifi[i][1,1], self.Ifi[i][2,2], self.Ifi[i][0,1], self.Ifi[i][0,2], self.Ifi[i][1,2] ] ).transpose() for i in range(1,self.dof+1) ]

    self.P = matrix(SR,0 ,1 )
    for i in range(0,self.dof):  self.P = self.P.stack(Pi[i])
    
    return self
  
  
  
  
  
  def gen_geom(self):
    
    self.Tdhi = range(self.dof+1)
    self.Tdhi[0] = identity_matrix(SR,4 )
    self.Tdhi_inv = range(self.dof+1)
    self.Tdhi_inv[0] = identity_matrix(SR,4 )
    self.Rdhi = range(self.dof+1 )
    self.Rdhi[0] = identity_matrix(SR,3 )
    self.pdhi = range(self.dof+1 )
    self.pdhi[0] = zero_matrix(SR,3,1 )
    self.pdhfi = range(self.dof+1 )
    self.pdhfi[0] = zero_matrix(SR,3,1 )
    
    p_dhf =  ( self.T_dh[0 :3 ,0 :3 ].transpose() * self.T_dh[0 :3 ,3 ] ).simplify_trig()
    
    T_dh_inv = utils.inverse_T(self.T_dh).simplify_trig()
    
    for l in range(self.dof):
      
      D = utils.LoP_to_D(zip(self.params_dh , self.links_dh[l]))
      
      self.Tdhi[l+1] = self.T_dh.subs(D)
      self.Tdhi_inv[l+1] = T_dh_inv.subs(D)
      self.Rdhi[l+1] = self.Tdhi[l+1][0 :3 ,0 :3 ]
      self.pdhi[l+1] = self.Tdhi[l+1][0 :3 ,3 ]
      self.pdhfi[l+1] = p_dhf.subs(D)
    
    
    self.Ti = range(self.dof+1 )
    self.Ti[0] = identity_matrix(SR,4 )
    
    for l in range(1 , self.dof+1 ):
      
      self.Ti[l] = self.Ti[l-1 ] *  self.Tdhi[l]
          
      #Ti[l] = Ti[l].simplify_rational()
      #Ti[l] = trig_reduce(Ti[l])
      #Ti[l] = Ti[l].simplify()
    
    
    self.Ri = range(0 ,self.dof+1 )
    self.pi = range(0 ,self.dof+1 )
    self.zi = range(0 ,self.dof+1 )
    
    for l in range(0 ,self.dof+1 ):
      self.Ri[l] = self.Ti[l][0 :3 ,0 :3 ]
      self.pi[l] = self.Ti[l][0 :3 ,3 ]
      self.zi[l] = self.Ri[l][0 :3 ,2 ]


    return self



  def gen_kinem( self ):

    self.Jpi = range(0 ,self.dof+1 )
    for l in range(1 ,self.dof+1 ):
      self.Jpi[l] = matrix(SR,3 ,self.dof)
      for j in range(1 ,l+1 ):
          self.Jpi[l][0 :3 ,j-1 ] = utils.skew(self.zi[j-1 ]) * ( self.pi[l] - self.pi[j-1 ] )

      #Jpi[i] = Jpi[i].simplify_rational()
      #Jpi[i] = trig_reduce(Jpi[i])
      #Jpi[i] = Jpi[i].simplify()

    self.Joi = range(0 ,self.dof+1 )
    for l in range(1 ,self.dof+1 ):
      self.Joi[l] = matrix(SR,3 ,self.dof)
      for j in range(1 ,l+1 ):
          self.Joi[l][0 :3 ,j-1 ] = (self.zi[j-1 ])

      #Joi[i] = Joi[i].simplify_rational()
      #Joi[i] = Joi[i].simplify()


    self.Jcpi = range(0 ,self.dof+1 )
    self.Jcoi = self.Joi

    for l in range(1 ,self.dof+1 ):
      self.Jcpi[l] = self.Jpi[l] - utils.skew( self.Ri[l]*self.li[l] ) * self.Joi[l]


    return self


  def J(self,link=None):
    if link == None : link = self.dof
    return self.Jpi[link].stack(self.Joi[link])


  def gen_dynEL_L(self):

    self.Ki = range(0 ,self.dof+1 )
    for i in range(1 ,self.dof+1 ):
      self.Ki[i] = ( ( (1.0/2.0) * self.mi[i] * self.dq.transpose() * self.Jpi[i].transpose() * self.Jpi[i] * self.dq ) + ( self.dq.transpose() * self.Jpi[i].transpose() * self.Ri[i] * utils.skew( self.Ri[i].transpose() * self.Joi[i] * self.dq ) * self.mli[i] ) + (1.0/2.0) * self.dq.transpose() * self.Joi[i].transpose() * self.Ri[i] * self.Ifi[i] * self.Ri[i].transpose() * self.Joi[i] * self.dq )[0,0]

      #Ki[i] = Ki[i].trig_reduce()
      #Ki[i] = Ki[i].simplify_rational()
      #Ki[i] = Ki[i].trig_reduce()
      #Ki[i] = Ki[i].simplify()
 
    self.Ui = range(0 ,self.dof+1 )
    for i in range(1 ,self.dof+1 ):
      self.Ui[i] = ( - self.grav.transpose() * self.pi[i] * self.mi[i] - self.grav.transpose() * self.Ri[i] * self.mli[i] )[0 ,0 ]

      #Ui[i] = Ui[i].trig_reduce()
      #Ui[i] = Ui[i].simplify_rational()
      #Ui[i] = Ui[i].trig_reduce()
      #Ui[i] = Ui[i].simplify()

    self.K = sum(self.Ki)
    self.U = sum(self.Ui)
    self.L = self.K - self.U

    return self


  def gen_dynEL_L_tau(self):
    self.tau = matrix(SR,self.dof,1 )
    for i in range(0 ,self.dof):
      self.tau[i,0 ] = ( self.L.derivative(self.dq[i,0 ]).subs(self.D_q_v2f).derivative(self.t).subs(self.D_q_f2v) - self.L.derivative(self.q[i,0 ]) )
    return self

  def gen_dynEL_tau_g(self):
    self.g1 = self.tau( dict( utils.subsm(self.dq,matrix(self.dof,1 )).items() + utils.subsm(self.ddq,matrix(self.dof,1 )).items()   )  )
    return self

  def gen_dynEL_tau_v(self):
    self.v1 = self.tau( utils.subsm(self.ddq,matrix(self.dof,1 )) ) - self.g1
    return self

  def gen_dynEL_tau_M(self):
    self.M1 = copy(zero_matrix(SR,self.dof,self.dof))
    for i in range(0 ,self.dof) :
      self.M1[:,i] = self.tau( dict( utils.subsm( self.ddq , identity_matrix(self.dof)[:,i] ).items() + utils.subsm( self.dq , zero_matrix(self.dof,1 ) ).items()  ) ) - self.g1
    return self



  def gen_dynEL_g(self):
      self.g2 = copy(zero_matrix(SR,self.dof,1 ))
      for i in range(0 ,self.dof):
        self.g2[i,0 ] = self.U.derivative(self.q[i,0 ])
      return self

  def gen_dynEL_M(self):
      self.M2 = copy(zero_matrix(SR,self.dof,self.dof))
      for i in range(1 ,self.dof+1 ):
        self.M2 = self.M2 + self.mi[i]*self.Jcpi[i].transpose()*self.Jcpi[i] + self.Jcoi[i].transpose()*self.Ri[i]*( self.Ifi[i] - self.mi[i] * utils.skew( self.Ri[i]*self.li[i] ).transpose() * utils.skew( self.Ri[i]*self.li[i] )  )*self.Ri[i].transpose()*self.Jcoi[i]
      return self

  def gen_dynEL_h(self):
      self.h2 = copy(zero_matrix(SR,self.dof,self.dof))
      for i,j in [(i,j) for i in range(self.dof) for j in range(self.dof)]:
        self.h2[i,j] = 0 ;
        for k in range(self.dof):
          self.h2[i,j] += 1.0/2.0 * ( self.M2[i,j].derivative(self.q[k,0])  + self.M2[i,k].derivative(self.q[j,0]) - self.M2[j,k].derivative(self.q[i,0])   ) * self.dq[k,0]
      return self

  def gen_dynEL_v(self):
      self.v2 = self.h2 * self.dq
      return self



def gen_regressor(tau,P,verify=False):
    '''
    given symbolic vectors tau and P,
    returns matrix Y
    which verifies
    Y * P = tau
    '''

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



def lineq_base(Y,P):
    
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








