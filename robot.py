# -*- coding: utf-8 -*-

###############################################################################
#  Sage Robotics: Robotics Toolbox for Sage Mathematical Software
#
#      Copyright (C) 2010 CristÃ³vÃ£o Sousa <crisjss@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL),
#  version 2 or any later version.  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
###############################################################################


'''
TODO

prismatic joints

'''


import utils
import dyn_RNE
from sage.all import var, RR, SR, matrix, identity_matrix, zero_matrix, cos, sin, function


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
    
    
    (alpha , a , d , theta) = SR.var('alpha,a,d,theta',domain=RR)
    default_params_dh = [ alpha , a , d , theta  ]
    default_T_dh = matrix( \
    [[ cos(theta), -sin(theta)*cos(alpha), sin(theta)*sin(alpha), a*cos(theta) ], \
    [ sin(theta), cos(theta)*cos(alpha), -cos(theta)*sin(alpha), a*sin(theta) ], \
    [ 0.0, sin(alpha), cos(alpha), d ]
    ,[ 0.0, 0.0, 0.0, 1.0] ] )
    
    #g_a = SR.var('g_a',domain=RR)
    g_a = 9.81
    self.grav = matrix([[0.0],[0.0],[-g_a]])
    
    self.params_dh = default_params_dh
    self.T_dh = default_T_dh
    
    self.links_dh = []
    
    self.motors = []
    self.Imzi = []
    
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
  
  
  def _gensymbols(self, usefricdyn=True) :
    
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
    self.fvi = range(0 ,self.dof+1 )
    self.fci = range(0 ,self.dof+1 )
    self.Ici = range(0 ,self.dof+1 )
    self.Ifi = range(0 ,self.dof+1 )
    self.Ici_from_Ifi = range(0 ,self.dof+1 )
    self.Ifi_from_Ici = range(0 ,self.dof+1 )
    self.LoP_Ii_c2f = []
    self.LoP_Ii_f2c = []
    self.mli = range(0 ,self.dof+1 )
    self.mli_e = range(0 ,self.dof+1 )
    self.LoP_mli_v2e = []
    self.D_mli_v2e = {}

    for i in range(1 ,self.dof+1 ):

      self.mi[i] = var('m'+str(i),domain=RR)

      self.fvi[i] = var('fv'+str(i),latex_name=r'{fv}_{'+str(i)+'}',domain=RR)
      self.fci[i] = var('fc'+str(i),latex_name=r'{fc}_{'+str(i)+'}',domain=RR)
      
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
    
    self.D_mli_v2e = utils.LoP_to_D( self.LoP_mli_v2e )
    self.D_Ii_f2c = utils.LoP_to_D( self.LoP_Ii_f2c )
    self.D_Ii_c2f = utils.LoP_to_D( self.LoP_Ii_c2f )
    self.D_parm_lin2elem =  utils.LoP_to_D( self.LoP_mli_v2e + self.LoP_Ii_f2c )
    
    Pi = [ matrix([ self.mi[i], self.mli[i][0,0], self.mli[i][1,0], self.mli[i][2,0], self.Ifi[i][0,0], self.Ifi[i][0,1], self.Ifi[i][0,2], self.Ifi[i][1,1], self.Ifi[i][1,2], self.Ifi[i][2,2] ] ).transpose() for i in range(1,self.dof+1) ]

    return self


  def Parms(self, usemotordyn = True, usefricdyn=True):
    
    P = matrix(SR,0 ,1 )
    for i in range(0+1,self.dof+1):

      Pi = matrix([ self.mi[i] ])
      Pi = Pi.stack( matrix([ self.mli[i][0,0], self.mli[i][1,0], self.mli[i][2,0] ]).transpose() )
      Pi = Pi.stack( matrix([ self.Ifi[i][0,0], self.Ifi[i][0,1], self.Ifi[i][0,2], self.Ifi[i][1,1], self.Ifi[i][1,2], self.Ifi[i][2,2] ]).transpose() )
      if usefricdyn: Pi = Pi.stack( matrix([ [self.fvi[i]], [self.fci[i]] ]) )

      P = P.stack(Pi)

    if usemotordyn:
      for i in xrange(len(self.motors)):
        P = P.stack( matrix([[ self.motors[i][0] ]]) )

    return P


  def _D_grav_2_zero(self):
    variables = matrix(SR,self.grav).variables()
    Dg2z = {}
    for i in variables: Dg2z[i] = 0
    if self.grav.subs(Dg2z).is_zero():
        return Dg2z
    else:
        return None


  def add_motor(self,position,link,longaxis):
    im = len(self.motors)+1
    Imz = var( 'I_M'+str(im), latex_name=r'I_{M'+str(im)+'}', domain=RR )
    self.Imzi.append( Imz )
    self.motors.append( (Imz,position,link,longaxis) )
  
  
  
  
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
    self.Jpi[0] = zero_matrix(SR,3,self.dof)
    for l in range(1 ,self.dof+1 ):
      self.Jpi[l] = matrix(SR,3 ,self.dof)
      for j in range(1 ,l+1 ):
          self.Jpi[l][0 :3 ,j-1 ] = utils.skew(self.zi[j-1 ]) * ( self.pi[l] - self.pi[j-1 ] )

      #Jpi[i] = Jpi[i].simplify_rational()
      #Jpi[i] = trig_reduce(Jpi[i])
      #Jpi[i] = Jpi[i].simplify()

    self.Joi = range(0 ,self.dof+1 )
    self.Joi[0] = zero_matrix(SR,3,self.dof)
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
