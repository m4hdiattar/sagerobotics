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


import utils
from sage.all import SR, exp, cos, sin, matrix, zero_matrix


def gen_potential_energy( robot ):
    
    Ui = range(0 ,robot.dof+1 )
    for i in range(1 ,robot.dof+1 ):
        Ui[i] = ( - robot.grav.transpose() * robot.pi[i] * robot.mi[i] - robot.grav.transpose() * robot.Ri[i] * robot.mli[i] )[0 ,0 ]
    U = sum(Ui)
    
    return U


def gen_kinetic_energy( robot, usemotordyn = True ):
    
    Ki = range(0 ,robot.dof+1 )
    for i in range(1 ,robot.dof+1 ):
        Ki[i] = ( ( (1.0/2.0) * robot.mi[i] * robot.dq.transpose() * robot.Jpi[i].transpose() * robot.Jpi[i] * robot.dq ) + ( robot.dq.transpose() * robot.Jpi[i].transpose() * robot.Ri[i] * utils.skew( robot.Ri[i].transpose() * robot.Joi[i] * robot.dq ) * robot.mli[i] ) + (1.0/2.0) * robot.dq.transpose() * robot.Joi[i].transpose() * robot.Ri[i] * robot.Ifi[i] * robot.Ri[i].transpose() * robot.Joi[i] * robot.dq )[0,0]

    if usemotordyn:
      Kmi = range(len(robot.motors))
      for i,m in enumerate(robot.motors):
        Im = m[0]; qm = m[1]; l = m[2]; zm = m[3]
        dqm = qm.subs(robot.D_q_v2f).derivative(robot.t).subs(robot.D_q_f2v)
        ddqm = dqm.subs(robot.D_q_v2f).derivative(robot.t).subs(robot.D_q_f2v)
        wl = robot.Ri[l].transpose() * robot.Joi[l] * robot.dq
        Kmi[i] =  dqm * Im * ( zm.transpose() * wl )[0,0] + (1.0/2.0) * dqm * dqm * Im

    K = sum(Ki)
    if usemotordyn: K += sum(Kmi)
    
    return K


def gen_lagrangian( robot, usemotordyn = True ):
    return gen_kinetic_energy(robot, usemotordyn) - gen_potential_energy(robot)


def gen_tau_EL( robot, lagrangian = None, usemotordyn = True ):
    L = lagrangian if lagrangian else gen_lagrangian( robot, usemotordyn )
    tau = matrix(SR,robot.dof,1 )
    for i in range(0 ,robot.dof):
      tau[i,0 ] = ( L.derivative(robot.dq[i,0 ]).subs(robot.D_q_v2f).derivative(robot.t).subs(robot.D_q_f2v) - L.derivative(robot.q[i,0 ]) )
    return tau


def gen_gravterm_EL( robot, potential_energy = None ):
      U = potential_energy if potential_energy else gen_potential_energy( robot )
      from copy import copy
      g = copy(zero_matrix(SR,robot.dof,1 ))
      for i in range(0 ,robot.dof):
        g[i,0 ] = U.derivative(robot.q[i,0 ])
      return g


def _gen_massmatrix_EL( robot, usemotordyn = True  ):
      print "Warning: result uses elementary dynamic parameters (opposed to grouped parameters to which dynamic model is linear)."
      if usemotordyn: raise Exception('gen_massmatrix_EL has no implementation for usemotordyn=True yet')
      from copy import copy
      M = copy(zero_matrix(SR,robot.dof,robot.dof))
      for i in range(1 ,robot.dof+1 ):
        M = M + robot.mi[i]*robot.Jcpi[i].transpose()*robot.Jcpi[i] + robot.Jcoi[i].transpose()*robot.Ri[i]* robot.Ici[i] *robot.Ri[i].transpose()*robot.Jcoi[i]
      #for i in range(1 ,robot.dof+1 ):
        #M = M + robot.mi[i]*robot.Jcpi[i].transpose()*robot.Jcpi[i] + robot.Jcoi[i].transpose()*robot.Ri[i]*( robot.Ifi[i] - robot.mi[i] * utils.skew( robot.Ri[i]*robot.li[i] ).transpose() * utils.skew( robot.Ri[i]*robot.li[i] )  )*robot.Ri[i].transpose()*robot.Jcoi[i]
      return M


def _gen_ccfmatrix_EL( robot, massmatrix = None, usemotordyn = True ):
      M = massmatrix if massmatrix else gen_massmatrix_EL( robot, usemotordyn )
      from copy import copy
      C = copy(zero_matrix(SR,robot.dof,robot.dof))
      for i,j in [(i,j) for i in range(robot.dof) for j in range(robot.dof)]:
        C[i,j] = 0 ;
        for k in range(robot.dof):
          C[i,j] += 1.0/2.0 * ( M[i,j].derivative(robot.q[k,0])  + M[i,k].derivative(robot.q[j,0]) - M[j,k].derivative(robot.q[i,0])   ) * robot.dq[k,0]
      return C


def _gen_ccfterm_EL( robot, massmatrix = None, usemotordyn = True ):
      C = gen_ccfmatrix_EL( robot, massmatrix, usemotordyn )
      v = C * robot.dq
      return v
