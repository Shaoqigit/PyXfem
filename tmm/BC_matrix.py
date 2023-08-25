import numpy as np

def bcm_poro_fluid(porosity):
    phi = porosity
    bcm_pf = np.zeros((4,6),dtype=complex)
    bcm_pf[0,3] = 1
    bcm_pf[1,5] = 1
    bcm_pf[2,1] = 1-phi
    bcm_pf[2,2] = phi
    bcm_pf[3,4] = 1

    bcm_fp = np.zeros((4,2),dtype=complex)
    bcm_fp[0,0] = -(1-phi)
    bcm_fp[1,0] = -phi
    bcm_fp[2,1] = 1

    return bcm_pf, bcm_fp


def bcm_fluid_poro(porosity):
    """
    used for [v^s_1, v^s_3, v^t_3, sigma^s_3, sigma^s_1, sigma^f_3]"""
    phi = porosity
    bcm_fp = np.zeros((4,2),dtype=complex)
    bcm_fp[0,1] = -1
    bcm_fp[1,0] = phi
    bcm_fp[2,0] = 1-phi

    bcm_pf = np.zeros((4,6),dtype=complex)
    bcm_pf[0,1] = 1-phi
    bcm_pf[0,2] = phi
    bcm_pf[1,5] = 1
    bcm_pf[2,3] = 1
    bcm_pf[3,4] = 1

    return bcm_fp, bcm_pf

def bcm_fluid_poro2(porosity):
    """
    used for porous up formulation"""
    phi = porosity
    bcm_fp = np.zeros((4,2),dtype=complex)
    bcm_fp[0,1] = -1
    bcm_fp[1,0] = -1

    bcm_pf = np.zeros((4,6),dtype=complex)
    bcm_pf[0,2] = 1
    bcm_pf[1,4] = 1
    bcm_pf[2,0] = 1
    bcm_pf[3,3] = 1

    return bcm_fp, bcm_pf

def bcm_poro_rigid_wall():
    bcm_pw = np.zeros((4,6),dtype=complex)
    bcm_pw[0,3] = 1
    bcm_pw[0,5] = 1
    bcm_pw[1,1] = 1
    bcm_pw[2,2] = 1

    bcm_wp = np.zeros((4,2),dtype=complex)
    bcm_wp[0,0] = -1

    return bcm_pw, bcm_wp

def bcm_rigid_wall():
    bcm_w = np.zeros((3,6),dtype=complex)
    bcm_w[0,0] = 1
    bcm_w[1,1] = 1
    bcm_w[2,2] = 1

    return bcm_w

def bcm_rigid_wall2():
    """
    used for porous up formulation"""
    bcm_w = np.zeros((3,6),dtype=complex)
    bcm_w[0,1] = 1
    bcm_w[1,2] = 1
    bcm_w[2,5] = 1

    return bcm_w