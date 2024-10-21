import sympy as sp
import sympy.physics.mechanics as me


def min_rotglob(qa,qb):
    rot_mat = mulQuat(qa,qb)
    rot_mat_J = rot_mat.jacobian(qa+qb)
    rot_mat_np = sp.lambdify([list(qa)+list(qb)],rot_mat)
    rot_mat_J_np = sp.lambdify([qa+qb],rot_mat_J)

    return rot_mat, rot_mat_J, rot_mat_np, rot_mat_J_np

def mulQuat(qa,qb):
    res = sp.Matrix([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
                     [qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2]],
                     [qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1]],
                     [qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0]]])
    return res

def Qrm(q):
    w = q[0]
    x = q[1]
    y = q[2]
    z = q[3]
    res =  sp.Matrix([[1-2*(y**2+z**2), 2*(x*y-z*w), 2*(x*z+y*w),0],
                      [2*(x*y+z*w), 1-2*(x**2+z**2), 2*(y*z-x*w),0],
                      [2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x**2+y**2),0],
                      [0           ,0          ,0             ,1]])

    return res

def YZY_seq(phi1,phi2,phi3):
    res = R_y(phi1) * R_z(phi2) * R_y(phi3)

    return res

def T_x(x):
    trans_x = sp.Matrix([[1,0,0,x],
                         [0,1,0,0],
                         [0,0,1,0],
                         [0,0,0,1]])
    return trans_x

def T_y(y):
    trans_y = sp.Matrix([[1,0,0,0],
                         [0,1,0,y],
                         [0,0,1,0],
                         [0,0,0,1]])
    return trans_y

def T_y(z):
    trans_y = sp.Matrix([[1,0,0,0],
                         [0,1,0,0],
                         [0,0,1,z],
                         [0,0,0,1]])
    return trans_y

def R_x(phix):
    rot_phix = sp.Matrix([[1,0          ,0           ,0],
                         [0,sp.cos(phix),-sp.sin(phix),0],
                         [0,sp.sin(phix), sp.cos(phix),0],
                         [0,0          ,0           ,1]])
    return rot_phix

def R_y(phiy):
    rot_phiy = sp.Matrix([[sp.cos(phiy) ,0,sp.sin(phiy),0],
                         [0           ,1,0          ,0],
                         [-sp.sin(phiy),0,sp.cos(phiy),0],
                         [0           ,0,0           ,1]])
    return rot_phiy

def R_z(phiz):
    rot_phiz = sp.Matrix([[sp.cos(phiz),-sp.sin(phiz),0,0],
                        [sp.sin(phiz), sp.cos(phiz),0,0],
                        [0           ,0            ,1,0],
                        [0           ,0            ,0,1]])
    return rot_phiz

def position(x,y,z):
    r = sp.Matrix([x,y,z,1])
    return r