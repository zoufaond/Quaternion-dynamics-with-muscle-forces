import sympy as sp
import sympy.physics.mechanics as me
import numpy as np
import matplotlib.pyplot as plt

def plot_results(solution, vals , time, num_nodes):
    solution_mat = solution[:2*4*num_nodes].reshape(8,num_nodes)
    sol_glob = np.zeros([4,num_nodes])
    for i in range(num_nodes):
        sol_glob[:,i] = mulQuat_np(solution_mat[0:4,i],solution_mat[4:8,i])[:,0]
        # print(mulQuat_np(solution_mat[0:4,i],solution_mat[4:8,i])[:,0])

    vals_mat = vals.reshape(4,num_nodes)

    fig, axs = plt.subplots(4)
    for j in range(4):
        axs[j].plot(time,vals_mat[j,:],marker = 'o')
        axs[j].plot(time,sol_glob[j,:])
        fig.set_figheight(3)

    print(sol_glob-vals_mat)
    
    # return vals_mat

def min_rotglob(sym_vec, vals,num_nodes,interval_value):
    obj_sep = sp.zeros(len(vals),1)
    sym_mat = sym_vec.reshape(8,num_nodes)
    vals_mat = vals.reshape(4,num_nodes)

    for i in range(num_nodes):
        obj_sep[0+i*4:4+i*4,0] = mulQuat_sp(sym_mat[0:4,i],sym_mat[4:8,i]) - sp.Matrix(vals_mat[:,i])

    obj = sum(interval_value * (obj_sep.applyfunc(lambda x: x**2)))
    obj_np = (sp.lambdify(sym_vec, obj, cse = True))
    obj_jacobian = sp.Matrix([obj]).jacobian(sym_vec)
    obj_jacobian_np = sp.lambdify(sym_vec, obj_jacobian,cse=True)

    return obj_np, obj_jacobian_np

def min_rotglob2(num_coords,num_nodes,interval_value):

    for i in range(num_nodes):
        obj_sep[0+i*4:4+i*4,0] = mulQuat_sp(sym_mat[0:4,i],sym_mat[4:8,i]) - sp.Matrix(vals_mat[:,i])

    obj = sum(interval_value * (obj_sep.applyfunc(lambda x: x**2)))
    obj_np = (sp.lambdify(sym_vec, obj, cse = True))
    obj_jacobian = sp.Matrix([obj]).jacobian(sym_vec)
    obj_jacobian_np = sp.lambdify(sym_vec, obj_jacobian,cse=True)

    return obj_np, obj_jacobian_np

def mulQuat_sp(qa,qb):
    res = sp.Matrix([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
                     [qa[0]*qb[1] + qa[1]*qb[0] + qa[2]*qb[3] - qa[3]*qb[2]],
                     [qa[0]*qb[2] - qa[1]*qb[3] + qa[2]*qb[0] + qa[3]*qb[1]],
                     [qa[0]*qb[3] + qa[1]*qb[2] - qa[2]*qb[1] + qa[3]*qb[0]]])
    return res

def mulQuat_np(qa,qb):
    res = np.array([[qa[0]*qb[0] - qa[1]*qb[1] - qa[2]*qb[2] - qa[3]*qb[3]],
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