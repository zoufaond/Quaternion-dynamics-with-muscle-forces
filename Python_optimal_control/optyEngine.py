import sympy.physics.mechanics as me
import sympy as sp
import numpy as np
import scipy as sc
from scipy.spatial.transform import Rotation as spat
import matplotlib.pyplot as plt

def create_trajectory(num_nodes, duration, interval_value):
    trajectory = np.zeros([6,num_nodes])
    d_trajectory = np.zeros([6,num_nodes])
    time = np.linspace(0.0, duration, num=num_nodes)
    np.random.seed(5)
    # weight = np.random.randn(6)*0.4
    # weight = np.ones(6)*0.5
    weight = [0.0,0.5,0.0,0.0,0.0,0.0]
    offset = [0.2,0.0,-0.2,0.0,0.0,0.0]
    for i in range(6):
        icos = (-np.cos(time*np.pi)+1)*weight[i]+offset[i]
        trajectory[i,:] = icos
        d_icos = np.concatenate((np.zeros(1),np.diff(icos)/interval_value))
        d_trajectory[i,:] = d_icos
        
    return trajectory, d_trajectory

def eul2quat_traj(num_nodes, euler, interval):
    trajectory = np.zeros([8,num_nodes])
    w = np.zeros([8,num_nodes])
    
    for jnt in range(2):
        rot = spat.from_euler('YZY',euler[jnt*3:(jnt+1)*3].T)
        quat = rot.as_quat(scalar_first=True).T
        trajectory[jnt*4:(jnt+1)*4,:] = quat
        dquat = np.concatenate((np.zeros((4,1)),np.diff(quat)),axis=1)/interval
        for i in range(num_nodes):
            w[jnt*3:(jnt+1)*3,i] = 2*G(quat[:,i])@dquat[:,i]
        w[jnt+6,:] = dquat[0,:]
    
    return trajectory,w

def interpolate_results(num_nodes, results, num_vars, num_nodes_new, time):
    res = []
    time_new = np.linspace(0.0, time[-1], num=num_nodes_new)
    
    for i in range(num_vars):
        ivar = results[i*num_nodes:(i+1)*num_nodes]
        cs_ivar = sc.interpolate.CubicSpline(time,ivar)
        cs_vals = list(cs_ivar(time_new))
        res = res + cs_vals
        
    return np.array(res)


def objective_function(num_coords,interval_value):
    x = sp.Matrix(sp.symbols(f'x1:{num_coords+1}'))
    x_traj = sp.Matrix(sp.symbols(f'x_traj1:{num_coords+1}'))

    obj_traj = sp.Matrix([interval_value*sum((sp.Matrix(x)-sp.Matrix(x_traj)).applyfunc(lambda x: x**2))])
    obj_traj_jac = obj_traj.jacobian(x)[:]

    obj_traj_np = sp.lambdify((x,x_traj),obj_traj)
    obj_traj_jac_np = sp.lambdify((x,x_traj),obj_traj_jac)

    return obj_traj_np,obj_traj_jac_np

def objective_rotglob_quat(num_coords,interval_value):
    x = sp.Matrix(sp.symbols(f'x1:{num_coords+1}'))
    x_traj = sp.Matrix(sp.symbols(f'x_traj1:{num_coords+1}'))

    rotglob1 = sp.Matrix(x[0:4])
    rotglob2 = sp.Matrix([mulQuat_sp(x[0:4],x[4:8])])
    follow_point = Qrm(x[0:4])*sp.Matrix([1,0,0,0])
    print(follow_point)

    obj_traj1 = sp.Matrix([interval_value*sum((rotglob1-sp.Matrix(x_traj[0:4])).applyfunc(lambda x: x**2))])
    obj_traj2 = sp.Matrix([interval_value*sum((rotglob2-sp.Matrix(x_traj[4:8])).applyfunc(lambda x: x**2))])
    obj = obj_traj1 + obj_traj2
    obj_traj_jac = (obj).jacobian(x)[:]

    obj_traj_np = sp.lambdify((x,x_traj),obj)
    obj_traj_jac_np = sp.lambdify((x,x_traj),obj_traj_jac)

    return obj_traj_np,obj_traj_jac_np

def objective_rotglob_eul(num_coords,interval_value):
    x = sp.Matrix(sp.symbols(f'x1:{num_coords+1}'))
    x_traj = sp.Matrix(sp.symbols(f'x_traj1:{num_coords+1}'))

    rotglob1 = sp.Matrix(x[0:3])
    rotglob2_rotmat = R_y(x[0])*R_z(x[1])*R_x(x[2])*R_y(x[3])*R_z(x[4])*R_x(x[5])
    z_rot = sp.asin(rotglob2_rotmat[1,0])
    x_rot = sp.atan2(-rotglob2_rotmat[1,2],rotglob2_rotmat[1,1])
    y_rot = sp.atan2(-rotglob2_rotmat[2,0],rotglob2_rotmat[0,0])
    rotglob2 = sp.Matrix([x_rot,y_rot,z_rot])

    obj_traj1 = sp.Matrix([interval_value*sum((rotglob1-sp.Matrix(x_traj[0:3])).applyfunc(lambda x: x**2))])
    obj_traj2 = sp.Matrix([interval_value*sum((rotglob2-sp.Matrix(x_traj[3:6])).applyfunc(lambda x: x**2))])
    obj = obj_traj1 + obj_traj2
    obj_traj_jac = obj.jacobian(x)[:]

    obj_traj_np = sp.lambdify((x,x_traj),obj)
    obj_traj_jac_np = sp.lambdify((x,x_traj),obj_traj_jac)

    return obj_traj_np,obj_traj_jac_np


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

    # print(sol_glob-vals_mat)
    
    # return vals_mat

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

def G(quat):
    Q0 = quat[0];
    Q1 = quat[1];
    Q2 = quat[2];
    Q3 = quat[3];
    res = np.array([[-Q1, Q0, Q3, -Q2],
                    [-Q2,-Q3, Q0, Q1],
                    [-Q3, Q2, -Q1, Q0]])
    return res
        
