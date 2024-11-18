import sympy.physics.mechanics as me
import sympy as sp
import numpy as np
import scipy as sc
from scipy.spatial.transform import Rotation as spat

def create_eoms_moore():
    ### BY JASON MOORE ###
    # just different way how create equations of motion,
    # here one more speed is introduced, but this approach has a singularity when q0=0

    t = sp.symbols('t')
    q0A, q1A, q2A, q3A, q0B, q1B, q2B, q3B = me.dynamicsymbols('q0A, q1A, q2A, q3A, q0B, q1B, q2B, q3B')  # quaternion
    w1A, w2A, w3A, w1B, w2B, w3B = me.dynamicsymbols('w1A, w2A, w3A, w1B, w2B, w3B')  # angular velocities
    u0A, u0B = me.dynamicsymbols('u0A u0B')
    l = 5
    m = 10
    g = 9.81
    Ixx = 10
    Iyy = 1
    Izz = 5
    F = me.dynamicsymbols('F1:7')


    N = me.ReferenceFrame('frame_ground')
    N0 = me.Point('point_ground')
    N0.set_vel(N,0)

    A = me.ReferenceFrame('A')
    B = me.ReferenceFrame('B')
    mA = me.Point('mA')
    mB = me.Point('mB')
    ABj = me.Point('ABj')


    A.orient(N, 'Quaternion', [q0A, q1A, q2A, q3A])

    N_w_A = A.ang_vel_in(N)

    kinematical1 = sp.Matrix([
        u0A - q0A.diff(t),
        w1A - N_w_A.dot(A.x),
        w2A - N_w_A.dot(A.y),
        w3A - N_w_A.dot(A.z),
    ])

    B.orient(A, 'Quaternion', [q0B, q1B, q2B, q3B])

    A_w_B = B.ang_vel_in(A)

    kinematical2 = (sp.Matrix([
        u0B - q0B.diff(t),
        w1B - A_w_B.dot(B.x),
        w2B - A_w_B.dot(B.y),
        w3B - A_w_B.dot(B.z),
    ]))

    A.set_ang_vel(N, w1A*A.x + w2A*A.y + w3A*A.z)
    B.set_ang_vel(A, w1B*B.x + w2B*B.y + w3B*B.z)

    mA.set_pos(N0, -l/2 * A.z)
    mA.v2pt_theory(N0,N,A)
    FG1 = [(mA, -m * g * N.z)]

    ABj.set_pos(N0, -l * A.z)
    ABj.v2pt_theory(N0,N,A)


    mB.set_pos(ABj, -l/2 * B.z)
    mB.v2pt_theory(ABj,N,B)


    I1 = me.inertia(A, Ixx, Iyy, Izz)
    I2 = me.inertia(B, Ixx, Iyy, Izz)

    BODY = []
    BODY.append(me.RigidBody('Abody', mA, A, m, (I1, mA)))
    BODY.append(me.RigidBody('Bbody', mB, B, m, (I2, mB)))

    kinematical = sp.Matrix([[kinematical1],[kinematical2]])

    FG2 = [(mB, -m * g * N.z)]
    Torque1 = [(A, F[0]*A.x+F[1]*A.y+F[2]*A.z)]
    Torque2 = [(B, F[3]*B.x+F[4]*B.y+F[5]*B.z)]

    holonomic = sp.Matrix([[q0A**2 + q1A**2 + q2A**2 + q3A**2 - 1],
                           [q0B**2 + q1B**2 + q2B**2 + q3B**2 - 1]])
    kane = me.KanesMethod(
        N,
        [q1A, q2A, q3A, q1B, q2B, q3B],
        [w1A, w2A, w3A, w1B, w2B, w3B],
        kd_eqs=kinematical,
        q_dependent=[q0A,q0B],
        u_dependent=[u0A,u0B],
        configuration_constraints=holonomic,
        velocity_constraints=holonomic.diff(t),
    )
    (fr, frstar) = kane.kanes_equations(BODY, (FG1+FG2+Torque1+Torque2))
    
    eoms = sp.Matrix(kinematical).col_join(fr+frstar).col_join(holonomic)
    state_symbols = [q0A, q1A, q2A, q3A, q0B, q1B, q2B, q3B, w1A, w2A, w3A, w1B, w2B, w3B, u0A, u0B]
    
    return eoms, state_symbols, F

def create_trajectory(num_nodes, duration, interval_value):
    trajectory = np.array([])
    d_trajectory = np.array([])
    time = np.linspace(0.0, duration, num=num_nodes)
    np.random.seed(5)
    # weight = np.random.randn(6)*0.4
    weight = np.ones(6)*0.3
    # weight = [1.5,1.2,0.2,1.1,2.1,1]
    for i in range(6):
        icos = (-np.cos(time*np.pi)+1)*weight[i]-0.11
        trajectory = np.append(trajectory,icos)
        d_icos = np.concatenate((np.zeros(1),np.diff(icos)/interval_value))
        d_trajectory = np.append(d_trajectory,d_icos)
        
    return trajectory, d_trajectory

def eul2quat_traj(num_nodes, trajectory, interval):
    split_traj = np.split(trajectory,2)
    resh_traj = []
    quat_res = np.zeros([2,4,num_nodes])
    w = np.zeros([2,3,num_nodes])
    u0 = []
    
    for jnt in range(2):
        resh_traj.append(split_traj[jnt].reshape(3,num_nodes).T)
        rot = spat.from_euler('YZY',resh_traj[jnt])
        quat = rot.as_quat(scalar_first=True).T
        quat_res[jnt,:,:] = quat
        dquat = np.concatenate((np.zeros((4,1)),np.diff(quat)),axis=1)/interval
        u0.append(dquat[0,:])
        for i in range(num_nodes):
            w[jnt,:,i] = 2*G(quat[:,i])@dquat[:,i]
    
    Q_my = quat_res.flatten()
    init_guess_my = np.block([quat_res.flatten(),w.flatten()])
    Q_moore = Q_my
    init_guess_moore = np.concatenate((init_guess_my,np.array(u0).flatten()))
    
    return Q_my,init_guess_my,Q_moore,init_guess_moore

def interpolate_results(num_nodes, results, num_vars, num_nodes_new, time):
    res = []
    time_new = np.linspace(0.0, time[-1], num=num_nodes_new)
    
    for i in range(num_vars):
        ivar = results[i*num_nodes:(i+1)*num_nodes]
        cs_ivar = sc.interpolate.CubicSpline(time,ivar)
        cs_vals = list(cs_ivar(time_new))
        res = res + cs_vals
        
    return np.array(res)

def G(quat):
    Q0 = quat[0];
    Q1 = quat[1];
    Q2 = quat[2];
    Q3 = quat[3];
    res = np.array([[-Q1, Q0, Q3, -Q2],
                    [-Q2,-Q3, Q0, Q1],
                    [-Q3, Q2, -Q1, Q0]])
    return res
        
