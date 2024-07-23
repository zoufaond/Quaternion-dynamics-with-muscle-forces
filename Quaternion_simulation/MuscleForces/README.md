Here the muscles forces and external torques are derived.

`GenFM_seq.m` derives muscle length Jacobian when the joints are represented by YZY rotation (just to use something different than XYZ rotation).

`GenFM_quat.m` derives muscle length Jacobian in terms of all quaternion elements, when genEq=1, the function will generate FM_quat.m.

`GenFM_quat_Cnst.m` derives muscle length Jacobian in terms of q2,q3 and q4 elements, when genEq=1, the function will generate FM_quat.m.

### Equations to mapping_analysis.m:

In this [paper](https://arxiv.org/abs/0811.2889) (most of the equations in this text can be found here) the mapping between 'quaternion' forces and external torques for a rigid body (when the rotation is described by a quaternion) is derived:

$$\vec{\tau_G} = 2\textbf{G}^T \vec{\tau_E}$$

where $\vec{\tau_G}$ is a $4 \times 1$ vector of 'quaternion' torques and $\vec{\tau_E}$ is a $3 \times 1$ vector of external torque and $\textbf{G}$ is the matrix that maps between those two

$$\textbf{G} = \begin{bmatrix}-Q_1 & Q_0 & Q_3 & -Q_2\\\ -Q_2&-Q_3& Q_0& Q_1 \\\ -Q_3& Q_2& -Q_1& Q_0\end{bmatrix}.$$

Equations of motion for 1 joint represented by quaternion is:

$$\textbf{M}(Q) \vec{\dot{q}}=\vec{f(q)}$$

where $\vec{q} = (Q_0,Q_1,Q_2,Q_3,\omega_x,\omega_y,\omega_z)^T$, $\textbf{M}(Q)$ is joint-space mass matrix and $\vec{f(q)}$ are the external forces ($\omega_x$, $\omega_y$, $\omega_z$ are angular velocities).
To add a contribution of muscle forces to the equations of motion, muscle length needs to be calculated ( $l_m = f(Q_0,Q_1,Q_2,Q_3)$ ) , and then the jacobian $\textbf{J}_Q$ w.r.t. the quaternion coordinates can be calculated as:

$$\textbf{J}_Q = \frac{\partial l_m(\vec{Q})}{\partial Q_i}$$

from this external torques can be obtained:

$$\vec{\tau_E} = \frac{1}{2} \textbf{G} \textbf{J}_Q F_m$$

where $F_m$ is a scalar muscle force. Then the finished equations of motion can be created:

$$
\textbf{M}(Q) \vec{\dot{q}}=\vec{f(q)}+
\begin{pmatrix}
\vec{0}^{4 \times 1} \\
\vec{\tau_E}^{3 \times 1}
\end{pmatrix}
$$

As shown in the `mapping_analysis.m`, this method is valid only mapping 'quaternion' forces into external torques, mapping from external torque to 'quaternion' forces is not unique - 4 of 'quaternion' forces (4 DOF) to represent 3 external torques (3 DOF) - new constraint must be introduced:

## Constrained muscle lengths
We have this equations for mapping time derivative of quaternion to angular velocities (paper cited, eq.17)

$$\vec{\omega} = 2\textbf{G} \vec{\dot{Q}}$$

We have a constraint for quaternion to be unit (for rotation). We can take a derivative of this constraint w.r.t. time, so the two constraints are like this:

$$Q_0^2 + Q_1^2 + Q_2^2 + Q_3^2 = 1, \quad \dot{Q_0}Q_0 + \dot{Q_1}Q_1 + \dot{Q_2}Q_2 + \dot{Q_3}Q_3 = 0 $$

from the seconds contraint we obtain"

$$
\underbrace{\begin{bmatrix} \dot{Q_0}\\\ \dot{Q_1} \\\ \dot{Q_2} \\\ \dot{Q_3} \end{bmatrix}}\_{\dot{Q}} = 
\underbrace{\begin{bmatrix} \frac{-Q_1}{Q_0} & \frac{-Q_1}{Q_0} & \frac{-Q_1}{Q_0} \\\ 1&0&0 \\\ 0&1&0 \\\ 0&0&1 \end{bmatrix}}_\textbf{T}
\underbrace{\begin{bmatrix}  \dot{Q_1} \\\ \dot{Q_2} \\\ \dot{Q_3} \end{bmatrix}}\_{\dot{Q}^*}
$$

So now we can write:

$$\vec{\omega} = \underbrace{2\textbf{G} \textbf{T}}\_\textbf{E} \vec{\dot{Q}}^*$$

Now we can map between external torques and 'quaternion forces' as follows

$$
\vec{\tau_{G}}^* = \textbf{E}^T \vec{\tau_{E}}, \quad \vec{\tau_{E}} = (\textbf{E}^T)^{-1} \vec{\tau_{G}}^*
$$

where $\vec{\tau_G}^*$ is now a $3\times1$ vector of 'quaternion' torques acting only on $Q1,Q2$ and $Q3$ coordinates. To calculate these torques from muscle forces, length of the muscle can be calculated the same way as in the unconstrained case (s.t. $l_m = f(Q_0,Q_1,Q_2,Q_3)$ ), then do a substitution $Q_0\sqrt{1-Q_1^2-Q_2^2-Q_3^2}$ (so now $l_m = f(Q_1,Q_2,Q_3)$ ) and then calculate Jacobian w.r.t. to $Q_1,Q_2,Q_3$

```math
\textbf{J}^*_Q = \frac{\partial l_m(\vec{Q}^*)}{\partial Q_{1..3}}
```

And map it to external torques

```math
\vec{\tau_{E}} = (\textbf{E}^T)^{-1} \textbf{J}^*_Q F_M
```
