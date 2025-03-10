# Muscle action in dynamical system with quaternions

## [EquationsOfMotion](EquationsOfMotion) folder
- Contains the derivation of equations of motion using Sympy Mechanics.
- [double_pend_quat.ipynb](EquationsOfMotion/double_pend_quat.ipynb) is a Jupyter notebook that derives Matlab functions that store space-joint mass matrix and external forces ([mm_python.m](EquationsOfMotion/mm_python.m) and [fo_python.m](EquationsOfMotion/fo_python.m)).

## [MuscleForces](MuscleForces) folder
- Contains Matlab functions that derives the muscle length Jacobians and external torques from muscle forces symbolically.

## [double_3D_pend_Quat.slx](double_3D_pend_Quat.slx)
- Simulink file with 3D double pendulum created in Simscape. Muscles are represented by point-to-point forces.
- Includes the simulation using the equations of motion derived in Sympy.
- Calls [initFunction.m](initFunction.m) that stores the muscle parameters and initial condition for simulation.

## [plot_jacobian.m](plot_jacobian.m)
- This runs the Simulink file and plots the Jacobians.

## [generate_model_struct.m](generate_model_struct.m)
- Creates [model_struct.mat](model_struct.mat) that stores the parameters for the dynamical system and muscles insertions.
