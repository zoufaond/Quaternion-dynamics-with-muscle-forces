load("parameters.mat");

% script called in simulink before simulation

% initial position for joints
q1init = [1,0.05,0.1,0.2];
q1init = q1init/norm(q1init);
q2init = [1,0.01,0.01,0.01];
q2init = q2init/norm(q2init);

% l0m of muscles
l0m = [2,2.3,1.2,1.5,2.5,2.1];

% muscle maximum isometric force
force = [3,20,2,12,10,9];

%activations (constant during the simulation)
activation = [0.1,0.5,0.1,0.1,0.1,0.1];