load("model_struct.mat");

% script called in simulink before simulation

% initial position for joints
q1init = [1,0.1,-0.05,0.05];
q1init = q1init/norm(q1init);
q2init = [1,0.1,0,0.1];
q2init = q2init/norm(q2init);

rotm1init = quat2rotm(q1init);
rotm2init = quat2rotm(q2init);
eul1init = rotm2eul(rotm1init,"YZY");
eul2init = rotm2eul(rotm2init,"YZX");

% l0m of muscles
lceopt = [2,2.3,1.2,1.5,2.5,2.1];

% muscle maximum isometric forces
fmax = [10,25,12,8,7,5]*1.2;
% rng(4)
% fmax = rand(1,6)*20;

%activations (constant during the simulation)
activation = [0.1,1,0.1,1,0.1,0.1];