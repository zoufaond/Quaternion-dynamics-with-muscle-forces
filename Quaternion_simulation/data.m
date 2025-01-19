% some costants for the model
% this is runned every time simulink model is runned, so every simulation
% in
m = 1;
g = 9.81;
IU = [8.1667, 8.1667, 1.66667,0,0,0];
IL = [8.1667, 8.1667, 1.66667,0,0,0];
Ixx = 8.1667;
Iyy = 8.1667;
Izz = 1.66667;
l = 1;
% Child and parent coordinates of joints (EOMs and muscle length are created for this parameters, so don't change)
rigid1C = [0,0,0.5];
rigid1P = [0,0,-0.5];
rigid2C = [0,0,0.5];

% random initial condition, 
% q1init = [rand(1,4)];
q1init = [1,0.05,0.1,0.2];
q1init = q1init/norm(q1init);
% % q2init = [rand(1,4)];
q2init = [1,0.01,0.01,0.01];
q2init = q2init/norm(q2init);

% damping in the joints
c = 0;



% l0m of muscles
l0m = [2,2.3,1.2,1.5,2.5,2.1];
rng(2)
%random force and activation of muscles
% force = rand(1,6)*50;
force = [3,20,2,12,10,9];
%activations are constant during the simulation
activation = [0.1,0.5,0.1,0.1,0.1,0.1];