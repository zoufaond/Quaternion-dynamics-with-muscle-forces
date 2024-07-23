% some costants for the model
% this is runned every time simulink model is runned, so every simulation
% in
m = 5;
g = 9.81;
IU = [8.1667, 8.1667, 1.66667,0,0,0];
IL = [8.1667, 8.1667, 1.66667,0,0,0];
Ixx = 8.1667;
Iyy = 8.1667;
Izz = 1.66667;
l = 10;
% Child and parent coordinates of joints (EOMs and muscle length are created for this parameters, so don't change)
rigid1C = [0,0,5];
rigid1P = [0,0,-5];
rigid2C = [0,0,5];

% random initial condition, 
q1init = [rand(1,4)];
q1init = q1init/norm(q1init);
q2init = [rand(1,4)];
q2init = q2init/norm(q2init);

% damping in the joints
c = 0;



% l0m of muscles
l0m = [20,30,30,30,30,30];

%random force and activation of muscles (in this model we don't care if muscle pushes or pulls, it's just for check if the equations are correct)
force = rand(1,6)*500-250;
%activations are constant during the simulation
activation = rand(6,1);