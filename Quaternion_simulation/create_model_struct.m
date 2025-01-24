clear all
% create struct with the model parameters
model.m = 1;
model.g = 9.81;
model.Ixx = 8.1667;
model.Iyy = 8.1667;
model.Izz = 1.66667;
model.l = 1;
model.radius = 0.1;
model.c = 1;

% Child and parent coordinates from center of mass
model.rigid1C = [0,0,model.l/2];
model.rigid1P = [0,0,-model.l/2];
model.rigid2C = [0,0,model.l/2];


% muscle parameters
%insertion_body = 1 ... Ground
%                 2 ... First body
%                 3 ... Second body

model.muscle1.origin = [-1, 1.2, 0];
model.muscle1.origin_body = 1;
model.muscle1.insertion = [0.1*sqrt(2)/2 0.1*sqrt(2)/2 -1];
model.muscle1.insertion_body = 2;

model.muscle2.origin = [-1, 1.2, 0];
model.muscle2.origin_body = 1;
model.muscle2.insertion = [-0.1*sqrt(2)/2 0.1*sqrt(2)/2 -1];
model.muscle2.insertion_body = 3;

model.muscle3.origin = [-0.1*sqrt(2)/2 0.1*sqrt(2)/2 -0.6];
model.muscle3.origin_body = 2;
model.muscle3.insertion = [0.1 0 -0.3];
model.muscle3.insertion_body = 3;

model.muscle4.origin = [-0.6, -0.6, 0];
model.muscle4.origin_body = 1;
model.muscle4.insertion = [-0.1,0,-0.5];
model.muscle4.insertion_body = 3;

model.muscle5.origin = [-0.7, 0.5, 0];
model.muscle5.origin_body = 1;
model.muscle5.insertion = [0,-0.1,-0.3];
model.muscle5.insertion_body = 3;

model.muscle6.origin = [0.8, -0.4, 0];
model.muscle6.origin_body = 1;
model.muscle6.insertion = [0,0.1,-0.4];
model.muscle6.insertion_body = 3;

save('parameters.mat','model')