function jac = genFm_quat(genEq)

% two spherical joints = 2*4 quaternion elements
q = sym('q',[1 8],'real');
% 6 muscles
F_iso = sym('F_iso',[1 6],'real');
l0m = sym('l0m', [1 6],'real');
akt = sym('akt',[6 1],'real');
t = sym('t','real');

% muscles with random attachments (but needs to be the same as in the SimScape model)
muscle_len(1) = muscle_length('Thorax','Clavicle',[-1, 1.2, 0],[0.1, 0.1, -1],q);
muscle_len(2) = muscle_length('Thorax','Scapula',[-1, 1.2, 0],[-0.1, 0.1, -1],q);
muscle_len(3) = muscle_length('Clavicle','Scapula',[-0.4 0.5 0],[0.1 -0.6 0.3],q);
muscle_len(4) = muscle_length('Thorax','Scapula',[-0.6, -0.6, 0],[0.2,-0.2,0.5],q);
muscle_len(5) = muscle_length('Thorax','Scapula',[-0.7, 0.5, 0],[-0.3,-0.3,0],q);
muscle_len(6) = muscle_length('Thorax','Scapula',[0.8, -0.4, 0],[-0.4,0.4,-0.4],q);

% 
jac = -(jacobian(muscle_len,q)');

if genEq == 1

)
for i=1:6
    muscle_forces(i) = muscle_force(muscle_len(i),F_iso(i),akt(i),l0m(i));
end

% calculate [8x1] vector of 'quaternion-space' forces 
fe = [jac*muscle_forces'];
% map fe to external torques via G matrix (Quaternion and Dynamics, Basile Graf, page 14, eq. 28 ... you just need to rewrite the equation to calculate external torque)
% each joint has to be treated seperately
F1 = 1/2*G(q(1:4))*fe(1:4);
F2 = 1/2*G(q(5:8))*fe(5:8);
FE = ([F1;F2]);
% generate optimized function
matlabFunction(FE,'file','FM_quat','vars',{t,q,F_iso,l0m,akt});
end

function force = muscle_force(length, F_iso, akt, l0m)
    f_gauss = 0.25;
    force = (((length / l0m)^3) * exp(8 * length / l0m - 12.9) + (exp(-(length / l0m - 1)^2 / f_gauss)) * akt) * F_iso;
end

function length = muscle_length(origin, insertion, O_pos, I_pos, q)
    % length of the muscle calculated as distance between two points

    if strcmp(origin, 'Thorax') && strcmp(insertion, 'Clavicle')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        I = Qrm(q(1:4)) * position(I_pos(1), I_pos(2), I_pos(3));
        
    elseif strcmp(origin, 'Thorax') && strcmp(insertion, 'Scapula')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        RW_C = Qrm(q(1:4));
        TC_S = T_z(-1);
        RC_S = Qrm(q(5:8));
        I = RW_C * TC_S * RC_S * position(I_pos(1), I_pos(2), I_pos(3));

    elseif strcmp(origin, 'Clavicle') && strcmp(insertion, 'Scapula')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        TC_S = T_z(-1);
        RC_S = Qrm(q(5:8));
        I = TC_S * RC_S * position(I_pos(1), I_pos(2), I_pos(3));
    end

    length = sqrt((O(1) - I(1))^2 + (O(2) - I(2))^2 + (O(3) - I(3))^2);
end

function trans_x = T_x(x)
    trans_x = [1,0,0,x;
               0,1,0,0;
               0,0,1,0;
               0,0,0,1];
end

function trans_y = T_y(y)
    trans_y = [1,0,0,0;
               0,1,0,y;
               0,0,1,0;
               0,0,0,1];
end

function trans_z = T_z(z)
    trans_z = [1,0,0,0;
               0,1,0,0;
               0,0,1,z;
               0,0,0,1];
end

function res = Qrm(q)
    % rotation matrix from quaternion
    w = q(1);
    x = q(2);
    y = q(3);
    z = q(4);
    Rq =  [1-2*(y^2+z^2), 2*(x*y-z*w), 2*(x*z+y*w);
     2*(x*y+z*w), 1-2*(x^2+z^2), 2*(y*z-x*w);
     2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x^2+y^2)];
    res = [Rq,sym(zeros(3,1));
            sym(zeros(1,3)),1];
end


function r = position(x,y,z)
    r = [x;y;z;1];
end

function res = G(Q)
    % mapping function, Quaternions and Dynamics, page 9
    Q0 = Q(1);
    Q1 = Q(2);
    Q2 = Q(3);
    Q3 = Q(4);
    res = [-Q1, Q0, Q3, -Q2;
            -Q2,-Q3, Q0, Q1;
            -Q3, Q2, -Q1, Q0];
end

end