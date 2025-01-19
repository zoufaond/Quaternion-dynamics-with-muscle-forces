function jac = genFm_quat_Cnst(genEq)
% function create external torques from the constrained muscle length 
% based only on the q2,q3,q4 elements (q1 = sqrt(1-q2^2-q3^2-q4^2))

% two spherical joints = 2*4 quaternion elements
q = sym('q',[1 8],'real');
% 6 muscles
F_iso = sym('F_iso',[1 6],'real');
l0m = sym('l0m', [1 6],'real');
akt = sym('akt',[6 1],'real');
t = sym('t','real');

% now the muscles are not functions of q1 (q(1) and q(5) in this case)

q1Cst = [sqrt(1-q(2)^2-q(3)^2-q(4)^2),q(2),q(3),q(4)];
q2Cst = [sqrt(1-q(6)^2-q(7)^2-q(8)^2),q(6),q(7),q(8)];
qCnst = [q1Cst,q2Cst];
muscle_len(1) = muscle_length('Thorax','Clavicle',[-1, 1.2, 0],[0.1*sqrt(2)/2 0.1*sqrt(2)/2 -1],qCnst);
muscle_len(2) = muscle_length('Thorax','Scapula',[-1, 1.2, 0],[-0.1*sqrt(2)/2 0.1*sqrt(2)/2 -1],qCnst);
muscle_len(3) = muscle_length('Clavicle','Scapula',[-0.1*sqrt(2)/2 0.1*sqrt(2)/2 -0.6],[0.1 0 -0.3],qCnst);
muscle_len(4) = muscle_length('Thorax','Scapula',[-0.6, -0.6, 0],[-0.1,0,-0.5],qCnst);
muscle_len(5) = muscle_length('Thorax','Scapula',[-0.7, 0.5, 0],[0,-0.1,-0.3],qCnst);
muscle_len(6) = muscle_length('Thorax','Scapula',[0.8, -0.4, 0],[0,0.1,-0.4],qCnst);

% Jacobian of muscle lengths

jac = -(jacobian(muscle_len,[q(2:4),q(6:8)])'); %the lenths are not functions of q(1) and q(5)
% calculate the muscle (Thelen 2003 without velocity function a with rigid tendon - so no differential equation needed)
for i=1:6
    muscle_forces(i) = muscle_force(muscle_len(i),F_iso(i),akt(i),l0m(i));
end

if genEq == 1
% calculate [6x1] vector of 'quaternion-space' forces (based on the quaternion constraint)
% explained in readme
fe = [jac*muscle_forces'];
F1 = invJtrans(q(1:4))*fe(1:3);
F2 = invJtrans(q(5:8))*fe(4:6);
FE = ([F1;F2]);

% generate function
matlabFunction(FE,'file','FM_quat_Cnst','vars',{t,q,F_iso,l0m,akt});
end

function force = muscle_force(length, F_iso, akt, l0m)
    f_gauss = 0.25;
    force = (((length / l0m)^3) * exp(8 * length / l0m - 12.9) + (exp(-(length / l0m - 1)^2 / f_gauss)) * akt) * F_iso;
end

function length = muscle_length(origin, insertion, O_pos, I_pos, q)
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
        O =position(O_pos(1), O_pos(2), O_pos(3));
        TC_S = T_z(-1);
        RC_S = Qrm(q(5:8));
        I =TC_S*RC_S * position(I_pos(1), I_pos(2), I_pos(3));
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
    x = q(2);
    y = q(3);
    z = q(4);
    w = q(1);
    Rq =  [1-2*(y^2+z^2), 2*(x*y-z*w), 2*(x*z+y*w);
     2*(x*y+z*w), 1-2*(x^2+z^2), 2*(y*z-x*w);
     2*(x*z-y*w), 2*(y*z+x*w), 1-2*(x^2+y^2)];
    res = [Rq,sym(zeros(3,1));
            sym(zeros(1,3)),1];

    % a = q(1);
    % b = q(2);
    % c = q(3);
    % d = q(4);
    % Rq = [a^2 + b^2 - c^2 - d^2, 2*b*c + 2*a*d, 2*b*d - 2*a*c;
    %        2*b*c - 2*a*d, a^2 + c^2 - d^2 - b^2, 2*c*d + 2*a*b;
    %        2*b*d + 2*a*c, 2*c*d - 2*a*b, a^2 + d^2 - b^2 - c^2];
    % res = [Rq,sym(zeros(3,1));
    %         sym(zeros(1,3)),1];
end


function r = position(x,y,z)
    r = [x;y;z;1];
end

function res = T(quat)
    a = quat(1);
    b = quat(2);
    c = quat(3);
    d = quat(4);
    res = [-b/a, -c/a, -d/a;eye(3)];
end

function res = Jtrans(quat)
    res = 2*G(quat)*T(quat);
    res = res';
end

function res = G(Q)
Q0 = Q(1);
Q1 = Q(2);
Q2 = Q(3);
Q3 = Q(4);
res = [-Q1, Q0, Q3, -Q2;
        -Q2,-Q3, Q0, Q1;
        -Q3, Q2, -Q1, Q0];
end

function res = invJtrans(quat)
    q1 = quat(1);
    q2 = quat(2);
    q3 = quat(3);
    q4 = quat(4);
    res = [ q1/2,  q4/2, -q3/2;
            -q4/2,  q1/2,  q2/2;
            q3/2, -q2/2,  q1/2];
end

end