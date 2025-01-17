function jac = genFM_seq()
% this one just calculates the jacobian (here moment arms w.r.t. the coordinates)
% we can use whatever sequence we want, here I use YZY just to use on
% something unusual

phi = sym('phi',[1 6]);
dq = sym('dq', [1 6]);
F_iso = sym('F_iso',[1 6]);
l0m = sym('l0m', [1 6]);
akt = sym('akt',[1 6]);
t = sym('t');
muscle_len(1) = muscle_length('Thorax','Clavicle',[-1, 1.2, 0],[0.1, 0.1, -1],phi);
muscle_len(2) = muscle_length('Thorax','Scapula',[-1, 1.2, 0],[-0.1, 0.1, -1],phi);
muscle_len(3) = muscle_length('Clavicle','Scapula',[-0.4 0.5 0],[0.1 -0.6 0.3],phi);
muscle_len(4) = muscle_length('Thorax','Scapula',[-0.6, -0.6, 0],[0.2,-0.2,0.5],phi);
muscle_len(5) = muscle_length('Thorax','Scapula',[-0.7, 0.5, 0],[-0.3,-0.3,0],phi);
muscle_len(6) = muscle_length('Thorax','Scapula',[0.8, -0.4, 0],[-0.4,0.4,-0.4],phi);

jac = -(jacobian(muscle_len,phi)');

function force = muscle_force(length, F_iso, akt, l0m)
    f_gauss = 0.25;
    force = (((length / l0m)^3) * exp(8 * length / l0m - 12.9) + (exp(-(length / l0m - 1)^2 / f_gauss)) * akt) * F_iso;
end

function length = muscle_length(origin, insertion, O_pos, I_pos, q)
    if strcmp(origin, 'Thorax') && strcmp(insertion, 'Clavicle')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        I =  R_y(q(1)) * R_z(q(2)) * R_y(q(3)) * position(I_pos(1), I_pos(2), I_pos(3));
        
    elseif strcmp(origin, 'Thorax') && strcmp(insertion, 'Scapula')
        O = position(O_pos(1), O_pos(2), O_pos(3));
        RW_C = R_y(q(1)) * R_z(q(2)) * R_y(q(3));
        TC_S = T_z(-1);
        RC_S = R_y(q(4)) * R_z(q(5)) * R_y(q(6));
        I = RW_C * TC_S * RC_S * position(I_pos(1), I_pos(2), I_pos(3));

    elseif strcmp(origin, 'Clavicle') && strcmp(insertion, 'Scapula')
        O =position(O_pos(1), O_pos(2), O_pos(3));
        TC_S = T_z(-1);
        RC_S = R_y(q(4)) * R_z(q(5)) * R_y(q(6));
        I = TC_S*RC_S * position(I_pos(1), I_pos(2), I_pos(3));
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

function rot_phix = R_x(phix)
    rot_phix = [1,0        , 0        ,0;
                0,cos(phix),-sin(phix),0;
                0,sin(phix), cos(phix),0;
                0,0        , 0        ,1];
end

function rot_phiy = R_y(phiy)
    rot_phiy = [cos(phiy),0,sin(phiy),0;
                0        ,1,0        ,0;
               -sin(phiy),0,cos(phiy),0;
                0        ,0,0        ,1];
end

function rot_phiz = R_z(phiz)
    rot_phiz = [cos(phiz),-sin(phiz),0,0;
                sin(phiz), cos(phiz),0,0;
                0           ,0      ,1,0;
                0           ,0      ,0,1];
end

function r = position(x,y,z)
    r = [x;y;z;1];
end

end
