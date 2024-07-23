genEq = 1; % don't generate external torques equations, we already have them (it also takes quite a long time)

q = sym('q',[1 8],'real');
phi = sym('phi',[1 6]);

% generate the jacobians of muscle lengths using all quaternions, length
% using the constraint and only 3 elements of quaternion and muscle length
% using YZY rotations for spherical joint
JQsym = genFm_quat(genEq); % unconstrained muscle lengths jacobian
JQCnstsym = genFm_quat_Cnst(genEq); % constrained muscle lengths jacobian
Jseqsym = genFM_seq(); % muscle lengths jacobian for YZY sequence of rotations (both joints have YZY)

% create anonymous functions
JQ = matlabFunction(JQsym,'Vars',{q});
JQCnst = matlabFunction(JQCnstsym,'Vars',{q});
Jseq = matlabFunction(Jseqsym,'Vars',{phi});
%% so we don't recalculate the symbolic computations every run
% ctrl+enter to run only this part
clc

% random configuration of our muscle jacobians in full range of YZY
% rotations
phiVal1 = rand(1,3)*2*pi-pi;
phiVal2 = rand(1,3)*2*pi-pi;
phiVal = [phiVal1,phiVal2];
quat1 = eul2quat(phiVal1,'YZY'); % get quaternion from YZY sequence
quat2 = eul2quat(phiVal2,'YZY'); % get quaternion from YZY sequence
quat = [quat1,quat2];

% calculate jacobians only for unconstrained muscle length and YZY sequence
% muscle lengths only, we'll get to the constrained one eventually
JQval = JQ(quat);
Jseqval = Jseq(phiVal);

disp('UNCONSTRAINED MUSCLE ANALYSIS:')
% UNCONSTRAINED MUSCLE LENGTH MODEL (muscle length is function of all quaternion elements)
% this maps the quaternion muscle length jacobians to YZY sequence muscle
% lengths jacobians, JseqQerror is the sum of errors of all elements, note
% that if we map from quaternion elements
SpatJfromQ = quat2spatial(JQval,quat,2);
JseqfromJQ = spatial2eul(SpatJfromQ,phiVal,2);
JseqQerror = sum(JseqfromJQ-Jseqval,'all');
% this works perfectly
disp(['Error sum: (J seq) - (J quat from J seq) = ', num2str(JseqQerror)])

% this maps YZY sequence muscle length Jacobians to quaternion Jacobians of
% nonconstrained muscle
SpatJfromseq = eul2spatial(Jseqval,phiVal,2);

% note that we now have both Jacobians in spatial coords and they are the
% same
SpatJseqSpatJQerr = sum(SpatJfromQ-SpatJfromseq,"all");
disp(['Error sum: (spatial J seq) - (spatial J quat) = ', num2str(SpatJseqSpatJQerr)])
% so we have the same result of Jacobians for spatial coords when mapping
% from sequence and from quaternions


% However, when mapping this spatial Jacobian to back quaternion jacobians
JQfromJseq = spatial2quat(SpatJfromseq,quat,2);
JQfromJseqerr = sum(JQfromJseq-JQval,"all");
disp(['Error sum: (quat J mapped from J seq) - (quat J computed) = ', num2str(JQfromJseqerr)])
% now you can see this doesn't work, but when you map this 'incorrect'
% quat Jacobian (mapped from YZY seq) - you'll get correct spatial J
back2spatial = quat2spatial(JQfromJseq,quat,2);
back2spaterror = sum(back2spatial-SpatJfromseq,'all');
disp(['Error sum: (spatial J quat from "incorrect" J quat) - (spatial J quat) = ', num2str(back2spaterror)])
% so as you can see, we were able to get to the correct spatial jacobians
% even from these incorrect quaternion Jacobians. This might be because there are no
% constraints in terms of the quaternion forces (3 spatial torques are represented by 4 quaternion 'forces')
% So it is possible (I don't really know that) that there are infinite
% possible configuration of these quaternion forces to represent the same
% external (spatial forces). So if we calculate excact muscle length, we
% can use this approach to calculate external torques but we can't take
% moment arms from opensim and map into quaternion Jacobians because from
% these examples it does not match.

disp('---------------------------')
disp('CONSTRAINED MUSCLE ANALYSIS:')
% the procedure of this analysis will be the same as in the unconstrained
% muscle lengths, just the mapping for spatial J to quat J and vice versa
% is different, this is derived in readme

%due to the constraint, this mapping doesn't work when q1 < 0, so you need
%to multiply the quaternion by -1 if q1 <0
% if you uncomment these 2 conditions, you'll se that the mapping doesn't
% work
if quat1(1) < 0
    quat1 = quat1*(-1);
end

if quat2(1) < 0
    quat2 = quat2*(-1);
end

quat = [quat1,quat2];
JQCnstval = JQCnst(quat);

% map constrained muscle length Jacobian to YZY sequence Jacobian
SpatJfromQCnst = quat2spatialCnst(JQCnstval,quat,2);
JseqfromJQCnst = spatial2eul(SpatJfromQCnst,phiVal,2);
JseqQerrorCnst = sum(JseqfromJQCnst-Jseqval,'all');
disp(['Error sum: (J seq) - (J quat from J seq) = ', num2str(JseqQerrorCnst)])

% we have already calculated SpatJfromseq, here we can show that it's the
% same as SpatJfromQcnst (I really struggle with the variables names here :D)
SpatJseqSpatJQCnsterr = sum(SpatJfromQCnst-SpatJfromseq,'all');
disp(['Error sum: (spatial J seq) - (spatial J quat) = ', num2str(SpatJseqSpatJQerr)])

%now we can map from spatial Jacobian to quat Jacobian
JQfromJseqCnst = spatial2quatCnst(SpatJfromseq,quat,2);
JQfromJseqCnsterr = sum(JQfromJseqCnst-JQCnstval,"all");
disp(['Error sum: (quat J mapped from J seq) - (quat J computed) = ', num2str(JQfromJseqCnsterr)])
% here the error is almost zero and it's unique, one joint external torques
% are represented by 3 'quaternion forces', not for as was in unconstrained
% case

%for polynomial fitting the second approach can be probably used. The
%moment arms would be mapped into the quaternion jacobians, then you would
%calculate the polynomials only in terms of q2,q3,q4 elements of
%quaternion. To calculate external forces you then to multiply the jacobian
%by invJtrans

function res = spatial2quat(jac,quat,Njoints)
    k = 1;
    kq = 1;
    for i=1:Njoints
        res(kq:kq+3,:) = 2 * G(quat(kq:kq+3))' * jac(k:k+2,:);
        k = k+3;
        kq = kq+4;
    end
end

function res = quat2spatial(jac,quat,Njoints)
    k = 1;
    kq = 1;
    for i=1:Njoints
        res(k:k+2,:) = 1/2 * G(quat(kq:kq+3)) * jac(kq:kq+3,:);
        k = k+3;
        kq = kq+4;
    end
end

function res = quat2spatialCnst(jac,quat,Njoints)
    k = 1;
    kq = 1;
    for i=1:Njoints
        res(k:k+2,:) = invEtrans(quat(kq:kq+3)) * jac(k:k+2,:);
        k = k+3;
        kq = kq+4;
    end
end

function res = spatial2quatCnst(jac,quat,Njoints)
    k = 1;
    kq = 1;
    for i=1:Njoints
        res(k:k+2,:) = Etrans(quat(kq:kq+3)) * jac(k:k+2,:);
        k = k+3;
        kq = kq+4;
    end
end

function res = eul2spatial(jac,phi,Njoints)
    k = 1;
    for i=1:Njoints
        res(k:k+2,:) = JG(phi(k:k+2))\jac(k:k+2,:);
        k = k+3;
    end
end

function res = spatial2eul(jac,phi,Njoints)
    k = 1;
    for i=1:Njoints
        res(k:k+2,:) = JG(phi(k:k+2))*jac(k:k+2,:);
        k = k+3;
    end
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

function res = E(Q)
    Q0 = Q(1);
    Q1 = Q(2);
    Q2 = Q(3);
    Q3 = Q(4);
    res = [-Q1, Q0, -Q3, Q2;
            -Q2,Q3, Q0, -Q1;
            -Q3, -Q2, Q1, Q0];
end

function res = JG(phi)
    % this is geometric jacobian of YZY sequence
    s2 = sin(phi(2));
    s3 = sin(phi(3));
    c2 = cos(phi(2));
    c3 = cos(phi(3));
    res = [s2*c3, c2, s2*s3; -s3, 0 ,c3; 0, 1, 0];
end

function res = T(quat)
    a = quat(1);
    b = quat(2);
    c = quat(3);
    d = quat(4);
    res = [-b/a, -c/a, -d/a;eye(3)];
end

function res = Etrans(quat)
    res = 2*G(quat)*T(quat);
    res = res';
end

function res = invEtrans(quat)
    % inverse of Jtrans - resulting matrix is actually this simple
    q1 = quat(1);
    q2 = quat(2);
    q3 = quat(3);
    q4 = quat(4);
    res = [ q1/2,  q4/2, -q3/2;
            -q4/2,  q1/2,  q2/2;
            q3/2, -q2/2,  q1/2];
end