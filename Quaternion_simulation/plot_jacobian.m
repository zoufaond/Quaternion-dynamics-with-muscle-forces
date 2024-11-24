clear all
clc
addpath MuscleForces\
genEq = 0; % don't generate external torques equations, we already have them (it also takes quite a long time)

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

%%
out = sim("double_3D_pend_Quat.slx");
%%

time = out.tout;
quat_1 = reshape(out.q_1.Data,[4,length(time)])';
quat_1 = positive_quat(quat_1);
eul_1 = reshape(out.yzy_1.Data,[3,length(time)])';
quat_2 = reshape(out.q_2.Data,[4,length(time)])';
quat_2 = positive_quat(quat_2);
eul_2 = reshape(out.yzy_2.Data,[3,length(time)])';

eul_full = [eul_1,eul_2];
quat_full = [quat_1,quat_2];
JQfromJseqCnst_mus1 = zeros(6,length(time));
JQCnstval_mus1 = zeros(6,length(time));
JQfromJseq_mus1 = zeros(8,length(time));
JQval_mus1 = zeros(8,length(time));

for i = 1:length(time)
    time_step = i;
    phiVal = eul_full(time_step,:);
    quatVal = quat_full(time_step,:);
    
    Jseqval = Jseq(phiVal);
    SpatJfromseq = eul2spatial(Jseqval,phiVal,2);
    JQfromJseqCnst = spatial2quatCnst(SpatJfromseq,quatVal,2);
    JQfromJseqCnst_mus1(:,i) = JQfromJseqCnst(:,1);
    JQCnstval = JQCnst(quatVal);
    JQCnstval_mus1(:,i) = JQCnstval(:,1);
    
    JQval = JQ(quatVal);
    JQval_mus1(:,i) = JQval(:,1);
    JQfromJseq = spatial2quat(SpatJfromseq,quatVal,2);
    JQfromJseq_mus1(:,i) = JQfromJseq(:,1);
end


%%
figure
tiledlayout(2,12)
% tiledlayout(12,2,"vertical")
gap = 1;

for j = 1:4
    nexttile([1,3])
    plot(time(1:gap:end), JQfromJseq_mus1(j,(1:gap:end)),'b',time(1:gap:end),JQval_mus1(j,(1:gap:end)),'r--','LineWidth',1.7)
    title('$$\frac{\partial l_m(\vec{Q})}{\partial Q_'+string(j-1)+'}$$','Interpreter','latex','FontSize',12)
    xlabel('time [s]')
    % ylabel('$$\textbf{R}_{Q_'+string(j-1)+'}$$','Interpreter','latex','FontSize',12)
end

for j = 1:3
    nexttile([1,4])
    plot(time(1:gap:end), JQfromJseqCnst_mus1(j,1:gap:end),'b',time(1:gap:end),JQCnstval_mus1(j,1:gap:end),'r--','LineWidth',1.7)
    title('$$\frac{\partial l_m(\vec{Q}^*)}{\partial Q_'+string(j)+'}$$','Interpreter','latex','FontSize',12)
    xlabel('time [s]')
    % ylabel('$$\textbf{R}^*_{Q_'+string(j)+'}$$','Interpreter','latex','FontSize',12)
    
end


% sgtitle({'Comparison of quaternion muscle length jacobian calculated analytically ','and mapped from Euler coordinate system'}, 'Interpreter', 'none')
fig = gcf;
fig.Position(3) = fig.Position(3) + 250;
Lgnd = legend('Mapped from Euler','Analytical derivatives');
Lgnd.Position(1) = 0.3;
Lgnd.Position(2) = 0.425;

text1 = annotation('textbox', [0.08, 0.55, 0.1, 0.1], 'string', {'$$\textbf{R}_Q = 2 \textbf{G}^T \textbf{R}_S$$'},'Interpreter','latex','FontSize',15);
text1.Rotation = 90;

text3 = annotation('textbox', [0.0, 0.45, 0.1, 0.1], 'string',{'Mapping','used:'},'FontSize',15,'EdgeColor','none');


text2 = annotation('textbox', [0.08, 0.06, 0.1, 0.1], 'string', {'$$\textbf{R}_{Q}^* = \textbf{E}^T \textbf{R}_{S}$$'},'Interpreter','latex','FontSize',15);
text2.Rotation = 90;

exportgraphics(fig,'Comparison_jacobian.png','Resolution',600);

function res = spatial2quat(jac,quat,Njoints)
    k = 1;
    kq = 1;
    for i=1:Njoints
        res(kq:kq+3,:) = 2 * G(quat(kq:kq+3))' * jac(k:k+2,:);
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

function res = positive_quat(quat)
    numdata = length(quat(:,1));

    for i = 1:numdata
        if quat(i,1) < 0
            quat(i,:) = quat(i,:) * -1;
        end
    end

    res = quat;
end