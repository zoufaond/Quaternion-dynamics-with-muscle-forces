clear all
clc
addpath MuscleForces\
load("model_struct.mat");

% create symbolic coordinates
q = sym('q',[1 8],'real');
phi = sym('phi',[1 6],'real');


genEq = 0; % 1 .. generate functions of muscle forces
JQsym = genFm_quat(genEq,model); % unconstrained muscle lengths jacobian
JQCnstsym = genFm_quat_Cnst(genEq,model); % constrained muscle lengths jacobian
JEulsym = genFM_seq(genEq,model); % muscle lengths jacobian for YZY sequence of rotations
disp('Muscle Jacobians derived')

% anonymous functions
JQ = matlabFunction(JQsym,'Vars',{q});
JQCnst = matlabFunction(JQCnstsym,'Vars',{q});
JEul = matlabFunction(JEulsym,'Vars',{phi});

%%
% run forward dynamics in simulink (initFunction )
out = sim("double_3D_pend_Quat.slx");
plot_results(out)
%%
% out .. output from simulink
% out.time .. timespan of simulation
% out.q_1/q_2 .. resulting quaternions
% out.yzy_1/yzy_2 .. resulting yzy angles

time = out.tout';
quat_1 = reshape(out.q_1.Data,[4,length(time)])';
eul_1 = reshape(out.yzy_1.Data,[3,length(time)])';
quat_2 = reshape(out.q_2.Data,[4,length(time)])';
eul_2 = reshape(out.yzy_2.Data,[3,length(time)])';

eul_full = [eul_1,eul_2];
quat_full = [quat_1,quat_2];


JQfromJEulCnst_mus = zeros(6,length(time));
JQCnstval_mus = zeros(6,length(time));
JQfromJEul_mus = zeros(8,length(time));
JQval_mus = zeros(8,length(time));
JEulfromJQ = zeros(6,length(time));
JEulval_mus = zeros(6,length(time));

% analyze element 2
imus = 2;

for i = 1:length(time)
    % get current yzy and quat values
    phiVal = eul_full(i,:);
    quatVal = quat_full(i,:);

    % get the i-th muscle length Jacobian in Euler coordinates (moment arms)
    JEulval = JEul(phiVal);
    JEulval_mus(:,i) = JEulval(:,imus);

    % map Jacobian from Euler coordinates to spatial coordinates using
    % geometric Jacobian
    SpatJfromJEul = eul2spatial(JEulval,phiVal,2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                             SECTION 2.2                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the i-th muscle length Jacobian in quaternions (w.r.t. all quaternion elements)
    JQval = JQ(quatVal);
    JQval_mus(:,i) = JQval(:,imus);

    % map Jacobian from quaternion to spatial coordinates
    spatJfromJQ = quat2spatial(JQval_mus(:,i),quatVal,2);

    % map the quaternion Jacobian from spatial coordinates to Euler
    % coordinates (to see what moment arms the quaternion Jacobian represents)
    JEulfromJQ(:,i) = spatial2eul(spatJfromJQ,phiVal,2);

    % the moment arms back to quaternion Jacobian (this doesn't give us correct muscle length Jacobian)
    JQfromJEul = spatial2quat(SpatJfromJEul,quatVal,2);
    JQfromJEul_mus(:,i) = JQfromJEul(:,imus);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                         END SECTION 2.2                             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                       SECTION 2.3                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get the i-th muscle length Jacobian in quaternion coordinates (approach based on section Constraint)
    quatVal_pos = quatVal;
    if quatVal(1) < 0
        quatVal_pos(1:4) = quatVal_pos(1:4)*(-1);
    end

    if quatVal(5) < 0
        quatVal_pos(5:8) = quatVal_pos(5:8)*(-1);
    end
        
    JQCnstval = JQCnst(quatVal_pos);
    JQCnstval_mus(:,i) = JQCnstval(:,imus);
    
    % map the Euler muscle length Jacobian from spatial coordinates to
    % quaternion coordinates (this gives us Jacobian values that we calculated analytically)
    JQfromJEulCnst = spatial2quatCnst(SpatJfromJEul,quatVal_pos,2);
    JQfromJEulCnst_mus(:,i) = JQfromJEulCnst(:,imus);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   END SECTION 2.3                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end


% Moment arms
jnts = {'R_y^{U}','R_z^{U}','R_{yy}^{U}','R_y^{L}','R_z^{L}','R_{x}^{L}'};

figure
for j = 1:6
    subplot(2,3,j)
    plot(time(1:end), JEulval_mus(j,1:end),'b',time(1:end),JEulfromJQ(j,1:end),'r--','LineWidth',1.7)
    % title('$$'+jnts{j}+'$$','Interpreter','latex','FontSize',12)
    if j == 4
        title({newline,['\indent $',jnts{j},'$']},'Interpreter','latex','FontSize',15)
    else
        title({newline,['$',jnts{j},'$']},'Interpreter','latex','FontSize',15)
    end
    xlabel('time [s]')
    ylabel({'Moment';'arm [m]'})
    axis([-inf,inf,-inf,inf])
    % ylabel('$$\textbf{R}^*_{Q_'+string(j)+'}$$','Interpreter','latex','FontSize',12)
end
% sgtitle('Jacobian in Euler coordinates (moment arms)')
fig = gcf;
fig.Position(3:4)=[800,350];
Lgnd = legend({'Euler angles ML Jacobian',['Quaternion ML Jacobian mapped ' newline 'to Euler angles ML Jacobian']});
Lgnd.Position(1) = 0.0;
Lgnd.Position(2) = 0.365;
Lgnd.FontSize = 8;
Lgnd.IconColumnWidth = 17;
text1 = annotation('textbox', [0.3, 0.44, 0.1, 0.1], 'string','Upper joint' ,'FontSize',13,'EdgeColor','none');
text1.Rotation = 0;
text2 = annotation('textbox', [0.3, 0.355, 0.1, 0.1], 'string','Lower joint' ,'FontSize',13,'EdgeColor','none');
text2.Rotation = 0;
annotation('line',[0.25,0.9],[0.45,0.45],'LineWidth',2)
% 
exportgraphics(fig,'moment_arms_new.jpg','Resolution',600);



% unconstrained approach
figure
qcoord = 0;
body = {'^{U}','^{L}'};
ibody = 1;
for j = 1:8
    subplot(2,4,j)
    plot(time(1:end), JQfromJEul_mus(j,(1:end)),'b',time(1:end),JQval_mus(j,(1:end)),'r--','LineWidth',1.7)
    ylabel({'$$\textbf{R}_{Q_'+string(qcoord)+body(ibody)+'}$$','$$[-]$$'},'rotation',0,'Interpreter','latex','FontSize',15)
    xlabel('time [s]')
    axis([-inf,inf,-inf,inf])

    if j == 4 || j == 5 || j == 6
        title(newline)
    end

    qcoord = qcoord +1;
    if j == 4
        qcoord = 0;
        ibody = ibody+1;
    end
end
fig = gcf;
fig.Position(3:4)=[800,350];
Lgnd = legend({['Euler ML Jacobian mapped' newline 'to quaternion coordinates'],'Quaternion ML Jacobian'});
Lgnd.Position(1) = 0.01;
Lgnd.Position(2) = 0.375;
text1 = annotation('textbox', [0.27, 0.45, 0.1, 0.1], 'string','Upper joint' ,'FontSize',13,'EdgeColor','none');
text1.Rotation = 0;
text2 = annotation('textbox', [0.27, 0.375, 0.1, 0.1], 'string','Lower joint' ,'FontSize',13,'EdgeColor','none');
text2.Rotation = 0;
annotation('line',[0.23,0.9],[0.46,0.46],'LineWidth',2)
% % 
% % exportgraphics(fig,'JQ_comparison.png','Resolution',600);

% Constrained approach
figure
qcoord = 1;
body = {'^{U}','^{L}'};
ibody = 1;
for j = 1:6
    subplot(2,3,j)
    plot(time(1:end), JQfromJEulCnst_mus(j,1:end),'b',time(1:end),JQCnstval_mus(j,1:end),'r--','LineWidth',1.7)
    % ylabel({newline,'$$\frac{\partial l_m(\vec{Q}^*)}{\partial Q_'+string(qcoord)+body(ibody)+'}$$'},'Rotation',0,'Interpreter','latex','FontSize',12)
    ylabel({'$$\textbf{R}^*_{Q_'+string(qcoord)+body(ibody)+'}$$','$$[-]$$'},'rotation',0,'Interpreter','latex','FontSize',15)

    if j == 4 || j == 5 || j == 6
        title(newline)
    end

    if j == 1
        xlabel('time [s]','Position',[3,-1.6])
    else
        xlabel('time [s]')
    end
    
    axis([-inf,inf,-inf,inf])
    
    qcoord = qcoord +1;
    if j == 3
        qcoord = 1;
        ibody = ibody+1;
    end

end
fig = gcf;
fig.Position(3:4)=[800,350];
Lgnd = legend({['Euler ML Jacobian mapped ' newline ' to quaternion coordinates'],'Quaternion ML Jacobian'});
Lgnd.Position(1) = 0.0;
Lgnd.Position(2) = 0.41;
text1 = annotation('textbox', [0.3, 0.47, 0.1, 0.1], 'string','Upper joint' ,'FontSize',13,'EdgeColor','none');
text1.Rotation = 0;
text2 = annotation('textbox', [0.3, 0.375, 0.1, 0.1], 'string','Lower joint' ,'FontSize',13,'EdgeColor','none');
text2.Rotation = 0;
annotation('line',[0.23,0.9],[0.475,0.475],'LineWidth',2)

% exportgraphics(fig,'JQCnst_comparison.png','Resolution',600);

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
    joints_seq = {'YZY','YZX'};
    for i=1:Njoints
        res(k:k+2,:) = JG(phi(k:k+2),joints_seq{i})\jac(k:k+2,:);
        k = k+3;
    end
end

function res = spatial2eul(jac,phi,Njoints)
    k = 1;
    joints_seq = {'YZY','YZX'};
    for i=1:Njoints
        res(k:k+2,:) = JG(phi(k:k+2),joints_seq{i})*jac(k:k+2,:);
        k = k+3;
    end
end

function res = JG(phi,seq)
    % this is geometric jacobian of YZY sequence
    s2 = sin(phi(2));
    s3 = sin(phi(3));
    c2 = cos(phi(2));
    c3 = cos(phi(3));
    if strcmp(seq,'YZY')
        res = [s2*c3, c2, s2*s3;
              -s3, 0 ,c3; 
               0, 1, 0];
    elseif strcmp(seq,'YZX')
        res = [s2,c2*c3,-s3*c2; 
               0, s3,c3; 
               1, 0, 0];
    end

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

function plot_results(out)
% disp('Number simulation steps:')
% disp(length(out.tout))

line_colors = [0, 0, 1, 1; 1, 0, 0, 1;0, 1, 0, 1; 0.4660, 0.6740, 0.1880, 1];
line_colors_simulink = [0, 0, 1, 1; 1, 0, 0, 1;0, 1, 0, 1; 0.4660, 0.6740, 0.1880, 1];
time = out.tout;
time_interp = linspace(0,time(end),30);
quat0_1 = out.q1_eq.Data(1,:);
quat1_1 = out.q1_eq.Data(2,:);
quat2_1 = out.q1_eq.Data(3,:);
quat3_1 = out.q1_eq.Data(4,:);
quat0_1_simulink = out.q_1.Data(1,:);
quat1_1_simulink = out.q_1.Data(2,:);
quat2_1_simulink = out.q_1.Data(3,:);
quat3_1_simulink = out.q_1.Data(4,:);
quat0_1_simulink_interp1 = interp1(time,quat0_1_simulink,time_interp);
quat1_1_simulink_interp1 = interp1(time,quat1_1_simulink,time_interp);
quat2_1_simulink_interp1 = interp1(time,quat2_1_simulink,time_interp);
quat3_1_simulink_interp1 = interp1(time,quat3_1_simulink,time_interp);

quat0_2 = out.q2_eq.Data(1,:);
quat1_2 = out.q2_eq.Data(2,:);
quat2_2 = out.q2_eq.Data(3,:);
quat3_2 = out.q2_eq.Data(4,:);
quat0_2_simulink = out.q_2.Data(1,:);
quat1_2_simulink = out.q_2.Data(2,:);
quat2_2_simulink = out.q_2.Data(3,:);
quat3_2_simulink = out.q_2.Data(4,:);
quat0_2_simulink_interp1 = interp1(time,quat0_2_simulink,time_interp);
quat1_2_simulink_interp1 = interp1(time,quat1_2_simulink,time_interp);
quat2_2_simulink_interp1 = interp1(time,quat2_2_simulink,time_interp);
quat3_2_simulink_interp1 = interp1(time,quat3_2_simulink,time_interp);

figure
end_point = 0;
subplot(2,1,1)
plot(time(1:end-end_point),quat0_1(1:end-end_point),'LineWidth',1.5,'Color',line_colors_simulink(1,:)); hold on
% plot(time_interp(1:end-end_point),quat0_1_simulink_interp1(1:end-end_point),'Color',line_colors_simulink(1,:),'Marker','o','LineStyle', 'none','MarkerSize',8); hold on
plot(time(1:end-end_point),quat1_1(1:end-end_point),'LineWidth',1.5,'Color',line_colors_simulink(2,:)); hold on
% plot(time_interp(1:end-end_point),quat1_1_simulink_interp1(1:end-end_point),'Color',line_colors_simulink(2,:),'Marker','o','LineStyle', 'none','MarkerSize',8); hold on
plot(time(1:end-end_point),quat2_1(1:end-end_point),'LineWidth',1.5,'Color',line_colors_simulink(3,:)); hold on
% plot(time_interp(1:end-end_point),quat2_1_simulink_interp1(1:end-end_point),'Color',line_colors_simulink(3,:),'Marker','o','LineStyle', 'none','MarkerSize',8); hold on
plot(time(1:end-end_point),quat3_1(1:end-end_point),'LineWidth',1.5,'Color',line_colors_simulink(4,:)); hold on
% plot(time_interp(1:end-end_point),quat3_1_simulink_interp1(1:end-end_point),'Color',line_colors_simulink(4,:),'Marker','o','LineStyle', 'none','MarkerSize',8); hold on

hold off
xlim([0 4])
ylim([-1 1])
% leg = legend({'Q_0','Q_0^{Sim}','Q_1','Q_1^{Sim}','Q_2','Q_2^{Sim}','Q_3','Q_3^{Sim}'},'NumColumns',4,'FontSize',8,'Location','northeast');
leg = legend({'Q_0^{Upper}','Q_1^{Upper}','Q_2^{Upper}','Q_3^{Upper}'},'NumColumns',2,'FontSize',8,'Location','northeast');
leg.IconColumnWidth = 10;
leg.Position(1) = [0.6];
leg.Position(2) = [0.58];
title('Quaternion-based model - Upper joint')
xlabel('Time [s]')
ylabel('[-]')

subplot(2,1,2)
plot(time(1:end-end_point),quat0_2(1:end-end_point),'LineWidth',1.5,'Color',line_colors_simulink(1,:)); hold on
% plot(time_interp(1:end-end_point),quat0_2_simulink_interp1(1:end-end_point),'Color',line_colors_simulink(1,:),'Marker','o','LineStyle', 'none','MarkerSize',8); hold on
plot(time(1:end-end_point),quat1_2(1:end-end_point),'LineWidth',1.5,'Color',line_colors_simulink(2,:)); hold on
% plot(time_interp(1:end-end_point),quat1_2_simulink_interp1(1:end-end_point),'Color',line_colors_simulink(2,:),'Marker','o','LineStyle', 'none','MarkerSize',8); hold on
plot(time(1:end-end_point),quat2_2(1:end-end_point),'LineWidth',1.5,'Color',line_colors_simulink(3,:)); hold on
% plot(time_interp(1:end-end_point),quat2_2_simulink_interp1(1:end-end_point),'Color',line_colors_simulink(3,:),'Marker','o','LineStyle', 'none','MarkerSize',8); hold on
plot(time(1:end-end_point),quat3_2(1:end-end_point),'LineWidth',1.5,'Color',line_colors_simulink(4,:)); hold on
% plot(time_interp(1:end-end_point),quat3_2_simulink_interp1(1:end-end_point),'Color',line_colors_simulink(4,:),'Marker','o','LineStyle', 'none','MarkerSize',8); hold on

hold off
xlim([0 4])
ylim([-1 1])

leg = legend({'Q_0^{Lower}','Q_1^{Lower}','Q_2^{Lower}','Q_3^{Lower}'},'NumColumns',2,'FontSize',8,'Location','northeast');
% leg = legend({'Q_0','Q_0^{Sim}','Q_1','Q_1^{Sim}','Q_2','Q_2^{Sim}','Q_3','Q_3^{Sim}'},'NumColumns',4,'FontSize',8,'Location','northeast');

leg.IconColumnWidth = 10;
leg.Position(1) = [0.6];
leg.Position(2) = [0.1];
title('Quaternion-based model - Lower joint')
xlabel('Time [s]')
ylabel('[-]')
end