function fo_python = fo_python(q,u,inputs,model)
l = model.l;
m = model.m;
g = model.g;
Ixx = model.Ixx;
Iyy = model.Iyy;
Izz = model.Izz;
c = model.c;
qm0 = q(1,:);
qm1 = q(2,:);
qm2 = q(3,:);
qm3 = q(4,:);
qm4 = q(5,:);
qm5 = q(6,:);
qm6 = q(7,:);
qm7 = q(8,:);
um1 = u(1,:);
um2 = u(2,:);
um3 = u(3,:);
um4 = u(4,:);
um5 = u(5,:);
um6 = u(6,:);
x0 = 0.5*um1;
x1 = 0.5*qm2;
x2 = 0.5*qm3;
x3 = 0.5*qm0;
x4 = 0.5*qm1;
x5 = 0.5*um4;
x6 = 0.5*qm6;
x7 = 0.5*qm7;
x8 = -qm5.*x5 - um5.*x6 - um6.*x7;
x9 = qm4.*x5 - um5.*x7 + um6.*x6;
x10 = 0.5*qm4;
x11 = 0.5*qm5;
x12 = qm7.*x5 + um5.*x10 - um6.*x11;
x13 = -qm6.*x5 + um5.*x11 + um6.*x10;
x14 = l.^2;
x15 = 2*qm0;
x16 = qm1.*x15 + 2*qm2.*qm3;
x17 = l.*m;
x18 = 2*qm4;
x19 = qm6.*x18;
x20 = 2*qm5.*qm7 - x19;
x21 = qm7.*x18;
x22 = 2*qm5.*qm6 - x21;
x23 = l/2;
x24 = x22.*x23;
x25 = qm5.*x18;
x26 = 2*qm7;
x27 = qm6.*x26 + x25;
x28 = qm5.^2;
x29 = qm7.^2;
x30 = -x29;
x31 = qm4.^2;
x32 = qm6.^2;
x33 = x31 - x32;
x34 = x28 + x30 + x33;
x35 = x23.*x34;
x36 = -x20.*x24 + x27.*x35;
x37 = l.*um1.^2 + l.*um2.^2;
x38 = m.*x37;
x39 = 2*qm5;
x40 = qm6.*x39 + x21;
x41 = -x28;
x42 = x30 + x31 + x32 + x41;
x43 = -l.*x34.*x42/2 + x24.*x40;
x44 = 2*qm1.*qm3 - qm2.*x15;
x45 = qm0.^2 - qm1.^2 - qm2.^2 + qm3.^2;
x46 = x16.*x40 + x20.*x45 + x34.*x44;
x47 = x16.*x42 + x22.*x44 + x27.*x45;
x48 = g.*m.*x47;
x49 = 2*qm6.*qm7 - x25;
x50 = um1.*x22 + um2.*x42 + um3.*x27 + um5;
x51 = um2.*x40;
x52 = um1.*x34;
x53 = um3.*x20 + um4 + x51 + x52;
x54 = x23.*x50.^2 + x23.*x53.^2;
x55 = qm7.*x39 + x19;
x56 = x50.*x53;
x57 = -Ixx.*x56 + Iyy.*x56 + Izz.*(-um4.*x50 + um5.*x53);
x58 = x29 + x33 + x41;
x59 = um1.*x55 + um2.*x49 + um3.*x58 + um6;
x60 = x53.*x59;
x61 = Ixx.*x60 + Iyy.*(um4.*x59 - um6.*x53) - Izz.*x60;
x62 = x13.*x18 + x26.*x8;
x63 = x39.*x9;
x64 = 2*qm6;
x65 = x12.*x64;
x66 = -x13.*x26 + x18.*x8;
x67 = -x23.*x60 - x23.*(um1.*(2*qm5.*x12 + 2*qm6.*x9 - x62) + um2.*(-x63 + x65 + x66) + um3.*(x12.*x26 + x13.*x64 + x18.*x9 + x39.*x8));
x68 = x50.*x59;
x69 = Ixx.*(-um5.*x59 + um6.*x50) - Iyy.*x68 + Izz.*x68;
x70 = -x23.*x68 + x23.*(um1.*(x63 - x65 + x66) + um2.*(x12.*x39 + x62 + x64.*x9) + um3.*(2*qm5.*x13 + 2*qm7.*x9 - x12.*x18 - x64.*x8));
x71 = m.*x70;
x72 = um1.*um3;
x73 = m.*x14;
x74 = x23.*x40;
x75 = -x20.*x23.*x42 + x27.*x74;
x76 = m.*x23;
x77 = g.*x47.*x76;
x78 = x70.*x76;
x79 = um3.*x73/2;
fo_python = [-qm1.*x0 - um2.*x1 - um3.*x2; qm0.*x0 - um2.*x2 + um3.*x1; qm3.*x0 + um2.*x3 - um3.*x4; -qm2.*x0 + um2.*x4 + um3.*x3; x8; x9; x12; x13; Iyy.*um2.*um3 - Izz.*um2.*um3 - c.*um1 + g.*l.*m.*x22.*x46/2 - 3*g.*x16.*x17/2 - l.*m.*um2.*um3.*x43 + l.*m.*x22.*x67/2 + 5*m.*um2.*um3.*x14/4 - x17.*x40.*x67 - x17.*x42.*x70 - x17.*x49.*x54 - x22.*x61 - x34.*x69 - x35.*x48 - x35.*x71 - x36.*x38 - x55.*x57; -Ixx.*x72 + Izz.*um1.*um3 - c.*um2 + g.*l.*m.*x42.*x46/2 + 3*g.*l.*m.*x44/2 + l.*m.*um1.*um3.*x43 + l.*m.*x22.*x70 + l.*m.*x34.*x67 + l.*m.*x42.*x67/2 + l.*m.*x54.*x55 - x38.*x75 - x40.*x69 - x42.*x61 - x48.*x74 - x49.*x57 - x71.*x74 - 5*x72.*x73/4; Ixx.*um1.*um2 - Iyy.*um1.*um2 - c.*um3 + g.*l.*m.*x27.*x46/2 - l.*m.*um1.*um3.*x36 - l.*m.*um2.*um3.*x75 + l.*m.*x27.*x67/2 - x20.*x69 - x20.*x77 - x20.*x78 - x27.*x61 - x57.*x58; -c.*um4 + m.*um1.*um3.*x14.*x22/2 + m.*um2.*um3.*x14.*x42/2 - x23.*x27.*x38 - x69 - x77 - x78; -c.*um5 + g.*l.*m.*x46/2 + l.*m.*x20.*x37/2 + l.*m.*x67/2 - x51.*x79 - x52.*x79 - x61; -c.*um6 - x57];
