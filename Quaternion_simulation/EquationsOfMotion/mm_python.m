function mm_python = mm_python(q,u,inputs,model)
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
x0 = qm5.^2;
x1 = qm6.^2;
x2 = -x1;
x3 = qm4.^2;
x4 = qm7.^2;
x5 = x3 - x4;
x6 = x0 + x2 + x5;
x7 = x6.^2;
x8 = 2*qm4;
x9 = qm7.*x8;
x10 = 2*qm5.*qm6 - x9;
x11 = x10.^2;
x12 = qm6.*x8;
x13 = 2*qm5;
x14 = qm7.*x13 + x12;
x15 = qm6.*x13 + x9;
x16 = l/2;
x17 = x15.*x16;
x18 = -x0;
x19 = x1 + x18 + x5;
x20 = -l.*x19.*x6/2 + x10.*x17;
x21 = l.^2;
x22 = x21/4;
x23 = x21/2;
x24 = x19.*x23;
x25 = x15.*x23;
x26 = -x10.*x25 + x21 + x24.*x6;
x27 = m.*x22;
x28 = Ixx + x27;
x29 = qm5.*x8;
x30 = 2*qm6.*qm7 - x29;
x31 = Izz.*x14;
x32 = Ixx.*x6;
x33 = Iyy.*x10;
x34 = x10.*x22;
x35 = x22.*x6;
x36 = m.*(x15.*x35 + x19.*x34) + x15.*x32 + x19.*x33 + x30.*x31;
x37 = 2*qm6.*qm7 + x29;
x38 = 2*qm5.*qm7 - x12;
x39 = x34.*x37 + x35.*x38;
x40 = x18 + x2 + x3 + x4;
x41 = x31.*x40 + x32.*x38 + x33.*x37;
x42 = m.*(x24 + x35) + x32;
x43 = m.*(-x25 + x34) + x33;
x44 = x15.^2;
x45 = x19.^2;
x46 = Iyy + x27;
x47 = x19.*x22;
x48 = x15.*x22.*x38 + x37.*x47;
x49 = Ixx.*x15;
x50 = Iyy.*x19;
x51 = Izz.*x30;
x52 = x37.*x50 + x38.*x49 + x40.*x51;
x53 = x10.*x23;
x54 = m.*(x15.*x21/4 - x53) + x49;
x55 = x23.*x6;
x56 = m.*(x47 + x55) + x50;
x57 = x38.^2;
x58 = x37.^2;
x59 = Ixx.*x38 + x27.*x38;
x60 = Iyy.*x37 + x27.*x37;
x61 = Izz.*x40;
mm_python = [1 0 0 0 0 0 0 0 0 0 0 0 0 0; 0 1 0 0 0 0 0 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 0 0 0 0; 0 0 0 0 1 0 0 0 0 0 0 0 0 0; 0 0 0 0 0 1 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 1 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 Ixx.*x7 + Iyy.*x11 + Izz.*x14.^2 + m.*(-l.*x20 + x11.*x22 + x22.*x7 + x26) + x28 x36 m.*(l.*(l.*x19.*x38/2 - x17.*x37) + x39) + x41 x42 x43 x31; 0 0 0 0 0 0 0 0 x36 Ixx.*x44 + Iyy.*x45 + Izz.*x30.^2 + m.*(-l.*x20 + x22.*x44 + x22.*x45 + x26) + x46 m.*(-l.*(l.*x10.*x38/2 - x16.*x37.*x6) + x48) + x52 x54 x56 x51; 0 0 0 0 0 0 0 0 m.*(x24.*x38 - x25.*x37 + x39) + x41 m.*(x37.*x55 - x38.*x53 + x48) + x52 Ixx.*x57 + Iyy.*x58 + Izz.*x40.^2 + Izz + m.*(x22.*x57 + x22.*x58) x59 x60 x61; 0 0 0 0 0 0 0 0 x42 x54 x59 x28 0 0; 0 0 0 0 0 0 0 0 x43 x56 x60 0 x46 0; 0 0 0 0 0 0 0 0 x31 x51 x61 0 0 Izz];
