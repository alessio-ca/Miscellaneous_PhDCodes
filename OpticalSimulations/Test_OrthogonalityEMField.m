%% Test orthogonality EM field by EM.Field classes

%% Plane wave
E0 = 1;
k = Vector(0,0,0,0,0,1);
e = ComplexVector(0,0,0,1,0,0);
ef = EFieldPlaneWave(E0,k,e);


[x,y,z] = meshgrid(-1e-6:1e-7:1e-6,-1e-6:5e-8:1e-6,1e-7);
r = Point(x,y,z);
E = ef.E(r);
B = ef.B(r);

figure(1)

subplot(1,2,1)
title('Real part')
hold on
E.real().plot('scale',[1e+9 50/max(max(max(abs(norm(E)))))],'color','k')
B.real().plot('scale',[1e+9 50/max(max(max(abs(norm(B)))))],'color','b')
hold off
axis equal
view(3)
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')

subplot(1,2,2)
title('imaginary part')
hold on
E.imag().plot('scale',[1e+9 50/max(max(max(abs(norm(E)))))],'color','r')
B.imag().plot('scale',[1e+9 50/max(max(max(abs(norm(B)))))],'color','g')
hold off
axis equal
view(3)
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')

%% Beam wave

w0 = 1e-6;
Ex0 = 1;
Ey0 = 0;
R = 4e-6;
Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);
f = 1.075*R;
ef = EFieldFocus(bg,f);

[x,y,z] = meshgrid(-2e-6:0.5e-7:2e-6,-2e-6:0.5e-7:2e-6,0);
r = Point(x,y,z);

E=ef.E(r);
B=ef.B(r);

figure(2)
subplot(2,2,1)
title('Real part')
hold on
E.real().plot('scale',[1e+9 50/max(max(max(abs(norm(E)))))],'color','k')
B.real().plot('scale',[1e+9 50/max(max(max(abs(norm(B)))))],'color','b')
hold off
axis equal
view(3)
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')

subplot(2,2,2)
title('imaginary part')
hold on
E.imag().plot('scale',[1e+9 50/max(max(max(abs(norm(E)))))],'color','r')
B.imag().plot('scale',[1e+9 50/max(max(max(abs(norm(B)))))],'color','g')
hold off
axis equal
view(3)
xlabel('x [nm]')
ylabel('y [nm]')
zlabel('z [nm]')


subplot(2,2,3)
title('Ex, z=0')
surf(x,y,sqrt(E.Vx.*conj(E.Vx)))
shading interp
%axis equal
view(3)
xlabel('x [nm]')
ylabel('y [nm]')

subplot(2,2,4)
title('Ez, z=0')
surf(x,y,sqrt(E.Vz.*conj(E.Vz)))
shading interp
%axis equal
view(3)
xlabel('x [nm]')
ylabel('y [nm]')

z=y;
y=x;
x=zeros(size(x));

r = Point(x,y,z);

E=ef.E(r);
B=ef.B(r);

figure(3)
surf(y,z,real(E.Vx))
shading interp
axis equal
view(2)
xlabel('x [nm]')
ylabel('y [nm]')

