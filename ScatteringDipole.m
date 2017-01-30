%% DEFINITION OF BEAM

w0 = 5e-3;
Ex0 = 1;
Ey0 = 0;
R = 5e-3;
Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr);


f = 1.075*R;
ef = EFieldFocus(bg,f);
%% DEFINITION OF INDUCEDDIPOLE
ep = 2.25;
a = 75e-9;
lambda0 = 633e-9;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0);
id = InducedDipole(alpharc,lambda0);
%%
close all

%% Plot of lateral focal field
[x,z,y] = meshgrid(-5e-7:1e-8:5e-7,-5e-7:1e-8:5e-7,0);
r = Point(x,y,z);
E=ef.E(r);
I = 0.5/PhysConst.Z0*sqrt(ef.er/ef.mr)*(E.*conj(E));
h=pcolor(x*1e6,z*1e6,I);
colormap bone
shading interp
colormap(flipud(colormap))

%% Plot of scattering and gradient force
[x,z,y] = meshgrid(-5e-7:1e-7:5e-7,-5e-7:1e-7:5e-7,0);
r = Point(x,y,z);
[F,Fgrad,Fscat,Fsc] = id.force(r,ef);
S1=1e6;
S2=1e-1/max(max(abs(norm(F))));
hold on
quiver(S1*Fgrad.X,S1*Fgrad.Z,real(S2*Fgrad.Vx),real(S2*Fgrad.Vz),0,'color','g')
quiver(S1*Fscat.X,S1*Fscat.Z,real(S2*Fscat.Vx),real(S2*Fscat.Vz),0,'color','r')
xlabel('x [nm]')
ylabel('z [nm]')
