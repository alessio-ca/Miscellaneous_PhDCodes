%% Initialization of the workspace
clear all;
close all;
clc;

%% Parameters
dx = 1e-9;
dy = 1e-9;
dz = 1e-9;
Lx = 2000e-9;
Ly = 2000e-9;
Lz = 20000e-9;

L = 32;

%% Dipole

ep = 1.60^2;
a = 0.527e-6/2;
lambda0 = 1064e-9;
nm = 1.33;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0,'em',nm^2);
id = InducedDipole(alpharc,lambda0,'er',nm.^2);


%% Setup in POM 0.26
% Beam waist = 10.7e-3 m
% f = 200 mm (tube length for infinity corrected) / 60 (obj magnification)
% Pmax = 60/100 (obj transmission) * 64/100 (AOD transmission) * 2 W

%% Beam & objective

w0 = 10.7e-3;
Ex0 = 1;
Ey0 = 0;
NA = 1.20;
f = 200e-3/60;
R = f*NA/nm;
f = nm*R/NA;

%% Focused field

Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr,'lambda0',lambda0,'er',nm^2);
bg = bg.normalize(60*64*2/(100*100));
ef = EFieldFocus(bg,f,'lambda0',1064e-9,'er',nm^2);

%% Field

[X,Y,Z] = meshgrid(-Lx:dx:Lx,-Ly:dy:Ly,0);
E= ef.E(Point(X,Y,Z),'lambda0',lambda0,'er',nm^2);
I = .5*PhysConst.c0/nm*(nm^2*PhysConst.e0)*norm(E).^2;
%%
figure(1)
surf(X,Y,I)
shading flat
axis auto
%% Forces

[X,Y,Z] = meshgrid(-Lx:10*dx:Lx,-Ly:10*dy:Ly,0);
r = Point(X,Y,Z);
[F,Fgrad,Fscat,Fsc] = id.force(r,ef);

%%
figure(2)
subplot(1,2,1)
plot(F.X((length(X)-1)/2,:),F.Vx((length(X)-1)/2,:))
hold on
plot(F.Y(:,(length(X)-1)/2),F.Vy(:,(length(X)-1)/2))
hold off

pp_x = csapi(F.X((length(X)-1)/2,:),F.Vx((length(X)-1)/2,:));
pp_y = csapi(F.Y(:,(length(X)-1)/2),F.Vy(:,(length(X)-1)/2));

ipp_x=fnint(pp_x);
ipp_y=fnint(pp_y);


subplot(1,2,2)
plot(F.X((length(X)-1)/2,:),-fnval(ipp_x,F.X((length(X)-1)/2,:)))
hold on
plot(F.Y(:,(length(X)-1)/2),-fnval(ipp_y,F.Y(:,(length(X)-1)/2)))
hold off
