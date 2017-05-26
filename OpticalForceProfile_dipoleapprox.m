%% Initialization of the workspace
clear all;
close all;
clc;

%% Parameters
dx = 5e-9;
dy = 5e-9;
dz = 5e-9;
Lx = 1500e-9;
Ly = 1500e-9;
Lz = 1500e-9;

L = 32;

%% Dipole

ep = 2.25;
a = 250e-9;
lambda0 = 1064e-9;
alpharc = InducedDipole.polarizability('radiative correction',a,ep,'lambda0',lambda0);
id = InducedDipole(alpharc,lambda0);


%% Setup in POM 0.26
% Beam waist = 10.7e-3 m
% f = 200 mm (tube length for infinity corrected) / 60 (obj magnification)
% Pmax = 60/100 (obj transmission) * 64/100 (AOD transmission) * 2 W

%% Beam & objective

w0 = 10.7e-3;
Ex0 = 1;
Ey0 = 0;
nm = 1.33;
NA = 1.20;
f = 200e-3/60;
R = f*NA/nm;
f = nm*R/NA;

%% Focused field

Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr,'lambda0',1064e-9);
bg = bg.normalize(60*64*2/10000);
ef = EFieldFocus(bg,f,'lambda0',1064e-9,'er',1.33^2);

%% Field

[X,Y,Z] = meshgrid([-Lx:dx:Lx],[-Ly:dy:Ly],0);
E= ef.E(Point(X,Y,Z),'er',1.33.^2);
I = .5*PhysConst.c0/nm*(nm^2*PhysConst.e0)*norm(E).^2;
%%
figure(1)
surf(X,Y,I)
shading flat
axis auto
%% Forces

[X,Y,Z] = meshgrid([-Lx:10*dx:Lx],[-Ly:10*dy:Ly],0);
r = Point(X,Y,Z);
[F,Fgrad,Fscat,Fsc] = id.force(r,ef);

%%
figure(2)
subplot(1,2,1)
hold on
plot(F.X((length(X)-1)/2,:),F.Vx((length(X)-1)/2,:))
hold off

pp = csapi(F.X((length(X)-1)/2,:),F.Vx((length(X)-1)/2,:));
ipp=fnint(pp);
dpp=fnder(pp);

subplot(1,2,2)
hold on
plot(F.X((length(X)-1)/2,:),fnval(dpp,F.X((length(X)-1)/2,:)))
hold off
