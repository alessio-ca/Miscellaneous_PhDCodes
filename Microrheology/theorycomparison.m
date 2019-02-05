% theorycomparison.m - Comparison between dipole approximation, geometrical optics and exact electromagnetic theory
%
% Calculation of the transverse trap stiffness produced by a focused
% Gaussian beam on a microscopic sphere. Results are obtained using:
% (1) the dipole approximation (DA), 
% (2) the geometrical optics approximation (GO) and
% (3) the exact electromagnetic theory (EMT).
% This code reproduces Fig. 5.1.
%
% Note that the running time of this code (EMT) is quite long (several days).
% Therefore, the values are pre-calculated and saved 
% in the files DA.mat, GO.mat and EMT.mat.
%
% See also BEAMGAUSS,EFIELDFOCUS, MIEPARTICLE.

%   Author: Giovanni Volpe
%   Revision: 1.0.0  
%   Date: 2015/01/01

%% Initialization of the workspace
clear all;
close all;
clc;

%% Particle Parameters
nm = 1.33;  % medium refractive index
em = real(nm^2);  % medium relative dielectric constant

np = 1.50;  % particle refractive index
ep = real(np^2);  % particle relative dielectric constant

lambda0 = 632e-9;  % light vacuum wavelength
k0 = 2*pi/lambda0;  % light vacuum wavenumber

kaDA = [.1:.1:3.5];  % sphere size parameter (DA)
aDA = kaDA/(nm*k0);  % sphere radius (DA) [m]

kaGO = [2:.1:5];  % sphere size parameter (GO)
aGO = kaGO/(nm*k0);  % sphere radius (GO) [m]

ka = [.1:.1:5];  % sphere size parameter (EMT)
a = ka/(nm*k0);  % sphere radius (EMT) [m]

Dx = 100e-9;  % particle displacement [m]
C = Point(Dx,0,0);

%% Beam
w0 = 100e-3;
Ex0 = 1;
Ey0 = 0;
R = 100e-3;
Nphi = 28;
Nr = 18;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr,'lambda0',lambda0);
bg = bg.normalize(1);

%% Objective
NA = 1.20;
f = R*nm/NA;

er1 = 1.00;
er2 = nm^2;
mr1 = 1.00;
mr2 = 1.00;

%% CALCULATION - DA
fileDA = 'DA.mat';
if ~exist(fileDA)
    
    ef = EFieldFocus(bg,f,'lambda0',lambda0,'er',em);
    
    FxDA = zeros(size(aDA));
    FyDA = zeros(size(aDA));
    FzDA = zeros(size(aDA));
    for i = 1:1:length(aDA)
        
        alpharc = InducedDipole.polarizability('Clausius-Mossotti',aDA(i),ep,'lambda0',lambda0,'em',em);
        id = InducedDipole(alpharc,lambda0);
        
        F = id.force(C,ef)
        FxDA(i) = F.Vx;
        FyDA(i) = F.Vy;
        FzDA(i) = F.Vz;
        
        figure(1)
        cla
        hold on
        plot(aDA*1e+6,-FxDA/Dx/100*1e+6)
        xlabel('a [um]')
        ylabel('k [fN/nm]')
        drawnow()
    end
    
    % save
    save(fileDA,'aDA','FxDA','FyDA','FzDA')  
end
load(fileDA)

%% CALCULATION - GO
fileGO = 'GO.mat';
if ~exist(fileGO)
    
    r = Ray.beam2focused(bg,f);
    
    FxGO = zeros(size(aGO));
    FyGO = zeros(size(aGO));
    FzGO = zeros(size(aGO));
    for i = 1:1:length(aGO)
        
        bead = ParticleSpherical(C,aGO(i),nm,np);
        F = bead.force(r)
        FxGO(i) = sum(sum( F.Vx( isfinite(F.Vx) ) ));
        FyGO(i) = sum(sum( F.Vy( isfinite(F.Vy) ) ));
        FzGO(i) = sum(sum( F.Vz( isfinite(F.Vz) ) ));
        
        figure(2)
        cla
        hold on
        plot(aGO*1e+6,-FxGO/Dx/100*1e+6)
        xlabel('a [um]')
        ylabel('k [fN/nm]')
        drawnow()
    end

    % save
    save(fileGO,'aGO','FxGO','FyGO','FzGO')
end
load(fileGO)

%% CALCULATION - EMT
fileEMT = 'EMT.mat';
if ~exist(fileEMT)

    amin = 1;
    astep = 1;
    amax = length(a);

    %% Field calculation

    % Integration surface (Omega)
    N = 100;
    radius = 3e-6;

    QF = (2*pi*radius/N)/(lambda0/nm)

    [Theta,Phi] = meshgrid([pi/(2*N):pi/N:pi-pi/(2*N)],[0:pi/N:2*pi-pi/N]);
    r = radius*ones(size(Theta));
    [X,Y,Z] = Transform.Sph2Car(Theta,Phi,r);
    X = X-C.X;
    Y = Y-C.Y;
    Z = Z-C.Z;
    P = Point(X,Y,Z);

    % Focused beam angles
    phi = bg.phi;
    dphi = phi(2,1)-phi(1,1);

    theta = asin(bg.r/f);
    dr = bg.r(1,2)-bg.r(1,1);
    R = dr*size(bg.r,2);
    dtheta = ones(size(bg.r,1),1)*(asin([dr:dr:R]/f)-asin([0:dr:R-dr]/f));

    k = 2*pi*sqrt(er2*mr2)/bg.lambda0;

    % Field calculations
    for i = amin:astep:amax

        % Mie particles
        mie{i} = MieParticle(nm,np,a(i),k0);    

        % Fields
        E{i} = ComplexVector(X,Y,Z,zeros(size(X)),zeros(size(X)),zeros(size(X)));
        B{i} = ComplexVector(X,Y,Z,zeros(size(X)),zeros(size(X)),zeros(size(X)));

    end
    tic
    for n = 1:1:numel(theta)

        disp('')
        disp(['CALCULATING RAY NUMBER ' int2str(n) '/' int2str(numel(theta)) ' ... ' num2str(toc) 's'])
        disp('')

        % Unit vectors
        utheta = Point(cos(phi(n)).*cos(theta(n)), sin(phi(n)).*cos(theta(n)), -sin(theta(n)));
        utheta = utheta.tovector();
        uphi = Point(-sin(phi(n)), cos(phi(n)), zeros(size(phi(n))));
        uphi = uphi.tovector();

        % Constants
        shift = exp(1i * (k*sin(theta(n)).*cos(phi(n))*C.X+k*sin(theta(n)).*sin(phi(n))*C.Y+k*cos(theta(n))*C.Z) );
        factor = 1i*k*f*exp(-1i*k*f)/(2*pi) * shift * sin(theta(n)) * dphi*dtheta(n);

        % Radial component
        Pr = P.zrotation(-phi(n)).yrotation(-theta(n));
        [theta_r,phi_r,r_r] = Transform.Car2Sph(Pr.X,Pr.Y,Pr.Z);
        multi_r = Multipole(theta_r,phi_r,r_r,k);

        % Azimuthal component
        Pphi = P.zrotation(-phi(n)).yrotation(-theta(n)).zrotation(-pi/2);
        [theta_phi,phi_phi,r_phi] = Transform.Car2Sph(Pphi.X,Pphi.Y,Pphi.Z);
        multi_phi = Multipole(theta_phi,phi_phi,r_phi,k);

        parfor i = amin:1:amax

            % Radial component
            [Er_t,Br_t] = mie{i}.total(theta_r,phi_r,r_r,'Multipoles',multi_r);
            Er_t = Er_t.yrotation(theta(n)).zrotation(phi(n));
            Br_t = Br_t.yrotation(theta(n)).zrotation(phi(n));

            % Azimuthal component
            [Ephi_t,Bphi_t] = mie{i}.total(theta_phi,phi_phi,r_phi,'Multipoles',multi_phi);
            Ephi_t = Ephi_t.zrotation(pi/2).yrotation(theta(n)).zrotation(phi(n));
            Bphi_t = Bphi_t.zrotation(pi/2).yrotation(theta(n)).zrotation(phi(n));

            n1 = sqrt(er1);
            n2 = sqrt(er2);
            Et = (bg.Ephi(n)*Ephi_t + bg.Er(n)*Er_t)*(sqrt(n1/n2)*sqrt(cos(theta(n))));
            Bt = (bg.Ephi(n)*Bphi_t + bg.Er(n)*Br_t)*(sqrt(n1/n2)*sqrt(cos(theta(n))));

            E{i} = E{i} + ComplexVector(E{i}.X,E{i}.Y,E{i}.Z, ...
                Et.Vx * factor, ...
                Et.Vy * factor, ...
                Et.Vz * factor ...
                );
            B{i} = B{i} + ComplexVector(B{i}.X,B{i}.Y,B{i}.Z, ...
                Bt.Vx * factor, ...
                Bt.Vy * factor, ...
                Bt.Vz * factor ...
                );
        end
    end

    %% Force calculation
    [ux,uy,uz] = Transform.Sph2CarVector(Theta,Phi,zeros(size(r)),zeros(size(r)),ones(size(r)));
    ur = Vector(E{amax}.X,E{amax}.Y,E{amax}.Z,ux,uy,uz);

    dOmega = sin(Theta) * radius * (Theta(1,2)-Theta(1,1)) * radius * (Phi(2,1)-Phi(1,1));

    Fx = zeros(size(a));
    Fy = zeros(size(a));
    Fz = zeros(size(a));
    for i = amin:astep:amax

        ExEr = ComplexVector(E{i}.X,E{i}.Y,E{i}.Z, ...
            E{i}.Vx .* (conj(E{i}).*ur), ...
            E{i}.Vy .* (conj(E{i}).*ur), ...
            E{i}.Vz .* (conj(E{i}).*ur) ...
            );
        BxBr = ComplexVector(B{i}.X,B{i}.Y,B{i}.Z, ...
            B{i}.Vx .* (conj(B{i}).*ur), ...
            B{i}.Vy .* (conj(B{i}).*ur), ...
            B{i}.Vz .* (conj(B{i}).*ur) ...
            );
        df = .5*PhysConst.e0*mie{i}.nm^2 .* ( ...
            ( ExEr + (PhysConst.c0./mie{i}.nm)^2*BxBr ) .* dOmega ...
            - .5*( E{i}.norm().^2 + (PhysConst.c0./mie{i}.nm)^2*B{i}.norm().^2 ) .* ur .* dOmega ...
            ) ;

        F = Vector(0,0,0,real( sum(sum(df.Vx)) ),real( sum(sum(df.Vy)) ),real( sum(sum(df.Vz)) ) )
        Fx(i) = F.Vx;
        Fy(i) = F.Vy;
        Fz(i) = F.Vz;

    end

    %% Save
    save(fileEMT,'a','Fx','Fy','Fz','Nphi','Nr')
end
load(fileEMT)

%% Overall plot
figure(1)
cla    
hold on
plot(aDA(1:33)*1e+6,-FxDA(1:33)/Dx/100*1e+6,'k:')
plot(aGO*1e+6,-FxGO/Dx/100*1e+6,'k--')
plot(a(Fx~=0)*1e+6,-Fx(Fx~=0)/Dx/100*1e+6,'r.')
xlabel('a [um]')
ylabel('k [fN/nm]')
hold off
figure(2)
cla
hold on
plot(aDA(1:33)*1e+6,-FxDA(1:33)*1e+12,'k:')
plot(aGO*1e+6,-FxGO*1e+12,'k--')
plot(a(Fx~=0)*1e+6,-Fx(Fx~=0)*1e+12,'r.')
xlabel('a [um]')
ylabel('F [pN]')