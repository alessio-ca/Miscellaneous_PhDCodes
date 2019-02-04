%% INITIALIZATION
clear all;
close all;
clc;

%Domain
Lx = 2000e-9;
Ly = 2000e-9;
Lz = 2000e-9;

L = 32;

ep = 2.25;
a = .527e-6/2;
lambda0 = 1064e-9;
alpharc = InducedDipole.polarizability('corrected',a,ep,'lambda0',lambda0);

%% Setup in POM 0.26
% Beam waist = 10.7e-3 m
% f = 200 mm (tube length for infinity corrected) / 60 (obj magnification)
% Pmax = 60/100 (obj transmission) * 64/100 (AOD transmission) * 2 W

% Beam & objective

w0 = 10.7e-3;
Ex0 = 1;
Ey0 = 0;
nm = 1.33;
NA = 1.20;
f = 200e-3/60;
R = f*NA/nm;
f = nm*R/NA;

% Focused field
Nphi = 16;
Nr = 10;
bg = BeamGauss(Ex0,Ey0,w0,R,Nphi,Nr,'lambda0',1064e-9);
bg = bg.normalize(60*64*2/10000);
E_inc = EFieldFocus(bg,f,'lambda0',1064e-9,'er',1.33^2);

% Field
meshStep = .05e-6;
[x,y,z] = meshgrid(-Lx - 3*meshStep:meshStep:Lx + 3*meshStep,-Ly - 3*meshStep:meshStep:Ly + 3*meshStep,0);
ef = E_inc.E(Point(x,y,z),'er',1.33.^2);
I = .5*PhysConst.c0/nm*(nm^2*PhysConst.e0)*norm(ef).^2;

figure(1)
surf(1e6*x,1e6*y,I)
shading flat
axis auto
xlabel('x [um]')
ylabel('y [um]')
zlabel('E [V/m]')



%% 2 INTERACTING DIPOLES: ITERATIVE SOLVER

%Fixed positions
x1=0;
y1=0;
z1=0;
r1=Point(x1,y1,z1);

phase_r1 = focus_phase(bg,f,x1,y1,z1,'er',1.33^2);


% Variable position: Coarse mesh (with 3 mesh steps space)
r = Point(x,y,z);

Ex=zeros(size(r,1),size(r,2));
Ey=zeros(size(r,1),size(r,2));
Ez=zeros(size(r,1),size(r,2));
Bx=zeros(size(r,1),size(r,2));
By=zeros(size(r,1),size(r,2));
Bz=zeros(size(r,1),size(r,2));
Fx=zeros(size(r,1),size(r,2));
Fy=zeros(size(r,1),size(r,2));
Fz=zeros(size(r,1),size(r,2));

rx=r.X(:);
ry=r.Y(:);
rz=r.Z(:);

%%
%Enable progress bar for parallel pool
try
    parpool;
catch ME
    if ~strcmp(ME.identifier,'parallel:convenience:ConnectionOpen')
        rethrow(ME)
    end
end

targetWorkCount = size(r,1)*size(r,2);
barWidth= int32( 30 );
p =  TimedProgressBar( targetWorkCount, barWidth, ...
    'Computing, wait for ', ', completed ', 'Concluded in ' );

%for k = 1:targetWorkCount
for k = 1:20
    r2=Point(rx(k),ry(k),rz(k));
    phase_r2 = focus_phase(bg,f,r2.X,r2.Y,r2.Z,'er',1.33^2);

    rP = [r1,r2];
    idP = InducedDipole(alpharc,'lambda0',lambda0,'er',1.33^2,'rd',rP(1));
    
    for i=2:length(rP)
        idP = [idP,InducedDipole(alpharc,'lambda0',lambda0,'er',1.33^2,'rd',rP(i))];
    end
    
    numDip=2;
    
    if (norm(Point(rx(k),ry(k),rz(k)) - r1)<(2*a - 6*meshStep)) 
        Ex(k)=0;
        Ey(k)=0;
        Ez(k)=0;
        Bx(k)=0;
        By(k)=0;
        Bz(k)=0;
        Fx(k)=0;
        Fy(k)=0;
        Fz(k)=0;
    else
        tol = 1;
        cc = 0;
        Ei_n_base=[E_inc.E(rP(1),'er',1.33.^2),E_inc.E(rP(2),'er',1.33.^2)];
        Bi_n_base=[E_inc.B(rP(1),'er',1.33.^2),E_inc.B(rP(2),'er',1.33.^2)];
        
        Ei_n=Ei_n_base;
        Bi_n=Bi_n_base;
        chk = zeros(1,numDip);

        while tol > 0.001 
            Ei_o = Ei_n;
            cc = cc + 1;        
            for i=1:length(rP)
                Ei_n(i) = Ei_n_base(i);
                Bi_n(i) = Bi_n_base(i);
            end
            for i=1:length(rP)
                for j=1:length(rP)
                    if i~=j
                        Ei_n(i) = Ei_n(i)+idP(j).E(rP(i),Ei_o(j));
                        Bi_n(i) = Bi_n(i)+idP(j).B(rP(i),Ei_o(j));
                    end
                end
            end
            
            for i=1:length(rP)
                chk(i) = norm(Ei_n(i) - Ei_o(i))./norm(Ei_o(i));
            end
            
            tol = max([real(chk),imag(chk)]);
        end
        
        %Field calculation
        
        Ex(k)=E_inc.E(r2,'er',1.33.^2).Vx;
        Ey(k)=E_inc.E(r2,'er',1.33.^2).Vy;
        Ez(k)=E_inc.E(r2,'er',1.33.^2).Vz;
        Bx(k)=E_inc.B(r2,'er',1.33.^2).Vx;
        By(k)=E_inc.B(r2,'er',1.33.^2).Vy;
        Bz(k)=E_inc.B(r2,'er',1.33.^2).Vz;
        
        for i=1:length(rP)-1
                Ex(k)=Ex(k)+idP(i).E(r2,Ei_n(i)).Vx;
                Ey(k)=Ey(k)+idP(i).E(r2,Ei_n(i)).Vy;
                Ez(k)=Ez(k)+idP(i).E(r2,Ei_n(i)).Vz;
                Bx(k)=Bx(k)+idP(i).B(r2,Ei_n(i)).Vx;
                By(k)=By(k)+idP(i).B(r2,Ei_n(i)).Vy;
                Bz(k)=Bz(k)+idP(i).B(r2,Ei_n(i)).Vz;
        end

        %Force calculation, 5 pts interpol.
        [x_temp,y_temp,z_temp] = meshgrid(r2.X-2*meshStep:meshStep:r2.X+2*meshStep,r2.Y-2*meshStep:meshStep:r2.Y+2*meshStep,0);
        r_temp = Point(x_temp,y_temp,z_temp);
        
        
        E_temp = E_inc.E(r_temp,'er',1.33.^2);
        B_temp = E_inc.B(r_temp,'er',1.33.^2);
        
        for i=1:length(rP)-1
                E_temp=E_temp + idP(i).E(r_temp,Ei_n(i));
                B_temp=B_temp + idP(i).B(r_temp,Ei_n(i));
        end
        
        [F_temp,~,~] = idP(2).force_general(r2,E_temp,B_temp,2);
        
        Fx(k)=F_temp.Vx;
        Fy(k)=F_temp.Vy;
        Fz(k)=F_temp.Vz;
        
        if (norm(Point(rx(k),ry(k),rz(k)) - r1)<(2*a-3*meshStep))
            Ex(k)=0;
            Ey(k)=0;
            Ez(k)=0;
            Bx(k)=0;
            By(k)=0;
            Bz(k)=0;
            Fx(k)=0;
            Fy(k)=0;
            Fz(k)=0;
        end
        
    end
    p.progress;

end
p.stop;

Et_PERF = ComplexVector(r.X,r.Y,r.Z,Ex,Ey,Ez);
Bt_PERF = ComplexVector(r.X,r.Y,r.Z,Bx,By,Bz);
F_PERF = ComplexVector(r.X,r.Y,r.Z,Fx,Fy,Fz);

% %% Plotting
% figure(1)
% subplot(1,2,1)
% surf(r.X,r.Y,norm(Et),'edgealpha',0)
% subplot(1,2,2)
% surf(norm(F),'edgealpha',0)
% 
% %% Optional: export of the calculated fields
% save('Force_beam_2.mat','F')
% 
% %% Optional: load the calculated fields
% load('Force_beam_2.mat')
% 
% %% Spline interpolant
% 
% F.Vx((F.X-x1).^2 + (F.Y-y1).^2 < (2*a-3*meshStep)^2)=NaN;
% F.Vy((F.X-x1).^2 + (F.Y-y1).^2 < (2*a-3*meshStep)^2)=NaN;
% 
% idxgood=~(isnan(F.Vx) | isnan(F.Vy));
% FxG = griddata(x(idxgood),y(idxgood),F.Vx(idxgood),x,y,'natural');
% FyG = griddata(x(idxgood),y(idxgood),F.Vy(idxgood),x,y,'natural');
% 
% spl_x=csapi({-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep:2e-6+3*meshStep},FxG');
% spl_y=csapi({-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep:2e-6+3*meshStep},FyG');
% 
% fun_x = @(X) fnval(spl_x,X);
% fun_y = @(X) fnval(spl_y,X);
% 
% [xx,yy,zz] = meshgrid(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,0);
% fineFx = fun_x({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep})';
% fineFy = fun_y({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep})';
% divF = divergence(xx,yy,fineFx,fineFy);
% spl_divF=csapi({-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep},divF');
% fun_divF = @(X) fnval(spl_divF,X);
% %% Optional: export of the calculated fields
% save('ForceGx_2_beam.mat','FxG')
% save('ForceGy_2_beam.mat','FyG')
% save('divF_2_beam.mat','divF')
% %% Testing fineFx & fineFy on y=0 and x=0 
% figure(2)
% subplot(2,2,1)
% plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFx(yy==0))
% hold on
% plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vx(F.Y==0),'o')
% hold off
% 
% subplot(2,2,2)
% plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFx(xx==0))
% hold on
% plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vx(F.X==0),'o')
% hold off
% 
% subplot(2,2,3)
% plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFy(yy==0))
% hold on
% plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vy(F.Y==0),'o')
% hold off
% 
% subplot(2,2,4)
% plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,fineFy(xx==0))
% hold on
% plot(-2e-6-3*meshStep:meshStep:2e-6+3*meshStep,F.Vy(F.X==0),'o')
% hold off
% 
% %% Testing divF on y=0 and x=-0.1979 (intersection point)
% divF((xx-x1).^2 + (yy-y1).^2 < (2*a)^2)=NaN; %Cut at particle diameter with an internal grid point
% edge=-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep;
% subplot(1,2,1)
% plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,divF(yy==0))
% hold on
% plot(edge,1e-4+fun_divF([edge;zeros(1,length(edge))]))
% hold off
% subplot(1,2,2)
% plot(-2e-6-3*meshStep:meshStep/10:2e-6+3*meshStep,divF(xx==0))
% hold on
% plot(edge,1e-4+fun_divF([zeros(1,length(edge));edge]))
% hold off
% 
% %%
% zmin = min(min(min(F.Vy)));
% zmax = max(max(max(F.Vy)));
% zinc = (zmax - zmin) / 20;
% zlevs = zmin:zinc:zmax;
% figure(3)
% g=subplot(1,2,1);
% %surf(x*1e+9,y*1e+9,norm(F),'edgealpha',0)
% contourf(x*1e+6,y*1e+6,F.Vy,zlevs);
% hold on
% h=streamslice(F.X(1:20:end,1:20:end)*1e+6,F.Y(1:20:end,1:20:end)*1e+6,F.Vx(1:20:end,1:20:end),F.Vy(1:20:end,1:20:end),4,'method','cubic');
% for i=1:length(h)
%     h(i).LineWidth=1;
% end
% set(h,'Color',[0.6 0.2 0.])
% colormap(g,'parula');
% hold off
% title('F_y and streamlines, xy plane')
% axis equal
% view(2)
% xlabel('x [\mum]')
% ylabel('y [\mum]')
% 
% magnitude=norm(F);
% subplot(1,2,2)
% surf(x*1e+6,y*1e+6,magnitude,'edgealpha',0)
% colormap jet
% title('Magnitude of Force, xy plane')
% axis equal
% view(2)
% xlabel('x [\mum]')
% ylabel('y [\mum]')
% 
% % magnitude=norm(E1);
% % figure(4)
% % surf(x*1e+9,y*1e+9,magnitude,'edgealpha',0)
% % title('xy plane')
% % axis equal
% % view(2)
% % xlabel('x [nm]')
% % ylabel('y [nm]')
% 
