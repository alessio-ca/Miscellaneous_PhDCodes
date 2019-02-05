function [omega,G,G_err,G_vec]=MSDtoG_Mason(t,data,varargin)
% MSDtoG_Mason MSD to G' and G'' conversion according to Mason & Dasgupta
% (1996 and 2002)
%
% [omega,G,G_err,G_vec]=MSDtoG_Mason(t,data,fps) calculates the
%  complex rheological modulus G from the Fourier transform of the MSD 
%  following the scheme of Mason & Dasgupta et al. (1997 & 2002), which 
%  assumes that the MSD follows a local power law. The procedure estimates
%  the power law behavior from the fist and second logarithmic time 
%  derivatives of the MSD. 
%  The scheme interpolates the MSD onto a spline to yield better results 
%  in the derivative estimation. 
%
%  DATA can be a matrix containing several signals (one per row, all the 
%  same length).
%  Gerr is the error calculated as standard deviation of the G value.
%  Gvec is the set of all G calculated
%
% [omega,G,G_err,G_vec]=MSDtoG_Mason(t,data,fps,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       T           -   Temperature (default = 298 K)
%       R           -   Bead radius (default = 500 nm)
%       Dim         -   Dimensionality of the system (default = 3)
%       CG          -   Coarse-graining factor (default = 1.45)
%       TauMax      -   Maximum time delay to be considered (default = +Inf)
%       Gsampling   -   Sampling style for G points (default = "log")
%       Gpoints     -   Number of Gpoints (default = 100)
%           Allowed styles: logarithmic (log), linear (lin)
%       width       -   Width of Gaussian window for polyfitw (default =
%                       0.7)
%       cutoff      -   Cutoff treshold for point discarding (default =
%                       0.03)

% CREATED: Alessio Caciagli, University of Cambridge, January 2018

T = 298;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'T')
        T = varargin{n+1};
    end
end
R = 115e-9;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'R')
        R = varargin{i+1};
    end
end
Dim = 3;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i},'Dim')
        Dim = varargin{i+1};
    end
end
taumax = +Inf;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'taumax')
        taumax = varargin{n+1};
    end
end
CG = 1.45;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'CG')
        CG = varargin{n+1};
    end
end
Gsampling = 'log';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Gsampling')
        Gsampling = varargin{n+1};
    end
end
Gpoints = 100;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Gpoints')
        Gpoints = varargin{n+1};
    end
end
width = 0.7;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'width')
        width = varargin{n+1};
    end
end
cutoff = 0.03;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'cutoff')
        cutoff = varargin{n+1};
    end
end

if t(1)==0
    %The first data point must NOT be at t=0
    t=t(2:end);
    data=data(2:end,:);
end

t = t(t<taumax);
data=data(1:length(t),:);

%Downsampling & spline interpolant
interval=1;
res=1;

while (interval <length(t))
    res=res+1;
    interval=ceil(CG^res);
end
res=res-1;
interval=ceil(CG.^(0:res));
interval=unique(interval);
plotinterval=unique(floor(logspace(0,log10(length(t)),1000)));
t_downsample=t(interval);
data_downsample=data(interval,:);
spline_fit=csape(t_downsample',data_downsample');

%Plotting section
fig1=figure(1);
clf
loglog(t(plotinterval),data(plotinterval,1),'o')
hold on
loglog(t(plotinterval),data(plotinterval(1),1)*data(plotinterval(end),1)./data(plotinterval,1),'o')
loglog(t_downsample,data_downsample(:,1),'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
spline_downsample = fnval(t_downsample,spline_fit);

loglog(t_downsample',spline_downsample(1,:),'--k');
%set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('Exp','Exp (inv)','Log-sampling','Location','northeast')
xlabel('Time (s)')
ylabel('MSD (t)')
hold off
drawnow

disp(' ')
disp('Check your MSD (t). If you are satisfied...')
disp('Press the RETURN key to GO ON (make sure the Figure is active!).') 
disp('Press any other key to STOP.')

%pause; % wait for a keypress
k = waitforbuttonpress;
currkey=get(gcf,'CurrentKey');
if strcmp(currkey, 'return')
    currkey=1;
else
    return;
end

disp(' ')
disp('Calculating...')

%Define omega vector 
switch Gsampling
    case 'log'
        t_G = logspace(log10(t_downsample(1)),log10(t_downsample(end)),Gpoints)';
        
    case 'lin'
        t_G = linspace(t_downsample(1),t_downsample(end),Gpoints)';
        
    otherwise
        error('Error: Unrecognized G sampling style.');
end
omega=flipud(1./t_G);
d0=fnval(t_G,spline_fit)';
if size(d0,1)==1
    d0=d0';
end


% use 2nd order local formula for G(s)-- good to 1% of G(s)
log_spline_fit=csape(log(t_G)',log(d0)');
d_log_spline_fit = fnder(log_spline_fit,1);
dd_log_spline_fit = fnder(log_spline_fit,2);

s_d0=fnval(t_G,spline_fit)';
s_d1=fnval(log(t_G),d_log_spline_fit)';
s_d2=fnval(log(t_G),dd_log_spline_fit)';
if size(s_d0,1)==1
    s_d0=s_d0';
    s_d1=s_d1';
    s_d2=s_d2';
end

[d0,d1,d2] = Mason_logderive(t_G,d0,width);

kB = 1.38e-23;
C = Dim*kB*T/(3*pi*R);  %Pre-factor
prefac = (pi/2)-1; 
s_Gs = flipud(C./((s_d0.*gamma(1+s_d1)).*(1 +(s_d2/2))));
Gs = flipud(C./((d0.*gamma(1+d1)).*(1 +(d2/2))));
% use 2nd order local formula for G'(w),G"(w)-- good to 2% of G(w)
log_spline_G=csape(log(omega)',log(s_Gs)');
d_log_spline_G = fnder(log_spline_G,1);
dd_log_spline_G = fnder(log_spline_G,2);

s_d1G=fnval(log(omega),d_log_spline_G)';
s_d2G=fnval(log(omega),dd_log_spline_G)';
if size(s_d1G,1)==1
    s_d1G=s_d1G';
    s_d2G=s_d2G';
end

[d0G,d1G,d2G] = Mason_logderive(omega,Gs,width);

s_Gp  = s_Gs.*(1./(1+s_d2G)).* ( cos( (pi/2)*s_d1G)   -     prefac*s_d1G.*s_d2G);     
s_Gpp = s_Gs.*(1./(1+s_d2G)).* ( sin( (pi/2)*s_d1G)   -     prefac*(1-s_d1G).*s_d2G);

Gp  = Gs.*(1./(1+d2G)).* ( cos( (pi/2)*d1G)   -     prefac*d1G.*d2G);     
Gpp = Gs.*(1./(1+d2G)).* ( sin( (pi/2)*d1G)   -     prefac*(1-d1G).*d2G);    

% clip off the suspicious (i.e. unreliable) data
for i = 1:size(data,2)
    w = find(Gp(:,i) < Gs(:,i)*cutoff);
    Gp(w,i)=0;
    w = find(Gpp(:,i) < Gs(:,i)*cutoff); 
    Gpp(w,i)=0;
end

if any((max(abs(d2)) - 0.15)>0) || any((max(abs(d2G)) - 0.15)>0)  %original value: 0.15
    warning ('Warning, high curvature in data, moduli may be unreliable!')
end



G = complex(mean(Gp,2),mean(Gpp,2));
G_err = complex(std(Gp,0,2),std(Gpp,0,2));
G_vec = complex(Gp,Gpp);

disp('Done!')

end