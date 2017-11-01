function [omega,G,G_err,G_vec,DeltaA]=MSDtoG_Evans_oversampling(t,data,fps,varargin)
% MSDtoG_Evans PACF to G' and G'' conversion according to Evans (2010)
%
% [omega,G,G_err,G_vec]=MSDtoG_Evans(t,data,fps) calculates the
%  complex rheological modulus G from the Fourier transform of a function 
%  f(t) using the scheme by Evans et al. (2010), which assumes that f(t) 
%  vanishes for t<0 and is sampled at a finite set of data points (T_k,DATA_k) 
%  which need not be equally spaced. f(t) has been obtained by a signal
%  sampled at frequency FPS. 
%  The scheme samples the MSD logarithmically and oversamples it to yield
%  best results.
%
%  DATA can be a matrix containing several signals (one per row, all the same length).
%  Gerr is the error calculated as standard deviation of the G value.
%  Gvec is the set of all G calculated
%
% [omega,G,G_err,G_vec]=MSDtoG_Evans(t,data,fps,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       TauMax      -   Maximum time delay to be considered (default = +Inf)
%       Nf0         -   Number of points to extrapolate f(t) at t=0 (default = 10)
%       Nfinf       -   Number of points to extrapolate f(t) at t=+Inf (default = 100)
%       CG          -   Coarse-graining factor (default = 1.45)
%       Beta        -   Oversampling factor (default = 100)
%       Jfactor     -   Conversion factor from MSD to J (default = 1)
%       Gsampling   -   Sampling style for G points (default = "log")
%           Allowed styles: logarithmic (log), linear (lin)
%       Eps         -   Filter epsilon for residuals (default = 0)
%       BM          -   Black Magic, secret weapon of Alessio to make
%                       pretty plots if they oscillate (default = 0 because
%                       we do not cheat... usually).

% CREATED: Alessio Caciagli, University of Cambridge, August 2017

taumax = +Inf;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'taumax')
        taumax = varargin{n+1};
    end
end
Nf0 = 10;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Nf0')
        Nf0 = varargin{n+1};
    end
end
Nfinf = 100;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Nfinf')
        Nfinf = varargin{n+1};
    end
end
CG = 1.45;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'CG')
        CG = varargin{n+1};
    end
end
beta = 100;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Beta')
        beta = varargin{n+1};
    end
end
Jfactor = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Jfactor')
        Jfactor = varargin{n+1};
    end
end
Gsampling = 'log';
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Gsampling')
        Gsampling = varargin{n+1};
    end
end
Eps = 0;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Eps')
        Eps = varargin{n+1};
    end
end
BM = 0;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'BM')
        BM = varargin{n+1};
    end
end

if t(1)==0
    %The first data point must NOT be at t=0
    t=t(2:end);
    data=data(2:end,:);
end


if Nf0==0
    %Specify J0
    %TO DO
end
if Nfinf==0
    %Specify Jinf
    %TO DO
end

t = t(t<taumax);
data=data(1:length(t),:);

%Convert MSD (<Dx^2(t)>) to creep (J(t))
data=Jfactor*data;

%Fit at t-> 0
%N.B! To make general for a matrix
t_fit_0 = t(1:Nf0);
data_fit_0 = data(1:Nf0,:);
%First do a power law
A_0_lin = [ones(length(t_fit_0),1),log(t_fit_0)];
bf_0_lin = A_0_lin\log(data_fit_0);
%Then do a nonlin fitting
A_0 = @(x,xdata) x(1) + x(2).*xdata.^x(3);
x0 = [data_fit_0(1),exp(bf_0_lin(1)),bf_0_lin(2)];
[bf_0,resnorm,~,exitflag,output] = lsqcurvefit(A_0,x0,t_fit_0,data_fit_0)



%Fit at t-> Inf
%N.B! To make general for a matrix
t_fit_inf = t(end-Nfinf:end);
data_fit_inf = data(end-Nfinf:end,:);
%First do a power law
A_inf_lin = [ones(length(t_fit_inf),1),log(t_fit_inf)];
bf_inf_lin = A_inf_lin\log(data_fit_inf);
%Then improve by nonlin fitting
A_inf = @(x,xdata) x(1) + x(2).*xdata.^x(3);
x0 = [0,exp(bf_inf_lin(1)),bf_inf_lin(2)];
[bf_inf,resnorm,~,exitflag,output] = lsqcurvefit(A_inf,x0,t_fit_inf,data_fit_inf)


%Show fit and ask if happy
%Plotting section
fig1=figure(1);
subplot(1,2,1)
plot(t_fit_0,data_fit_0,'o')
hold on
plot(t_fit_0,A_0(bf_0,t_fit_0),'--k');
legend('Exp','Fit','Location','northeast')
xlabel('Time (s)')
ylabel('J(t)')
hold off

subplot(1,2,2)
plot(t_fit_inf,data_fit_inf,'o')
hold on
%loglog(t_fit_inf,exp(A_inf*bf_inf),'--k');
plot(t_fit_inf,A_inf(bf_inf,t_fit_inf),'--k');
legend('Exp','Fit','Location','northeast')
xlabel('Time (s)')
ylabel('J(t)')
hold off

drawnow
disp(' ')
disp('Check your fits. If you are satisfied...')
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
disp(['J0 = ',num2str(bf_0(1))])
disp(['fpinf = ',num2str(bf_inf(3))])
if bf_0(1)<0 || abs(bf_inf(3)-1)>0.1
    disp('Warning: fit coefficients are ill defined.')
    if bf_0(1) < 0
        disp('Forcing J0 = 0.')
        bf_0 = [0,exp(bf_0_lin(1)),bf_0_lin(2)];
    end
    if abs(bf_inf(3)-1)>0.1
        disp('System does not display a long time eta along the measured times. Proceeding anyway...')
    end
end

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
data_spline=csape(t_downsample',data_downsample');

%Plotting section
fig1=figure(1);
clf
loglog(t(plotinterval),data(plotinterval,1),'o')
hold on
loglog(t(plotinterval),data(plotinterval(1),1)*data(plotinterval(end),1)./data(plotinterval,1),'o')
loglog(t_downsample,data_downsample(:,1),'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
spline_downsample = fnval(t_downsample,data_spline)';
loglog(t_downsample',spline_downsample(1,:)','--k');
%set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('Exp','Exp (inv)','Log-sampling','Location','northeast')
xlabel('Time (s)')
ylabel('J(t)')
hold off
drawnow

disp(' ')
disp('Check your J(t). If you are satisfied...')
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
resfreq=1/t(end);
minfreq=1/(resfreq*t_downsample(end));
maxfreq=fps/(2*resfreq); 
switch Gsampling
    case 'log'
        omega=2*pi*resfreq*logspace(log10(minfreq),log10(maxfreq),100)';
        
    case 'lin'
        omega=2*pi*resfreq*linspace(minfreq,maxfreq,100)';
        
    otherwise
            error('Error: Unrecognized G sampling style.');
end
    
t_oversample=(1/(beta*fps):1/(beta*fps):t_downsample(end))';
data_oversample=fnval(t_oversample,data_spline)';

if size(data_oversample,1)==1
    data_oversample=data_oversample';
end


%Evans_fft is composed of three single terms and a summatory (see Evans et
%al., (2010))

f0 = bf_0(1);
%fpinf = bf_inf(2);
fpinf = bf_inf(2)*bf_inf(3)*t_oversample(end).^(bf_inf(3)-1);


%Calculate the single terms
evans_fft = repmat(1i*omega,1,size(data_oversample,2))*f0;

A_Evans = diff(data_oversample)./diff(t_oversample);
if BM == 0
    A_Evans = [0;(data_oversample(1) - f0)/t_oversample(1);A_Evans;fpinf];
    t_oversample = [0;t_oversample];
else
    A_Evans = [0;A_Evans];
    t_oversample=t_oversample(1:end-1);
end
diffA_Evans = diff(A_Evans);
diffA_Evans(abs(diffA_Evans)<Eps)=0;

%Output residuals
datainterval=unique(floor(logspace(0,log10(length(t_oversample)),1000)));
DeltaA=abs(diffA_Evans(datainterval,1));


%Calculate the summatory
for i=1:length(omega)
    evans_fft(i,:)=evans_fft(i,:)+sum(diffA_Evans.*exp(-1i*omega(i)*t_oversample));
end

%Invert the result and multiply
G = repmat(1i*omega,1,size(data_oversample,2))./evans_fft;
G=G(omega<omega(end)/10);
omega=omega(omega<omega(end)/10); %Exclude the upper 10% due to artifacts
Gp = real(G);
Gpp = imag(G);
G_err = complex(std(Gp,0,2),std(Gpp,0,2));
G_vec = complex(Gp,Gpp);

disp('Done!')

end