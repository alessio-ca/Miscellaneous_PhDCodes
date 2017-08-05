function [omega,G,G_err,G_vec]=PACFtoG_Evans_oversampling(t,data,fps,varargin)
% PACFtoG_Evans PACF to G' and G'' conversion according to Evans (2010)
%
% [omega,G,G_err,G_vec]=PACFtoG_Evans(t,data,fps) calculates the
%  complex rheological modulus G from the Fourier transform of a function 
%  f(t) using the scheme by Evans et al. (2010), which assumes that f(t) 
%  vanishes for t<0 and is sampled at a finite set of data points (T_k,DATA_k) 
%  which need not be equally spaced. f(t) has been obtained by a signal
%  sampled at frequency FPS. 
%  The scheme samples the PACF logarithmically and oversamples it to yield
%  best results.
%
%  DATA can be a matrix containing several signals (one per row, all the same length).
%  Gerr is the error calculated as standard deviation of the G value.
%  Gvec is the set of all G calculated

% [omega,G,G_err,G_vec]=PACFtoG_Evans(t,data,fps,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       TauMax      -   Maximum time delay to be considered (default = +Inf)
%       f0          -   Value of f(t) at t=0 (default = 1)
%       fpinf       -   Value of f(t) at t->Inf (default = 0)
%       CG          -   Coarse-graining factor (default = 1.2)
%       Beta        -   Oversampling factor (default = 100)
%       Gfactor     -   Multiplication factor for G* (default = 1)

% CREATED: Alessio Caciagli, University of Cambridge, April 2017

taumax = +Inf;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'taumax')
        taumax = varargin{n+1};
    end
end
f0 = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'f0')
        f0 = varargin{n+1};
    end
end
fpinf = 0;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'fpinf')
        fpinf = varargin{n+1};
    end
end
CG = 1.2;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'CG')
        CG = varargin{n+1};
    end
end
beta = 100;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'beta')
        beta = varargin{n+1};
    end
end
Gfactor = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Gfactor')
        Gfactor = varargin{n+1};
    end
end

if t(1)==0
    %The first data point must NOT be at t=0
    t=t(2:end);
    data=data(2:end,:);
end

t = t(t<taumax);
data=data(1:length(t),:);

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
f=figure(1);
semilogx(t(plotinterval),data(plotinterval,1),'o')
hold on
semilogx(t_downsample,data_downsample(:,1),'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
spline_downsample = fnval(t_downsample,data_spline);
h=semilogx(t_downsample',spline_downsample(1,:),'--k');
%set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend('ACF','Log-sampling','Location','northeast')
xlabel('Time (s)')
ylabel('A(t)')
hold off
drawnow

disp(' ')
disp('Check your ACF. If you are satisfied...')
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
omega=2*pi*resfreq*logspace(log10(minfreq),log10(maxfreq),100)';
t_oversample=(1/(beta*fps):1/(beta*fps):t_downsample(end))';
data_oversample=fnval(t_oversample,data_spline)';

if size(data_oversample,1)==1
    data_oversample=data_oversample';
end


%Evans_fft is composed of three single terms and a summatory (see Evans et
%al., (2010))

%Calculate the single terms
evans_fft = repmat(1i*omega,1,size(data_oversample,2))*f0;
evans_fft = evans_fft + (1-exp(-1i*omega*t_oversample(1)))*(data_oversample(1,:)-f0)/t_oversample(1);
evans_fft = evans_fft + fpinf*repmat(exp(-1i*omega*t_oversample(end)),1,size(data_oversample,2));

%Calculate the summatory
for i=1:length(omega)
    evans_fft(i,:)=evans_fft(i,:)-sum(diff(exp(-1i*omega(i)*t_oversample)).*diff(data_oversample)./diff(t_oversample));
end

%Split into real and imaginary parts
fft_real = -real(evans_fft)./repmat(-omega.^2,1,size(data,2));
fft_imag = -(imag(evans_fft)./repmat(-omega.^2,1,size(data,2)) + 1./repmat(omega,1,size(data,2)));
fft_tot = fft_imag.^2 + fft_real.^2;

%Calculate G moduli (N.B! G is unitless unless you specify a Gfactor with suitable units)
Gp = -fft_imag./(repmat(omega,1,size(data,2)).*fft_tot); %Uncorrected for trap!
Gpp = -fft_real./(repmat(omega,1,size(data,2)).*fft_tot);

%Finalize
G = Gfactor*complex(mean(Gp,2),mean(Gpp,2));
G_err = Gfactor*complex(std(Gp,0,2),std(Gpp,0,2));
G_vec = Gfactor*complex(Gp,Gpp);

disp('Done!')

end