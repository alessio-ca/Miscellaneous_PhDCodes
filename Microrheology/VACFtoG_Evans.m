function [omega,G,G_err,G_vec]=VACFtoG_Evans(t,data,fps,varargin)
% PACFtoG_Evans VACF to G' and G'' conversion according to Evans (2010)
%
% [omega,G,G_err,G_vec]=VACFtoG_Evans(t,data,fps) calculates the
%  complex rheological modulus G from the Fourier transform of a function 
%  f(t) using the scheme by Evans et al. (2010), which assumes that f(t) 
%  vanishes for t<0 and is sampled at a finite set of data points (T_k,DATA_k) 
%  which need not be equally spaced. f(t) has been obtained by a signal
%  sampled at frequency FPS. The conversion follows the approach by
%  Yaganashima et al. (2011).
%  DATA can be a matrix containing several signals (one per row, all the same length).
%  Gerr is the error calculated as standard deviation of the G value.
%  Gvec is the set of all G calculated

% [omega,G,G_err,G_vec]=VACFtoG_Evans(t,data,fps,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       TauMax      -   Maximum time delay to be considered (default = +Inf)
%       Temp        -   Temperature (default = 300)
%       Diam        -   Particle diameter (default = 1e-6 m)


% CREATED: Alessio Caciagli, University of Cambridge, January 2017

taumax = +Inf;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'taumax')
        taumax = varargin{n+1};
    end
end
if isfinite(taumax)
    maxlags = ceil(taumax*fps);
else
    maxlags = size(data,1);
end
t=t(1:maxlags);
data=data(1:maxlags,:);

%Define constants
kB = 1.3806503 * 1E-23;
Temp = 300;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'temp')
        Temp = varargin{n+1};
    end
end
R=0.5e-6;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'diam')
        R = 0.5*varargin{n+1};
    end
end

const = kB*Temp/(6*pi*R);


%Define omega vector 
minfreq=1/t(end);
maxfreq=fps/2;
omega=2*pi*logspace(log10(minfreq),log10(maxfreq),100)';

%Evans_fft is composed of three single terms and a summatory (see Evans et
%al., (2010))

%Account for delta like first point in correlation (see Yaganashima et al.,
%(2011))
C0Corr = data(1,:) - data(2,:);
data(1,:)=data(2,:);

%Calculate the single terms
evans_fft = repmat(0.5 * C0Corr * t(2),length(omega),1);

%Calculate the summatory (with correction)
for i=1:length(omega)
    evans_fft(i,:)=evans_fft(i,:) - sum(exp(-1i*omega(i)*t(2:end-1)).*(diff(data(2:end,:))./diff(t(2:end)) - diff(data(1:end-1,:))./diff(t(1:end-1))))./repmat(omega(i).^2,1,size(data,2));
    evans_fft(i,:)=evans_fft(i,:) - 1i*data(1,:)./repmat(omega(i),1,size(data,2));
end

%Obtain D(s) dividing for s=iw
D=-1i*evans_fft./repmat(omega,1,size(data,2));

%Calculate G moduli
Gp = real(const./D); %Uncorrected for trap!
Gpp = imag(const./D); 

%Finalize
G = complex(mean(Gp,2),mean(Gpp,2));
G_err = complex(std(Gp,0,2),std(Gpp,0,2));
G_vec = complex(Gp,Gpp);


% % % % %calculate via FFT (just for testing)
% % % % sampling=[1:5,floor(logspace(log10(6),log10(length(t)/2),100-5))];
% % % % freq=1/t(end):fps/length(t):fps/2;
% % % % %fft_fft=fft(data)/(2*pi*length(t));
% % % % fft_fft=fft(data)/fps;
% % % % fft_fft=fft_fft(2:length(fft_fft)/2+1);
% % % % 
% % % % omega_fft=2*pi*freq(sampling)';
% % % % fft_fft=fft_fft(sampling);
% % % % D_fft=-1i*fft_fft./repmat(omega_fft,1,size(data,2));
% % % % Gp_fft = real(const./D_fft); %Uncorrected for trap!
% % % % Gpp_fft = imag(const./D_fft); 
% % % % 
% % % % G_fft = complex(mean(Gp_fft,2),mean(Gpp_fft,2));



end