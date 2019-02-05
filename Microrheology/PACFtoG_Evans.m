function [omega,G,G_err,G_vec]=PACFtoG_Evans(t,data,fps,varargin)
% PACFtoG_Evans PACF to G' and G'' conversion according to Evans (2010)
%
% [omega,G,G_err,G_vec]=PACFtoG_Evans(t,data,fps) calculates the
%  complex rheological modulus G from the Fourier transform of a function 
%  f(t) using the scheme by Evans et al. (2010), which assumes that f(t) 
%  vanishes for t<0 and is sampled at a finite set of data points (T_k,DATA_k) 
%  which need not be equally spaced. f(t) has been obtained by a signal
%  sampled at frequency FPS.
%  DATA can be a matrix containing several signals (one per row, all the same length).
%  Gerr is the error calculated as standard deviation of the G value.
%  Gvec is the set of all G calculated

% [omega,G,G_err,G_vec]=PACFtoG_Evans(t,data,fps,'PropertyName',PropertyValue)
%  permits to set the value of PropertyName to PropertyValue.
%  Admissible Properties are:
%       TauMax      -   Maximum time delay to be considered (default = +Inf)
%       f0          -   Value of f(t) at t=0 (default = 1)
%       fpinf       -   Value of f(t) at t->Inf (default = 0)
%       Gfactor     -   Conversion factor (default = 1)


% CREATED: Alessio Caciagli, University of Cambridge, January 2017

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

if t(1)==0
    %The first data point must NOT be at t=0
    t=t(2:end);
    data=data(2:end,:);
end

%Define omega vector 
minfreq=1/t(end);
maxfreq=fps/2;
omega=2*pi*logspace(log10(minfreq),log10(maxfreq),100)';

%Evans_fft is composed of three single terms and a summatory (see Evans et
%al., (2010))

%Calculate the single terms
evans_fft = repmat(1i*omega,1,size(data,2))*f0;
evans_fft = evans_fft + (1-exp(-1i*omega*t(1)))*(data(1,:)-f0)/t(1);
evans_fft = evans_fft + fpinf*repmat(exp(-1i*omega*t(end)),1,size(data,2));

%Calculate the summatory
for i=1:length(omega)
    evans_fft(i,:)=evans_fft(i,:)-sum(diff(exp(-1i*omega(i)*t)).*diff(data)./diff(t));
end

%Split into real and imaginary parts
fft_real = real(evans_fft)./repmat(omega.^2,1,size(data,2));
fft_imag = -(imag(evans_fft)./repmat(-omega.^2,1,size(data,2)) + 1./repmat(omega,1,size(data,2)));
fft_tot = fft_imag.^2 + fft_real.^2;

%Calculate G moduli
Gp = -fft_imag./(repmat(omega,1,size(data,2)).*fft_tot); %Uncorrected for trap!
Gpp = -fft_real./(repmat(omega,1,size(data,2)).*fft_tot);

%Finalize
G = complex(mean(Gp,2),mean(Gpp,2));
G_err = complex(std(Gp,0,2),std(Gpp,0,2));
G_vec = complex(Gp,Gpp);



end