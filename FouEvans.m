function [omega,Gp,Gpp]=FouEvans(tau,data,fps,f0,fpinf,cutoff)

%Calculate the Fourier transform of a function f(t) using the scheme by
%Evans et al. (2010), which assumes that f(t) vanishes for t<0 and is
%sampled at a finite set of data points (tk,fk) which need not be equally
%spaced. It assumes 2D data (x(t) and y(t) series).

%Parameters:
%tau is the time vector
%data is the f(t) vector (contains ACFx(t) and ACFy(t))
%f0 is the value f(0) (need to be the same for x(t) and y(t))
%fpinf is the gradient of f extrapolated to infinite time (need to be the
%same for x(t) and y(t)).

%Optional:
%cutoff is a max value on time (useful to eliminate noise in the tails of
%the sample f(t))

if nargin == 6
    %Apply cutoff
    data=data(tau<cutoff,:);
    tau=tau(tau<cutoff);
end

if tau(1)==0
    %The first data point must NOT be at t=0
    tau=tau(2:end);
    data=data(2:end,:);
end

%Define omega vector 
minfreq=1/tau(end);
maxfreq=fps/2;
omega=2*pi*logspace(log10(minfreq),log10(maxfreq),100)';

%Evans_fft is composed of three single terms and a summatory (see Evans et
%al., (2010))

%Calculate the single terms
evans_fft = repmat(1i*omega,1,2)*f0;
evans_fft = evans_fft + [(1-exp(-1i*omega*tau(1)))*(data(1,1)-f0)/tau(1),(1-exp(-1i*omega*tau(1)))*(data(1,2)-f0)/tau(1)];
evans_fft = evans_fft + fpinf*repmat(exp(-1i*omega*tau(end)),1,2);

%Calculate the summatory
for i=1:length(omega)
    evans_fft(i,:)=evans_fft(i,:)-[sum(diff(exp(-1i*omega(i)*tau)).*diff(data(:,1))./diff(tau)),sum(diff(exp(-1i*omega(i)*tau)).*diff(data(:,2))./diff(tau))];
end

%Split into real and imaginary parts
evans_fft = [omega,real(evans_fft),imag(evans_fft)];
evans_fft(:,2:3) = evans_fft(:,2:3)./repmat(evans_fft(:,1).^2,1,2);
evans_fft(:,4:5) = -(evans_fft(:,4:5)./repmat(-evans_fft(:,1).^2,1,2) + 1./repmat(evans_fft(:,1),1,2));
evans_fft(:,6:7) = evans_fft(:,2:3).^2 + evans_fft(:,4:5).^2;

%Calculate G moduli
omega=evans_fft(:,1);
Gp = -evans_fft(:,4:5)./(repmat(evans_fft(:,1),1,2).*evans_fft(:,6:7)); %Uncorrected for trap!
Gpp = -evans_fft(:,2:3)./(repmat(evans_fft(:,1),1,2).*evans_fft(:,6:7));

end