function avar=allan(data, fps, tau)

% Compute various Allan deviations for a constant-rate time series
% [osig,osigerr,tau1,tauerr]=allan(DATA, TAU) 
%
% INPUTS:
%  DATA    the time series measurements in arb. units
%  FPS     constant rate of time series in (Hz)
%  TAU is an array of the tau values for computing Allan deviations
%
% OUTPUTS: 
%  osig    = Allan deviation with overlapping estimate
%  osigerr = standard error of overlapping Allan deviation
%  tau1    = measurement interval in (s)
%  tauerr  = errors in tau that might occur because of initial rounding
%
% NOTES:
%
% No pre-processing of the data is performed.
% For constant-rate time series, the deviations are only calculated for tau
% values greater than the minimum time between samples and less than half
% the total time.


n=length(data); %Length of data series
jj=length(tau); %Length of tau series
m=floor(tau*fps); %Indeces of points used for allan variance calculation

osig    = zeros(jj,1);
osigerr = zeros(jj,1);


tic;

for j=1:jj
    %Nomenclature: 
    %   m(j) is the sample period (i.e. how many acquisitions separated two points)
    %   D is a vector containing, for each m(j), the data series averaged over a period of m(j) starting at each different acquisition.   
    
    cs = cumsum(data);
    D = (cs(m(j):end) - [0;cs(1:end-m(j))])./m(j);
    
    %overlapping Allan deviation
    z1=D(m(j)+1:n+1-m(j));
    z2=D(1:n+1-2*m(j));
    u=sum((z1-z2).^2);
    osig(j)=sqrt(u/(n+1-2*m(j))/2);
    osigerr(j)=osig(j)/sqrt(n-m(j));
        
end;

tau1=m'/fps;
tauerr=tau'-tau1;

avar=[tau1,osig,osigerr,tauerr];

toc;
end