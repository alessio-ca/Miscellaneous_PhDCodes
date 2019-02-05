function CGMuRheoF(N,F,M,scale,Dim,ShutT,srate,Temp,R,SwitchFifo,ParticleN,filename,header,out,Counter)

%CGMuRheoF.m
%Taiki Yanagishima & Alessio Caciagli, OE (formerly BSS), Cavendish Labs, University of Cambridge
%
%One-line function form of CGMuRheo
%
%Calculates a multiple-tau velocity autocorrelation function from
%Brownian displacement data, and converts to a complex modulus using
%a Fourier-Laplace transform scheme proposed by Evans (PRE 2010)
%
%Inputs are (N,F,M,scale,Dim,ShutT,srate,Temp,R,SwitchFifo)
%N - number of points per array
%F - factor by which arrays are coarse-grained
%M - number of arrays
%scale - factor to convert displacement series to physical units
%Dim - number of columns of data after index and time
%   Data should be formatted
%ShutT  - Shutter exposure time for camera, if visualising
%       - can be set to 0
%srate - Sampling rate
%Temp - Temperature
%R - radius of probe
%SwitchFifo - Set to 1 if reading from a named pipe, 0 otherwise
%ParticleN - index of particle to be studied.
%header - Set to 1 if your file has a one-line header (more is NOT allowed)
%Counter - optional for batch analysis

if nargin < 15
    Counter='';
end

fclose('all');
file = fopen(filename,'r');

if file == -1
    disp('Cannot open file, check address!');
    return;
end


%Initialise arrays to accelerate calculation
format long e
Q = zeros(M);

LevelX = zeros(M,N);
LevelY = zeros(M,N);
CorrX = zeros(N,M);
CorrY = zeros(N,M);
CorrScaledX = zeros(N,M);
CorrScaledY = zeros(N,M);
CorrTermsN = zeros(1,M);

x = zeros(1,4);
Xprev = zeros(1,4);
Count = 0;
% c = 0;
% Snap = 1;
Valid = 0;
STARTED = 0;

%Empty arrays for lag, omega, Real and Imag G
lag = zeros(N,M);
Omega = zeros(N/2,M);

delT = ones(M);
timestep = 1 / srate;

%Various object handles
VxPlot = ones(M,1);
VyPlot = ones(M,1);
GxPlot = ones(2,1);
GyPlot = ones(2,1);


%Constants used in calculation
%Boltzmann constant and constant
kB = 1.3806503 * 1E-23;
const = kB*Temp/(6*pi*R);
%Important functions

%1. UpdateBaseFile
%Enters new data into the first array from a file

    function UpdateBaseFile
        
        %Next point to be filled in ungrained level
        %Index Q(1) moved forward (cyclic buffer)
        Q(1) = mod(Q(1),N) + 1;

        %First point treatment - After Xprev is set, becomes redundant.
        while Valid == 0
            if header == 1
                fgetl(file);
            end
            if Dim <= 3
                [x,Valid] = fscanf(file,'%f %f %f %f',[1,Dim+1]);
                
            else
                error('Wrong dimensionality of your system.')
            end
            
            
            Xprev = x;
           
            if STARTED == 1
                if Valid == 0
                    break;
                end
            end
            
            if Valid > 0
                STARTED = 1;
            elseif Valid == 0
                EXIT = -1;
                return;
            end
        end

        %Read in new point
        %Only fails once file read has finished.
        if Dim <= 3
            [x,Valid] = fscanf(file,'%f %f %f %f',[1,Dim+1]);
        else
            error('Wrong dimensionality of your system.')
        end
        
        if Valid == 0
            EXIT = -1;
            return;
        end
        
        
        %Update first level with new point
        if SwitchFifo == 1
            LevelX(1,Q(1)) = (x(1,2) - Xprev(1,2)) * scale;
            LevelY(1,Q(1)) = (x(1,3) - Xprev(1,3)) * scale;
        else
            LevelX(1,Q(1)) = (x(1,2 + (ParticleN -1) * 2) - Xprev(1,2 + (ParticleN -1) * 2)) * scale;
            LevelY(1,Q(1)) = (x(1,3 + (ParticleN -1) * 2) - Xprev(1,3 + (ParticleN -1) * 2)) * scale;
        end
        
        
        Xprev = x;

    end


%2. UpdateBaseFifo
%Enters new data into the first array from a file


    function UpdateBaseFifo
        
        Q(1) = mod(Q(1),N) + 1;
        
        %First point treatment
        %After Xprev is set with a non-zero Valid value,
        %becomes redundant.
        while Valid == 0
            [Ct,Valid] = fread(file,1,'uint64');
            X = fread(file,[1,2],'double');
            Xprev = X;
            disp([Ct,X])
        end
        
        %Read in values from fifo
        %Will wait for Valid return of values.
        ValidRun = 0;
        while ValidRun == 0
            [Ct,ValidRun] = fread(file,1,'uint64');
            X = fread(file,[1,2],'double');
        end
        
        LevelX(1,Q(1)) = (X(1,1) - Xprev(1,1)) * scale;
        LevelY(1,Q(1)) = (X(1,2) - Xprev(1,2)) * scale;
        
        testX = abs(log(X(1,1)) - log(Xprev(1,1)));
        testY = abs(log(X(1,2)) - log(Xprev(1,2)));
        
        if testX > 1
            disp('Pipe Error in X')
        end
        
        if testY > 1
            disp('Pipe Error in Y')
        end
        
        Xprev = X;
        
    end


%3. CoarseGrain(L)
%Coarse grains data in level L into level L+1

    function CoarseGrain(L)
        
        Q(L+1) = mod(Q(L+1),N) + 1;
        
        LevelX(L + 1,Q(L+1)) = sum(LevelX(L,Q(L)-int32(N/F)+1:Q(L)));
        LevelY(L + 1,Q(L+1)) = sum(LevelY(L,Q(L)-int32(N/F)+1:Q(L)));
        
    end


%4. CorrUpdate(L)
%Updates the level L correlation function

    function CorrUpdate(L)
        
        for j = 1:N
            CorrX(j,L) = CorrX(j,L) ...
                + LevelX(L,mod(Q(L) - N,N)+1)*LevelX(L,mod(Q(L) - N + (j - 1),N) + 1);
            
            CorrY(j,L) = CorrY(j,L) ...
                + LevelY(L,mod(Q(L) - N,N)+1)*LevelY(L,mod(Q(L) - N + (j - 1),N) + 1);
        end
        
        CorrTermsN(L) = CorrTermsN(L) + 1;
        CorrScaledX(:,L) = CorrX(:,L)/(CorrTermsN(L)*(delT(L)^2));
        CorrScaledY(:,L) = CorrY(:,L)/(CorrTermsN(L)*(delT(L)^2));
        
    end


%5. FLMulti
%Carries out the Fourier-Laplace transform following Evans et al (2010)
%To explain the range of arrays here:
%   Corr is the original, multiple tau set of correlation functions
%   Om is the corresponding set of conmensurate angular frequencies
%   Cc is the concatenated full correlation function
%   Ct is the time axis for Cc
%   COmega are the angular frequencies conmensurate with Cc.
%       Note that Cc is not evenly spaced.
%   Cw is the FT of Cc (Cwpart are increments in the loop)
%   Dw is the FT of the frequency dependent diffusion 'constant' D
%   Gstar is the complex modulus.
%
%   L is the level to be transformed
%   Lstart is the level to start with.
%       i.e. the function can be set to omit the first few levels.

    function [GMod] = FLMulti(L,XY,Lstart)
        
        
        if XY == 1
            Corr = CorrScaledX(:,Lstart:L);
        else
            Corr = CorrScaledY(:,Lstart:L);
        end
        
        lagCorr = lag(:,Lstart:L);
        Om = Omega(:,Lstart:L);
        
        Cc = zeros(N + (L-1)*(N - F),1);
        Ct = zeros(N + (L-1)*(N - F),1);
        COmega = zeros((N/2-(F/2))*L,1);
        Gstar = zeros((N/2-(F/2))*L,1);
        
        Cc(1:N,1) = Corr(1:N,1);
        Ct(1:N,1) = lagCorr(1:N,1);
        COmega((L-1)*(N/2-F/2)+1:L*(N/2-F/2),1) = Om(F/2+1:N/2,1);
        
        if L > 1
            
            for i = 1:L-1
                
                Cc(N + 1 + (N - F) * (i - 1): N + ( N - F )*i) = Corr(F+1:N,i+1);
                Ct(N + 1 + (N - F) * (i - 1): N + ( N - F )*i) = lagCorr(F+1:N,i+1);
                COmega((L - i - 1)*(N/2 - F/2) + 1: ( N/2 - F/2 )*(L-i)) = Om(F/2+1:N/2,i+1);
                
            end
            
        end
        
        %Account for delta like first point in correlation
        C0Corr = Corr(1,1) - Corr(2,1) ;
        Cc(1,1) = Cc(2,1);
        
        
        %ONLY if D(t --> 0) condition is explicitly applied
        %Not recommended for accurate eta measurements
        %Constant correction - Make sure D = 0 at long time
        %Put the first point equal to second point AFTER constant correction
        
        % DCorrect = Trapz(Ct,Cc);
        % dt = Ct(2,1) - Ct(1,1);
        % C0Corr = Corr(1,1) - (2 * DCorrect / dt) - Corr(2,1) ;
        % Cc(1,1) = Cc(2,1); % + (2 * DCorrect / dt);
        
        
        Nw = size(COmega,1);
        Nt = size(Cc);
        Dw = zeros(Nw,1);
        Cw = zeros(Nw,1);
        
        
        for k = 1 : Nw
            
            for j =  2 : Nt - 1
                
                Cwpart = ((Cc(j+1) - Cc(j))/(Ct(j+1)-Ct(j)) ...
                    - (Cc(j) - Cc(j-1))/(Ct(j) - Ct(j-1)))*exp(- Ct(j) * COmega(k) * 1i) / (-(COmega(k))^2);
                
                Cw(k,1) = Cw(k,1) + Cwpart;
                
            end
            
            Cw(k,1) = Cw(k,1) - 1i * Cc(1,1) / COmega(k,1) + 0.5 * C0Corr * Ct(2,1);
            Dw(k,1) = Cw(k,1) / (1i * COmega(k,1));
            Gstar(k,1) = const / Dw(k,1);
            
        end
        
        %Account for finite shutter speed, when necessary.
        if ShutT > 0
            Gstar = Gstar .* (sin(COmega .* ShutT / 2) ./ (COmega * ShutT / 2)).^2;
        end
        
        GMod = [COmega real(Gstar) imag(Gstar)];
        
    end


%6. subtitle
%Creates a overall central title for a group of subplots
%Returns a handle to the title and the handle to an axis.
%[ax,h]=subtitle(text)
%returns handles to both the axis and the title.
%ax=subtitle(text)
%returns a handle to the axis only.

    function [ax,h]=subtitle(text)

        ax=axes('Units','Normal','Position',[.075 .075 .85 .85],'Visible','off');
        
        set(get(ax,'Title'),'Visible','on')
        
        title(text);
        
        if (nargout < 2)
            
            return
            
        end
        
        h=get(ax,'Title');
    end
        
        
       

%Now set up graphs necessary for plotting functions.

warning off all
set(0,'defaultaxesfontsize',15)

figure %It was figure(1)
subplot(2,2,1);
axis auto
ColorV = get(gca,'ColorOrder');
set(gca,'Xscale','log','Yscale','linear');
xlabel('t (s)')
ylabel('C(t)(m^2s^{-1})')
title('C_x(t)')

subplot(2,2,3);
axis auto
set(gca,'Xscale','log','Yscale','linear');
xlabel('t (s)')
ylabel('C(t) (m^2s^{-1})')
title('C_y(t)')

subplot(2,2,2);
axis auto
ColorG = get(gca,'ColorOrder');
set(gca,'Xscale','log','Yscale','log');
xlabel('\omega (Hz)')
ylabel('G'',G" (Pa)')
title('G^*_x')

subplot(2,2,4);
axis auto
set(gca,'Xscale','log','Yscale','log');
xlabel('\omega (Hz)')
ylabel('G'',G" (Pa)')
title('G^*_y')

%Graphs need some initial input, can't delete root object during first run.

for l = 1:M
    subplot(2,2,1)
    VxPlot(l) = line(lag(:,l),CorrScaledX(:,l),'LineStyle','none','Marker','x');
    
    subplot(2,2,3)
    VyPlot(l) = line(lag(:,l),CorrScaledY(:,l),'LineStyle','none','Marker','x');
end

subplot(2,2,2)
GxPlot(1) = line(1,1,'LineStyle','none','Marker','x');
GxPlot(2) = line(1,1,'LineStyle','none','Marker','o');

legend(GxPlot,'G''','G"','Location','NorthWest')

subplot(2,2,4)
GyPlot(1) = line(1,1,'LineStyle','none','Marker','x');
GyPlot(2) = line(1,1,'LineStyle','none','Marker','o');

legend(GyPlot,'G''','G"','Location','NorthWest')


%Initialise lag and omega axes

for i = 1:M
    delT(i) = timestep * ((N/F) ^ (i-1));
    lag(:,i) = linspace(0,N-1,N)' * delT(i);
    Omega(:,i) = linspace(0,N/2-1,N/2)' * 2 * pi / (delT(i) * N);
end



%Program starts here:

EXIT = 1;

while EXIT > 0
    
    %500 updates of correlation functions.
    %Can be tailored to anything, 1E4 is arbitrary
    %MUST BE TAILORED TO STH LESS THAN DATA LENGTH
    
    for CountInner = 1:500
        
        Count = Count + 1;
        
        %First, update base array. Always necessary.
        if SwitchFifo == 1
            UpdateBaseFifo;
        else
            UpdateBaseFile;
        end
        
        if EXIT == 0
            break;
        end
    
        
        %Perform correlation update only if base array has been filled once.
        
        if Count > N
            CorrUpdate(1);
        end
        
        %Now get onto coarse-grained arrays.
        
        for n = 1:M-1
            
            if rem(Count,(N/F)^(n)) == 0
                %Only Coarse grain when the right number of
                %updates have been performed.
                
                CoarseGrain(n);
                
                %Corr update for these only occurs when
                %CoarseGrain is activated, and thus a new point
                %is introduced. Only perform when these have
                %also been filled at least once.
                
                if Count > N*((N/F)^(n))
                    CorrUpdate(n+1);
                end
                
            end
            
        end
    end
    

    %All data updates have been completed for VACF, now transform if
    %necessary - only carry out updates of G if plots are updating, VACFs
    %are cumulative but G, derived from them, are not.
    
    disp(Count);
    Lact = ceil(log(Count/N)/log(N/F));
    
    if Lact > M
        Lact = M;
    end
    
    GX = FLMulti(Lact,1,1);
    GY = FLMulti(Lact,2,1);
    
    dlmwrite(strcat(out,'GstarX.dat'),GX);
    dlmwrite(strcat(out,'GstarY.dat'),GY);
    
    %Now plot everything!!
    
    %figure(1) %It wasn't uncommented
    
    for l = 1:M
        
        subplot(2,2,1)
        delete(VxPlot(l))
        VxPlot(l) = line(lag(:,l),CorrScaledX(:,l),'LineStyle','none','Marker','x','Color',ColorV(mod(l,7)+1,:));
        
        subplot(2,2,3)
        delete(VyPlot(l))
        VyPlot(l) = line(lag(:,l),CorrScaledY(:,l),'LineStyle','none','Marker','x','Color',ColorV(mod(l,7)+1,:));
        
    end
    
    subplot(2,2,2)
    delete(GxPlot(1))
    GxPlot(1) = line(GX(:,1),GX(:,2),'LineStyle','none','Marker','x','Color',ColorG(1,:));
    delete(GxPlot(2))
    GxPlot(2) = line(GX(:,1),GX(:,3),'LineStyle','none','Marker','o','Color',ColorG(2,:));
    
    legend(GxPlot,'G''','G"','Location','NorthWest')
    
    
    subplot(2,2,4)
    delete(GyPlot(1))
    GyPlot(1) = line(GY(:,1),GY(:,2),'LineStyle','none','Marker','x','Color',ColorG(1,:));
    delete(GyPlot(2))
    GyPlot(2) = line(GY(:,1),GY(:,3),'LineStyle','none','Marker','o','Color',ColorG(2,:));
    
    legend(GyPlot,'G''','G"','Location','NorthWest')
    
    if nargin < 15
        subtitle('MicroRheo');
    else
        subtitle(['MicroRheo ',num2str(Counter)]);
    end
    
    dlmwrite(strcat(out,'VacfX.dat'),CorrScaledX);
    dlmwrite(strcat(out,'VacfY.dat'),CorrScaledY);
    
end


%Analysis complete
%Close file / pipe and exit

fclose('all');
disp('Analysis complete');

end