%Input: enter mother directory name (it must end with a "\" sign)
folder1='/Users/AleCac/Documents/Python_stuff/';

%Input: enter the last stamp of the track_output files (a list if more than
%one is present)

cell_folder_list={'RheoDataFree'};

listlen=length(cell_folder_list);

%Input: enter frame rate
fps=1000;
%Input: FFT calculation
FFT=1;

%Input: microrheo param
Dbead=1.2e-6;
kB = 1.3806503 * 1E-23;
Temp = 273+25;
t_0=5e-3;
t_end=0.5;
t_max=1;


close all
%%
for frame=1:listlen
    filename=['track_particle_output',char(cell_folder_list(frame)),'.csv'];
    try MSD=dlmread([folder1,'MSD_',filename],',');
        disp('MSD file found.')
    catch ME
        if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
            error('No MSD file found! Please run the DataAnalysis_Standard routine.')
        else
            rethrow(ME)
        end
    end
    %Conversion to SI units 
    Jfactor=1e-12*3.*pi*(Dbead/2)/(2*kB*Temp);
    indx=find(MSD(:,1)<=t_max);
    J=[MSD(indx,1),MSD(indx,2)*Jfactor];
    loglog(J(:,1),J(:,2))
    
    %Extrapolation
    indx=find(J(:,1)<=t_0);
    if length(indx)<2
        error(['Invalid timeframe t_0: ',num2str(t_0)])
    end
    x=J(indx,1);
    y=J(indx,2);
    fit=polyfit(x,y,1);
    J0=fit(2);
    
    if J0 >= 0
        disp(['Fitted J0: ',num2str(J0)])
    else
        error('Invalid J0. Please check your MSD.')
    end
    
    indx=find(J(:,1)>=t_end);
    if length(indx)<2
        error(['Invalid timeframe t_end: ',num2str(t_end)])
    end
    x=J(indx,1);
    y=J(indx,2);
    fit=polyfit(x,y,1);
    eta=1/fit(1);
    
    disp(['Fitted eta: ',num2str(eta)])
    
    %Log sampling
    interval=1;
    res=1;
    while (interval < length(J))
        res=res+1;
        interval=ceil(1.45^res);
    end
    res=res-1;
    interval=ceil(1.45.^(0:res));
    
    %IS THIS REALLY USEFUL???
    
    plotinterval=unique(floor(logspace(0,log10(length(J)),1000)));
    tau_downsample=J(interval,1);
    J_downsample=J(interval,2);
    J_spline=csape(tau_downsample,J_downsample,'second');
    
    %Plotting section
    figure(1)
    ax1=subplot(1,2,1);
    loglog(J(plotinterval,1),J(plotinterval,2),'o')
    hold on
    loglog(tau_downsample,J_downsample,'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
    loglog(tau_downsample,fnval(tau_downsample,J_spline),'--k');
    xlabel('Time (s)')
    ylabel('J(t)')
    hold off
    %%
    if FFT==0
        try J_fft=dlmread([folder1,'FFT-J_',filename],',');
            disp('FFT-J file found.')
        catch ME
            if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                error('FFT-J file not found. Re-run with the FFT calculation enabled if you wish to continue with the microrheology analysis.');
            else
                rethrow(ME)
            end
        end
    elseif FFT==1
        if J(1,1)>0
            new_J=[0,J0 ; J];
        else
            error('Invalid starting point of t! (<= 0.0)')
        end
        
        %POINT OF DEVELOPMENT!
     
        Beta=1;
        
        resfreq=fps/length(J(:,1));
        minfreq=1/(resfreq*tau_downsample(end));
        maxfreq=fps/(2);
        omega=resfreq*logspace(log10(minfreq),log10(maxfreq));
        oversample_tau=(1/(Beta*fps):1/(Beta*fps):tau_downsample(end));
       %%

        %Enable progress bar for parallel pool
        try
            parpool;
        catch ME
            if ~strcmp(ME.identifier,'parallel:convenience:ConnectionOpen')
                rethrow(ME)
            end
        end
        warning('off','MATLAB:Java:DuplicateClass')
        pctRunOnAll javaaddpath java
        progressStep=ceil(length(omega)/100);
        
        disp('Analysis in progress...')
        oversample_J=fnval(oversample_tau,J_spline);
        
        fft_J = 1i*2*pi*omega;
        fft_J = fft_J + (1-exp(-1i*2*pi*omega*oversample_tau(1)))*(oversample_J(1)-1)/oversample_tau(1);
        ppm = ParforProgMon('Progress: ', length(omega), progressStep, 300, 80);
        
        parfor i=1:length(omega)
            fft_J(i)=fft_J(i)-sum(diff(exp(-1i*2*pi*omega(i)*oversample_tau)).*diff(oversample_J)./diff(oversample_tau));
            if mod(i,progressStep)==0
                ppm.increment();
            end
        end
        ppm.delete()
        clear oversample_J
        J_fft = [omega',real(fft_J)',imag(fft_J)'];
        J_fft(:,2) = J_fft(:,2)./(-4*pi*pi*J_fft(:,1).^2);
        J_fft(:,3) = J_fft(:,3)./(-4*pi*pi*J_fft(:,1).^2);
        clear fft_J
        
        dlmwrite([folder1,'FFT-J_',filename],J_fft, 'delimiter', ',', 'precision', 9);
    else
        error('Wrong FFT parameter!')
    end
    J_fft(:,2)=-J_fft(:,2);
    J_fft(:,3)=-(J_fft(:,3) + 1./(2*pi*J_fft(:,1)));
    
    ampl=J_fft(:,2).^2 + J_fft(:,3).^2;
    omega=J_fft(:,1);
    Gp = -Jfactor(1)*(J_fft(:,3)./(2*pi*J_fft(:,1).*ampl(:,1))); %Uncorrected for trap for visualization purposes
    Gpp = -Jfactor(1)*J_fft(:,2)./(2*pi*J_fft(:,1).*ampl(:,1));
    clear J_fft
    %%
    ax2=subplot(1,2,2);
    semilogx(omega,Gp,'*')
    hold on
    semilogy(omega,Gpp,'x')
    hold off
    xlabel('Frequency(Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''','Location','northwest')
    dlmwrite([folder1,'Gmoduli_',filename],[omega,Gp,Gpp], 'delimiter', ',', 'precision', 9);

end