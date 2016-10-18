%Input: enter mother directory name (it must end with a "\" sign)
folder1='/Users/AleCac/Documents/Python_stuff/';

%Input: enter the last stamp of the track_output files (a list if more than
%one is present)

cell_folder_list={'RheoDataFree.csv'};

listlen=length(cell_folder_list);

%Input: enter frame rate
fps=1000;

%Input: FFT calculation
FFT=0;

%Input: microrheo param
Dbead=1.2e-6;
kB = 1.3806503 * 1E-23;
Temp = 273+25;
t_0=1e-1;
t_end=50;


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
    J=[MSD(:,1),MSD(:,2)*Jfactor];
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
    
    
    plotinterval=unique(floor(logspace(0,log10(length(J)),1000)));
    tau_downsample=J(interval,1);
    J_downsample=J(interval,2);
    J_spline=csape(tau_downsample,J_downsample,'second');
    
    %Plotting section
    figure(1)
    loglog(J(plotinterval,1),J(plotinterval,2),'o')
    hold on
    loglog(tau_downsample,J_downsample,'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
    loglog(tau_downsample,fnval(tau_downsample,J_spline),'--k');
    xlabel('Time (s)')
    ylabel('J(t)')
    hold off
    %%
    if FFT==0
        try acf_fft=dlmread([folder1,'FFT-J_',filename],',');
            disp('FFT-J file found.')
        catch ME
            if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                error('FFT-J file not found. Re-run with the FFT calculation enabled if you wish to continue with the microrheology analysis.');
            else
                rethrow(ME)
            end
        end
    elseif FFT==1
        
        Beta=5000;
        
        resfreq=fps/length(J(:,1));
        minfreq=1/(resfreq*tau_downsample(end));
        maxfreq=fps/(2*resfreq);
        omega=resfreq*logspace(log10(minfreq),log10(maxfreq));
        oversample_tau=(1/(Beta*fps):1/(Beta*fps):tau_downsample(end));
       

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
        
        %X analysis
        disp('Analysis of X direction in progress...')
        oversample_acf=fnval(oversample_tau,acfx_spline);
        
        fft_acf = 1i*2*pi*omega;
        fft_acf = fft_acf + (1-exp(-1i*2*pi*omega*oversample_tau(1)))*(oversample_acf(1)-1)/oversample_tau(1);
        ppm = ParforProgMon('Progress: ', length(omega), progressStep, 300, 80);
        
        parfor i=1:length(omega)
            fft_acf(i)=fft_acf(i)-sum(diff(exp(-1i*2*pi*omega(i)*oversample_tau)).*diff(oversample_acf)./diff(oversample_tau));
            if mod(i,progressStep)==0
                ppm.increment();
            end
        end
        ppm.delete()
        clear oversample_acf
        acf_fft = [omega',real(fft_acf)',imag(fft_acf)'];
        acf_fft(:,2) = acf_fft(:,2)./(-4*pi*pi*acf_fft(:,1).^2);
        acf_fft(:,3) = acf_fft(:,3)./(-4*pi*pi*acf_fft(:,1).^2);
        clear fft_acf

        %Y analysis
        disp('Analysis of Y direction in progress...')
        oversample_acf=fnval(oversample_tau,acfy_spline);
        
        fft_acf = 1i*omega*2*pi;
        fft_acf = fft_acf + (1-exp(-1i*2*pi*omega*oversample_tau(1)))*(oversample_acf(1)-1)/oversample_tau(1);
        ppm = ParforProgMon('Progress: ', length(omega), progressStep, 300, 80);
        
        parfor i=1:length(omega)
            fft_acf(i) = fft_acf(i)-sum(diff(exp(-1i*2*pi*omega(i)*oversample_tau)).*diff(oversample_acf)./diff(oversample_tau));
            if mod(i,progressStep)==0
                ppm.increment();
            end
        end
        ppm.delete()
        clear oversample_tau 
        clear oversample_acf
        acf_fft(:,4:5) = [real(fft_acf)',imag(fft_acf)'];
        acf_fft(:,4)=acf_fft(:,4)./(-4*pi*pi*acf_fft(:,1).^2);
        acf_fft(:,5)=acf_fft(:,5)./(-4*pi*pi*acf_fft(:,1).^2);
        clear fft_acf
        
        dlmwrite([folder1,'FFT-ACF_',filename],acf_fft, 'delimiter', ',', 'precision', 9);
    else
        error('Wrong FFT parameter!')
    end
    acf_fft(:,[2,4])=-acf_fft(:,[2,4]);
    acf_fft(:,[3,5])=-(acf_fft(:,[3,5]) + 1./(2*pi*acf_fft(:,[1,1])));
    
    ampl=acf_fft(:,[2,4]).^2 + acf_fft(:,[3,5]).^2;
    omega=acf_fft(:,1);
    Gpx = -Jfactor(1)*(acf_fft(:,3)./(2*pi*acf_fft(:,1).*ampl(:,1))); %Uncorrected for trap for visualization purposes
    Gppx = -Jfactor(1)*acf_fft(:,2)./(2*pi*acf_fft(:,1).*ampl(:,1));
    Gpy = -Jfactor(2)*(acf_fft(:,5)./(2*pi*acf_fft(:,1).*ampl(:,2))); %Uncorrected for trap for visualization purposes
    Gppy = -Jfactor(2)*acf_fft(:,4)./(2*pi*acf_fft(:,1).*ampl(:,2));
    clear acf_fft
    %%
    ax3=subplot(2,2,2);
    loglog(omega,Gpx,'*')
    hold on
    loglog(omega,Gppx,'x')
    loglog(omega,Jfactor(1)*ones(length(omega),1),'--k')
    hold off
    xlabel('Frequency(Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''',['Trap: ',num2str(Jfactor(1),3),' Pa'],'Location','northwest')
    
    ax4=subplot(2,2,4);
    loglog(omega,Gpy,'*')
    hold on
    loglog(omega,Gppy,'x') 
    loglog(omega,Jfactor(2)*ones(length(omega),1),'--k')
    hold off
    xlabel('Frequency (Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''',['Trap: ',num2str(Jfactor(2),3),' Pa'],'Location','northwest')
    linkaxes([ax3,ax4],'x')
    
    dlmwrite([folder1,'Gmoduli_',filename],[omega,Gpx,Gpy,Gppx,Gppy], 'delimiter', ',', 'precision', 9);


end