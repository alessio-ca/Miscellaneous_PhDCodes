%This routine calculates the visco-elastic modulus G* from position series
%of single bead(s). It employs a position autocorrelation (ACF) calculation
%with subsequent conversion to G moduli via a Fourier Transform routine for
%non-equally spaced points. For a description of the trategy, see Tassieri
%et al. (2012).

%The procedure is as follows:
%   1) The data series is converted into the ACF.
%   2) The ACF is sampled at logaritmic points to reduce the stastical
%      noise at large lag times.
%   3) A natural cubic spline is fitted to the log-sampled ACF.
%   4) A largeer, oversampled ACF is obtained via the spline fit.
%   5) The oversampled ACF is Fourier-Transformed to obtain G*.

%Some tecnical details:
%   - The ACF tail is improved via a smoothing procedure (shown in plots)
%   - The ACF log-sampling has a very high impact on the Gp
%     appearance at high frequencies. Recommended range is 1+ < logsampling
%     < 1.5, but you may want to play with the parameter for better-looking
%     results.
%   - The oversampling factor (Beta) is typically chosen so that Beta*fps
%     is roughly 5 MHz. Values higher than that are reported not to yield
%     better results (see Tassieri et al., (2012)).
%   - The FFT calculation is optional. If FFT is activated, the program 
%   - will save the generated FFT-ACF in an output file. If you wish to
%     just use the calculated FFT-ACF for data 
%     ACF in an output file if 
%   - Parallelization is used to speed-up the FT operation. If a parallel
%     pool is not already active, switching it ON may require some time.
%     Please DON'T PANIC! It's normal, everything is gonna be alright.
%   - The code is optimized and makes heavy use of vectorization. If you
%     (for some obscure reason) find yourself memory-constrained (99% of
%     the times it is because your oversampled ACF is too big -> reduce
%     that damn Beta and log-sample!), you can split the X and Y analysis
%     and do them serially. The code will be slower, but you can free up
%     memory used in the FT operation.
%   - The code loads a track_output*.csv file. It then looks for either
%     linked or dedrifted versions (a.k.a. post-processed data). If it does not
%     find them, it prompt you post-process your raw data. If you need to,
%     you CAN modify the routine to accept also raw data BUT DO IT IF YOU
%     KNOW WHAT YOU ARE DOING. THE CODE WANTS TO BE PEDAGOGICAL IN
%     PROMPTING YOU TO POST-PROCESS YOUR DATA BEFORE DOING MICRORHEOLOGY.
%     Post-processed data can avoid very nasty problems such as
%     code-crashing beacuse of some bad raw data or data with drift, with
%     consequent bullshit/like G* values. 
%     
%     REMEMBER: 
%     MICRORHEOLOGY IS AN ART, LIKE MUSIC. YOU CAN'T PLAY JAZZ IF YOU DON'T 
%     FIRST REFINE YOUR SKILLS. SAME HERE: YOU CAN'T DO MICRORHEOLOGY IF 
%     YOU DON'T HAVE APPROPRIATELY CLEAN DATA.


%Input: enter the last stamp of the track_output files. It can be list if more than
%one is present, but they must refer to same global settings (such as bead
%diameter, temperature, etc.).
out=uipickfiles('FilterSpec','track_particle_output*.csv','Output','Struct');
stamp=cell(1,length(out));
for i = 1:length(out)
    [fullpath,stamp{i},~]=fileparts(out(i).name);
end
fullpath=[fullpath,'/'];
if ~iscell(stamp)
    cell_folder_list={stamp};
end
listlen=length(stamp);

%Input: enter frame rate, bead diameter, temperature (in °C) and
%log-sampling parameter.
%Set security commands for wrong user input
loop=-1;
while loop<=0
    loop=input('Insert frame rate (fps): ');
    if loop>0
        fps=loop;
    else
        disp('Incorrect input from user.')
        loop = -1;
    end
end
loop=-1;
while loop<=0
    loop=input('Insert bead diameter (in um): ');
    if loop>0
        Dbead=loop*1e-6;
    else
        disp('Incorrect input from user.')
        loop = -1;
    end
end
while loop<=0
    loop=input('Insert temperature (in °C): ');
    if loop>0
        Temperature=loop+273;
    else
        disp('Incorrect input from user.')
        loop = -1;
    end
end
loop=-1;
while loop<=0
    loop=input('Insert coarse graining factor: ');
    if loop>0
        CG=loop;
        if loop<=1 
            disp('CG factor must be > 1')
            loop=-1;
        end
    else
        disp('Incorrect input from user.')
        loop = -1;
    end
end
%Input: FFT calculation (1 for yes, 0 for no)
loop=2;
while (loop ~= 0) && (loop ~= 1)
    loop=input('FFT calculation (0 for no, 1 for yes): ');
    if ~isscalar(loop)
        loop = 2;
    end
    if loop==1 || loop==0
        FFT=loop;
    else
        disp('Incorrect input from user.')
        loop = 2;
    end
end

%If FFT calculation is activated, the following parameters are also
%required:
if FFT ==1
    %Input: insert oversampling factor. 
    %Set security commands for wrong user input
    loop=-1;
    while loop<=0
        loop=input('Insert oversampling factor: ');
        if loop>0
            Beta=loop;
        else
            disp('Incorrect input from user.')
            loop = -1;
        end
    end
    %Input: insert time cut-off (set it to 0 if you don't want any time cut-off)
    %Set security commands for wrong user input
    loop=-1;
    while loop<0
        loop=input('Insert time cut-off (set it to 0 if not needed): ');
        if loop>=0
            taucutoff=loop;
        else
            disp('Incorrect input from user.')
            loop = -1;
        end
    end
end




kB = 1.3806503 * 1E-23;


close all
disp(' ')
%%
for frame=1:listlen
    filename=[char(stamp(frame)),'.csv'];
    try traj=dlmread([fullpath,'De-drifted_',filename],',');
        disp('De-drifted trajectory file found.')
    catch ME
        if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
            disp('De-drifted trajectory file not found. Proceeding with un-drifted...')
            try traj=dlmread([fullpath,'Linked_',filename],',');
                traj(:,4)=[];
                traj=circshift(traj,[0 1]);
                disp('Linked trajectory file found.')
            catch ME
                if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                    error('No linked trajectory file found! Please run the DataAnalysis_Standard routine.')
                else
                    rethrow(ME)
                end
            end
            const=traj(1,2:3);
            traj(:,2:3)=traj(:,2:3)-repmat(const,length(traj),1);
        else
            rethrow(ME)
        end
    end
    rheofactor=1e12*kB*Temp./[mean(traj(:,2).^2),mean(traj(:,3).^2)];
    rheofactor=rheofactor./(3*pi*Dbead);
    [acfx,~,~,~] = acf(1/fps,traj(:,2));
    [acfy,tau,~,~] = acf(1/fps,traj(:,3));
    
    idx=find(acfx<0, 1 );
    idy=find(acfy<0, 1 );
    xx=smooth(acfx(idx:end),6*fps);
    yy=smooth(acfy(idy:end),6*fps);
    acfxfiltered=[acfx(1:idx-1);xx];
    acfyfiltered=[acfy(1:idy-1);yy];

    

    
    maxlag=max(length(tau)-50000,50);
    maxlag=min(floor(length(tau)/4),maxlag);
    %Log sampling
    interval=1;
    res=1;
    while (interval < maxlag/4)
        res=res+1;
        interval=ceil(CG^res);
    end
    res=res-1;
    interval=ceil(CG.^(0:res));
    
    plotinterval=unique(floor(logspace(0,log10(maxlag),1000)));
    if taucutoff>0
        interval=interval(tau(interval)<taucutoff);
    end
    tau_downsample=tau(interval);
    acfx_downsample=acfxfiltered(interval);
    acfy_downsample=acfyfiltered(interval);
    acfx_spline=csape(tau_downsample,acfx_downsample);
    acfy_spline=csape(tau_downsample,acfy_downsample);
    disp('ACF calculation done.')
    disp(' ')
    
    %Plotting section
    figure(1)
    ax1=subplot(2,2,1);
    semilogx(tau(plotinterval),acfx(plotinterval),'o')
    hold on
    semilogx(tau(plotinterval),acfxfiltered(plotinterval),'o')
    semilogx(tau_downsample,acfx_downsample,'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
    h=semilogx(tau_downsample,fnval(tau_downsample,acfx_spline),'--k');
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('ACF','Smooth ACF','Log-sampling','Location','northeast')
    xlabel('Time (s)')
    ylabel('A(t)')
    hold off
    
    ax2=subplot(2,2,3);
    semilogx(tau(plotinterval,:),acfy(plotinterval,:),'o');
    hold on
    semilogx(tau(plotinterval),acfyfiltered(plotinterval),'o')
    semilogx(tau_downsample,acfy_downsample,'s','MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',10)
    h=semilogx(tau_downsample,fnval(tau_downsample,acfy_spline),'--k');
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    legend('ACF','Smoothened ACF','Log-sampling','Location','northeast')
    xlabel('Time (s)')
    ylabel('A(t)')
    hold off
    
    linkaxes([ax1,ax2],'x')

    %%
    if FFT==0
        try fft_acf=dlmread([fullpath,'FFT-ACF_',filename],',');
            disp('FFT-ACF file found.')
        catch ME
            if strcmp(ME.identifier,'MATLAB:dlmread:FileNotOpened')
                error('FFT-ACF file not found. Re-run with the FFT calculation enabled if you wish to continue with the microrheology analysis.');
            else
                rethrow(ME)
            end
        end
    elseif FFT==1
                
        resfreq=fps/length(tau);
        minfreq=1/(resfreq*tau_downsample(end));
        maxfreq=fps/(2*resfreq);
        omega=resfreq*logspace(log10(minfreq),log10(maxfreq))';
        oversample_tau=(1/(Beta*fps):1/(Beta*fps):tau_downsample(end))';
        
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
        
        %Analysis section
        disp('FT Analysis in progress...')
        oversample_acf(:,1)=fnval(oversample_tau,acfx_spline);
        oversample_acf(:,2)=fnval(oversample_tau,acfy_spline);
        
        fft_acf = repmat(1i*2*pi*omega,1,2);
        sum_point=[(1-exp(-1i*2*pi*omega*oversample_tau(1)))*(oversample_acf(1,1)-1)/oversample_tau(1),(1-exp(-1i*2*pi*omega*oversample_tau(1)))*(oversample_acf(1,2)-1)/oversample_tau(1)];
        fft_acf = fft_acf + initial_point;
        ppm = ParforProgMon('Progress: ', length(omega), progressStep, 300, 80);
        
        parfor i=1:length(omega)
            sum_point=[sum(diff(exp(-1i*2*pi*omega(i)*oversample_tau)).*diff(oversample_acf(:,1))./diff(oversample_tau)),sum(diff(exp(-1i*2*pi*omega(i)*oversample_tau)).*diff(oversample_acf(:,2))./diff(oversample_tau))];
            fft_acf(i,:)=fft_acf(i,:)-sum_point;
            if mod(i,progressStep)==0
                ppm.increment();
            end
        end
        ppm.delete()
        clear oversample_acf
        fft_acf = [omega,real(fft_acf),imag(fft_acf)];
        fft_acf(:,2:3) = fft_acf(:,2:3)./repmat(-4*pi*pi*fft_acf(:,1).^2,1,2);
        fft_acf(:,4:5) = fft_acf(:,4:5)./repmat(-4*pi*pi*fft_acf(:,1).^2,1,2);
        
        dlmwrite([fullpath,'FFT-ACF_',filename],fft_acf, 'delimiter', ',', 'precision', 9);
    else
        error('Wrong FFT parameter!')
    end
    fft_acf(:,2:3)=-fft_acf(:,2:3);
    fft_acf(:,4:5)=-(fft_acf(:,4:5) + 1./repmat(2*pi*fft_acf(:,1),1,2));
    
    fft_acf(:,6:7)=fft_acf(:,2:3).^2 + fft_acf(:,4:5).^2;
    omega=fft_acf(:,1);
    Gpx = -rheofactor(1).*fft_acf(:,2)./(2*pi*fft_acf(:,1).*fft_acf(:,6)); %Uncorrected for trap for visualization purposes
    Gpy = -rheofactor(2).*fft_acf(:,3)./(2*pi*fft_acf(:,1).*fft_acf(:,7)); %Uncorrected for trap for visualization purposes
    Gppx = -rheofactor(1).*fft_acf(:,4)./(2*pi*fft_acf(:,1).*fft_acf(:,6)); 
    Gppy = -rheofactor(2).*fft_acf(:,5)./(2*pi*fft_acf(:,1).*fft_acf(:,7));
    clear fft_acf
    %%
    ax3=subplot(2,2,2);
    loglog(omega,Gpx,'*')
    hold on
    loglog(omega,Gppx,'x')
    loglog(omega,rheofactor(1)*ones(length(omega),1),'--k')
    hold off
    xlabel('Frequency(Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''',['Trap: ',num2str(rheofactor(1),3),' Pa'],'Location','northwest')
    
    ax4=subplot(2,2,4);
    loglog(omega,Gpy,'*')
    hold on
    loglog(omega,Gppy,'x') 
    loglog(omega,rheofactor(2)*ones(length(omega),1),'--k')
    hold off
    xlabel('Frequency (Hz)')
    ylabel('G'',G'''' (Pa)')
    legend('G''','G''''',['Trap: ',num2str(rheofactor(2),3),' Pa'],'Location','northwest')
    linkaxes([ax3,ax4],'x')
    
    dlmwrite([fullpath,'Gmoduli_',filename],[omega,Gpx,Gpy,Gppx,Gppy], 'delimiter', ',', 'precision', 9);


end