%This routine executes a single-particle tracking for bright field images
%in which the particle has a doughnut-like appearance. It applies an
%2D autocorrelation filter with a custom filter input by the user and
%paraboloid fitting to get sub-pixel resolution (see AnalyseTrack_Movie.m
%documentation for more details).
%It DOESN'T require pre-tracking procedures.
%This version processes serially multiple files.
close all

%Input: enter mother directory name (it must end with a "\" sign)
folder1='J:\DNAHydrogel_ParticleTracking_2019_01_24\L0_1.5%\';


%Input: tag of the desired movie files (type '*' to select all)
folder2='*';

list = dir([folder1,folder2,'.movie']);

for i=1:length(list)
    filename=[folder1,list(i).name];
    disp(list(i).name);
    [fileframes,filesize,filedata,fileheader]=movie2info(filename);
    
    disp(['File is ',Bytes2str(filesize)]);
    disp(['Number of frames: ',num2str(fileframes)]);
    
    %Change the offset if you wish to start from a different frame (but
    %you need to know its binary offset. It is: 
    %offset = #startframe*filedata + fileheader
    offset=0;
    img=squeeze(movie2frame(filename,offset,1));
    
    [filter{1,1},filter{1,2},r_max]=CreateFilter(img);
    
    
    %The movie2frame routine processes all the file by default. Modify as
    %follows if you wish to process only part of it:
    %1) Insert the desired nfr as a third argument of the first call to
    %   movie2frame
    %2) Insert nfr as a third argument of all the subsequent calls to
    %   movie2frame
    [video,offset,nfr]=movie2frame(filename,offset);
    [poslist,filter_list]=AnalyseTrack_Movie(video,filter,r_max*2,10);
    [video,offset,nfr]=movie2frame(filename,offset);
    while length(video)>1
        [poslist,filter_list]=AnalyseTrack_Movie(video,filter,r_max*2,10,poslist);
        [video,offset,nfr]=movie2frame(filename,offset);
    end
    
    %Input: by default, the program will create an output file in the mother
    %directory. Modify if necessary
    fullFileName=fullfile(folder1,['track_particle_output',regexprep(list(i).name,'.movie',''),'.csv']);
    if exist(fullFileName,'file')==2
        %Specify how to deal with pre-existing file
        loop='A';
        while (loop ~= 'Y') && (loop ~= 'N')
            loop=input('Output track file already exists. Do you want to overwrite it? [Y/N] ','s');
            if loop=='Y'
                delete(fullFileName);
            elseif loop=='N'
                while exist(fullFileName,'file')==2
                    fullFileName = [fullFileName(1:end-4),'bis.csv'];
                end
            else
                disp('Incorrect input from user.')
                loop = 'A';
            end
        end
    end
    %Write the poslist to file
    dlmwrite(fullFileName, poslist, 'delimiter', ',', 'precision', 9);
    
    %Plot trajectory at the end
    figure(2)
    plot(poslist(:,1),poslist(:,2))
    hold on
    plot(poslist(:,1),poslist(:,3))
    hold off
    title(['File ',num2str(i),': Coordinates over time'])
    xlabel('Time step')
    ylabel('Pixel')
    legend('X coordinate','Y coordinate')
end


