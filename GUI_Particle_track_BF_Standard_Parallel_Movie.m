%This routine executes a single-particle tracking for bright field images
%in which the particle has a doughnut-like appearance. It applies an
%2D autocorrelation filter with a custom filter input by the user and
%paraboloid fitting to get sub-pixel resolution (see AnalyseTrack_Movie.m
%documentation for more details).
%It DOESN'T require pre-tracking procedures.
%It supports parallelization.
%It DOESN'T adapt the filter (required for parallelization). Frames which
%can't be tracked are lost (reported at the end of the execution).

%It loads the frames from the .movie file in the 3D array "video". By default, 
%it stores only a number of frames specified by "maxframes" in movie2frame.m
%in order to control memory allocation.
%If the total number of frames is larger, video is flushed after
%"maxframes". 


%Filter syntax:
%rMax,rMid and rMin are mandatory and must be positive
%int1,int2 and int3 are optional
%Background intensity is always 1
%Outer rim is int1+1
%Mid rim is 1+int1-int2
%Inner rim is 1+int1-int2+int4
%
%Default is: 1 (outer rim = 2)
%            2 (mid rim = 0)
%            4 (inner rim = 4)
%
%If you wish an inverted image (black center and bright corona), simply
%type negative intensities (ex: -1 -2 -4 produces an inverted filter of
%the default one)



close all

%Input: enter file name
[filename,fullpath]=uigetfile('*.movie','Select a file','/Users/AleCac/Desktop/');


[fileframes,filesize,filedata,fileheader]=movie2info([fullpath,filename]);

disp(['File is ',Bytes2str(filesize)]);
disp(['Number of frames: ',num2str(fileframes)]);


%Change the offset if you wish to start from a different frame (but
%you need to know its binary offset. It is: 
%offset = #startframe*filedata + fileheader
offset=0;
img=squeeze(movie2frame([fullpath,filename],offset,1)); 
imshow(mat2gray(img));

[filter{1},r_max]=CreateFilter(img);

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
progressStep=5000;
ppm = ParforProgMon('Progress: ', fileframes, progressStep, 300, 80);
disp('Tracking in progress...');


%The movie2frame routine processes all the file by default. Modify as
%follows if you wish to process only part of it:
%1) Insert the desired nfr as a third argument of the first call to
%   movie2frame
%2) Insert nfr as a third argument of all the subsequent calls to
%   movie2frame
[video,offset,nfr]=movie2frame([fullpath,filename],offset);
[poslist,histostep]=AnalyseTrack_Parallel_Movie(video,filter,ppm,progressStep,r_max*2);
figure(1)
histogram(histostep);
[video,offset,nfr]=movie2frame([fullpath,filename],offset);
while length(video)>1
    [poslist,histostep]=AnalyseTrack_Parallel_Movie(video,filter,ppm,progressStep,r_max*2,poslist); 
    [video,offset,nfr]=movie2frame([fullpath,filename],offset);
end
ppm.delete()

%Input: by default, the program will create an output file in the mother
%directory. Modify if necessary
fullFileName=fullfile(fullpath,['track_particle_output',filename,'.csv']);
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
title('Coordinates over time')
xlabel('Time step')
ylabel('Pixel')
legend('X coordinate','Y coordinate')


