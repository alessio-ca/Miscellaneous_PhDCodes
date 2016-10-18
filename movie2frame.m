function [video,current_offset,nfr_residual]=movie2frame(filename,offset_tag,nfr)
%Converts a .movie file into a 3D array of frames, with "maxframes" number of
%frames stored in memory. 
%Read operations are performed in binary.
%
%INPUT: 
%filename:          name of .movie file
%offset:            starting offset for read operations (optional, default=0)
%nfr:               number of frames to read (optional, default=all)
%
%OUTPUT:
%video:             3D array storing the last "maxframes" frames read
%current_offset:    current position of file cursor for read operations
%nfr_residual:      residual number of frames to be read

maxframes=30000;

%Termination rule
if nargin == 3
    if nfr < 0
        error('nfr must be >=0!')
    elseif nfr > maxframes
        disp(['Nfr exceeds the number of frames which can be temporarily stored. The video output will be flushed after ',num2str(maxframes)]);
        nfr_residual=nfr-maxframes;
        nfr = maxframes;
    else
        nfr_residual=0;
    end
else
    nfr=maxframes;
    nfr_residual=maxframes;
end
if nargin < 2
    offset_tag=0;
end



%Camera values
CAMERA_MOVIE_MAGIC=1231906132;
CAMERA_MOVIE_VERSION=1;
CAMERA_TYPE_IIDC=1;
CAMERA_PIXELMODE_MONO_8=1;
CAMERA_PIXELMODE_MONO_16BE=2;
CAMERA_PIXELMODE_MONO_16LE=3;

%IIDC defines
IIDC_MOVIE_HEADER_LENGTH=172;
%Feature modes
IIDC_FEATURE_MODE_OFF=bin2dec('1');
IIDC_FEATURE_MODE_RELATIVE=bin2dec('10');
IIDC_FEATURE_MODE_ABSOLUTE=bin2dec('100');
IIDC_FEATURE_MODE_AUTO=bin2dec('1000');
IIDC_FEATURE_MODE_ONEPUSH=bin2dec('10000');
IIDC_FEATURE_MODE_ADVANCED=bin2dec('100000');
%Trigger
IIDC_TRIGGER_INTERNAL = -1;
IIDC_TRIGGER_EXTERNAL0 = 0;
IIDC_TRIGGER_EXTERNAL1 = 1;
IIDC_TRIGGER_EXTERNAL15 = 7;

offset=offset_tag;

%Open file
fileID=fopen(filename,'r');
status=fseek(fileID, offset, -1 );
if status~=0
    fclose(fileID);
    error('Offset value is not valid. EOF reached.')
end

%Read file header.
found=0;
[A,msg]=fread(fileID,1,'uint32','l');
while msg==1
    if A==CAMERA_MOVIE_MAGIC
        found = 1;
        break
    end
    offset=offset+1;
    fseek(fileID, offset, -1 );
    [A,msg]=fread(fileID,1,'uint32','l');
end
fseek(fileID, offset, -1 );

if found
    %Frame loading.
    %Preliminary routine on first frame.
    [camera_frame,msg]=fread(fileID,6,'uint32','l');
    %Camera_frame is a vector with entries:
    %   1 magic offset
    %   2 camera version
    %   3 camera type
    %   4 pixel mode
    %   5 header data (single frame)
    %   6 total data (single frame)
    
    if msg==6
        %Consistency checks
        if camera_frame(1)~=CAMERA_MOVIE_MAGIC
            offset=ftell(fileID);
            fclose(fileID);
            error(['Wrong magic at offset ',num2str(offset)])
        end
        if camera_frame(2)~=CAMERA_MOVIE_VERSION
            fclose(fileID);
            error(['Unsupported movie version ',num2str(A(2))])
        end
        %Rewind frame header
        fseek(fileID, -24, 0 );
        
        %Read first frame header
        if (camera_frame(3)==CAMERA_TYPE_IIDC)
            [A,msgA]=fread(fileID,6,'uint32','l');
            [B64,msgB64]=fread(fileID,1,'uint64','l');
            [B32,msgB32]=fread(fileID,2,'uint32','l');
            [C,msgC]=fread(fileID,2,'uint32','l');
            [D,msgD]=fread(fileID,1,'uint64','l');
            [E,msgE]=fread(fileID,6,'uint32','l');
            [F,msgF]=fread(fileID,3,'uint32','l');
            [G32,msgG32]=fread(fileID,1,'uint32','l');
            [G64,msgG64]=fread(fileID,1,'uint64','l');
            [H,msgH]=fread(fileID,17,'uint32','l');
            msgframe=msgA+msgB64+msgB32+msgC+msgD+msgE+msgF+msgG32+msgG64+msgH;
            
            %Check on frame validity
            if msgframe~=40
                offset=ftell(fileID);
                fclose(fileID);
                error(['Corrupted header at offset ',num2str(offset)])
            end
            
            %Size in pixels
            size_x=E(3);
            size_y=E(4);
            
        else
            fclose(fileID);
            error('Non-defined camera.')
        end
        
        fseek(fileID, offset, -1 );
        
    else
        fclose(fileID);
        error('Data is corrupted.')
    end
    
    %Allocate memory for the 3D array
    if (camera_frame(4)==CAMERA_PIXELMODE_MONO_8)
        video=zeros(nfr,size_y,size_x,'uint8');
    else
        video=zeros(nfr,size_y,size_x,'uint16');
    end
    
    %Read routine on frames
    index=1;
    [camera_frame,msg]=fread(fileID,6,'uint32','l');
    %Loop on frames
    while msg==6 && index<=nfr
        %Consistency checks
        if camera_frame(1)~=CAMERA_MOVIE_MAGIC
            offset=ftell(fileID);
            fclose(fileID);
            error(['Wrong magic at offset ',num2str(offset)])
        end
        if camera_frame(2)~=CAMERA_MOVIE_VERSION
            fclose(fileID);
            error(['Unsupported movie version ',num2str(A(2))])
        end
        %Rewind frame header
        fseek(fileID, -24, 0 );
        
        [~,msgframe]=fread(fileID,43,'uint32','l');
        %Check on frame validity
        if msgframe~=43
            offset=ftell(fileID);
            fclose(fileID);
            error(['Corrupted header at offset ',num2str(offset)])
        end
        
        %Get frame image
        if (camera_frame(4)==CAMERA_PIXELMODE_MONO_8)
            [data,~]=fread(fileID,[size_x,size_y],'uint8');
            data=uint8(data);
        else
            [data,~]=fread(fileID,[size_x,size_y],'uint16');
            data=uint16(data);
        end
        
        %Store video in the 3D array
        video(index,:,:)=data';
        current_offset=ftell(fileID);
        
        %Loop advance
        index=index+1;
        [camera_frame,msg]=fread(fileID,6,'uint32','l');
    end
    
    %Closure conditions
    index=index-1;
    if index < nfr
        video(index+1:nfr,:,:)=[];
    end
    if nfr_residual == 0
        fseek(fileID, 0, 1 );
        current_offset=ftell(fileID);
    end
        
else
    disp('EOF reached.')
    video=-1;
    current_offset=offset_tag;
end

fclose(fileID);
end
