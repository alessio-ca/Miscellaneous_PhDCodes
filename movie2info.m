function [num_frames,num_bytes,size_data,size_header]=movie2info(filename)
%Get info of a .movie file. 
%Read operations are performed in binary.
%
%INPUT: 
%filename:          name of .movie file
%
%OUTPUT:
%num_frames:    total number of frames in video     
%num_bytes:     total number of bytes
%size_data:     size of a single frame (header + image data) (in bytes)
%size_header:   size of the .movie file header (in bytes)

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


%Open file
fileID=fopen(filename,'r');
if fileID < 0
    error('File not found.')
end
status=fseek(fileID, 0, -1 );
if status~=0
    fclose(fileID);
    error('File is not valid. EOF reached.')
end

%Read file header
offset=0;
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
size_header=offset;
if found
    %First frame loading.
    %Preliminary routine
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
        
        %Read frame header
        if (camera_frame(3)==CAMERA_TYPE_IIDC)
            [~,msgframe]=fread(fileID,43,'uint32','l');
            %Check on frame validity
            if msgframe~=43
                offset=ftell(fileID);
                fclose(fileID);
                error(['Corrupted header at offset ',num2str(offset)])
            end
            
        else
            fclose(fileID);
            error('Non-defined camera.')
        end
        
        %Output quantities
        num_bytes=dir(filename);
        num_bytes=num_bytes.bytes;
        size_data=camera_frame(5)+camera_frame(6);
        num_frames=(num_bytes - size_header)/size_data;
        
        if round(num_frames)~= num_frames
            error('Data is corrupted.')
        end
        
    else
        fclose(fileID);
        error('Data is corrupted.')
    end
    
else
    error('File is not valid. EOF reached.')
end

fclose(fileID);
end
