function movie2tiff

desiredFolder = input('Folder: ');
cd(desiredFolder);
desiredFile = input('Movie files (type ''*'' to select all): ');
selectedFiles = dir([desiredFile '.movie']);
if isempty(selectedFiles)
    error('No input files.')
end
filenames={selectedFiles.name};
for i=1:length(filenames)
   movie2tiff_func(filenames{1});
end


function movie2tiff_func(filename)

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


%Create directory
myfolder=regexprep(filename,'.movie','');   
mkdir(myfolder);
disp(myfolder);

%Open movie file
fileID=fopen(filename,'r');

%Find the beginning of binary data
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

%Go to the beginning of frames
fseek(fileID, offset, -1 );
if found
    %Go through the movie
    [camera_frame,msg]=fread(fileID,6,'uint32','l');
    index=0;
    while msg==6
        if mod(index,2000)==0
            disp(num2str(index));
        end
        %Check for EOF
        if camera_frame(1)~=CAMERA_MOVIE_MAGIC
            offset=ftell(fileID);
            fclose(fileID);
            error(['Wrong magic at offset ',num2str(offset)])
        end
        if camera_frame(2)~=CAMERA_MOVIE_VERSION
            fclose(fileID);
            error(['Unsupported movie version ',num2str(A(2))])
        end
        
        %Return to the beginning of the frame
        fseek(fileID, -24, 0 );
        
        %Read header
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
            
            if msgframe~=40
                offset=ftell(fileID);
                fclose(fileID);
                error(['Corrupted header at offset ',num2str(offset)])
            end
            
            size_x=E(3);
            size_y=E(4);
            
            %Read data
            if (camera_frame(4)==CAMERA_PIXELMODE_MONO_8)
                bpp=1;
                [data,msg]=fread(fileID,[size_x,size_y],'uint8');
                data=uint8(data);
            else
                bpp=2;
                [data,msg]=fread(fileID,[size_x,size_y],'uint16');
                data=uint16(data);
            end
            if msg~=size_x*size_y
                offset=ftell(fileID);
                fclose(fileID);
                error(['Corrupted data at offset ',num2str(offset)])
            end
            
            imagename=[myfolder,'/',myfolder,'_',num2str(index,'%.6d'),'.tiff'];
            try
                image=Tiff(imagename,'w');
            catch ME
                error(['Could not open ',imagename,' for writing.']);
            end
            
            %Set values for image tag
            image.setTag('ImageWidth',size_x);
            image.setTag('ImageLength',size_y);
            image.setTag('Photometric',Tiff.Photometric.MinIsBlack);
            image.setTag('BitsPerSample',8*bpp);
            image.setTag('SamplesPerPixel',1);
            image.setTag('RowsPerStrip',size_y);
            image.setTag('Compression',Tiff.Compression.LZW);
            image.setTag('FillOrder',1);
            image.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
            
            %Write image
            writeEncodedStrip(image,1,data');
            image.close();
            
            index=index+1;
            [camera_frame,msg]=fread(fileID,6,'uint32','l');
        end
        
    end
    
    fclose(fileID);
end



