% My first Matlab GUI
% Loads up a volume image and allows you to scroll up and down inside it
% Exper12p48

% Added downsizing from 512x512 to 256x256 in OpenMenuItem_Callback and
% guivolviewih_export_OpeningFcn (IH 7/12/02)
% imoverlay(nslice, on, handles) function to do a live update of the overlay (IH
% 7/12/02)
% Picking up a zsp value when loading IH 7/12/02
% Stop autoloading data at start IH 7/12/02
% Add overlay on loading up droplet file IH 7/12/02
% Changed imoverlay so that it finds current slice IH 8/12/02
% Fit spheres achieved, implemented Clear Points and Fit Sphere button IH 8/12/02
% Added multi-(droplet file) load, sphere deleting, spherelist clearing and non-permanant pointlist display IH 14/12/02
% Also added droplet file save IH 14/12/02
% Added in the spherelist tidying routine IH 15/12/02
% Conversion to run under Matlab 6.1 IH 16/1/03
% Note: a pixel size is hardcoded in the function
% SaveDropletMenuItem_Callback IH 16/1/03
% Added find neighbours function and button in gui IH 20/1/03
% Added in overlay_mode - required update in cirpaint IH 4/3/03

%Callbacks where all the action happends:
%figure1_WindowButtonUpFcn
function varargout = guivolviewih_export(varargin)
% GUIVOLVIEWIH_EXPORT M-file for guivolviewih_export.fig
%      GUIVOLVIEWIH_EXPORT, by itself, creates a new GUIVOLVIEWIH_EXPORT or raises the existing
%      singleton*.
%
%      H = GUIVOLVIEWIH_EXPORT returns the handle to a new GUIVOLVIEWIH_EXPORT or the handle to
%      the existing singleton*.
%
%      GUIVOLVIEWIH_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIVOLVIEWIH_EXPORT.M with the given input arguments.
%
%      GUIVOLVIEWIH_EXPORT('Property','Value',...) creates a new GUIVOLVIEWIH_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before guivolviewih_export_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to guivolviewih_export_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help guivolviewih_export

% Last Modified by GUIDE v2.5 04-Mar-2003 10:04:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @guivolviewih_export_OpeningFcn, ...
                   'gui_OutputFcn',  @guivolviewih_export_OutputFcn, ...
                   'gui_LayoutFcn',  @guivolviewih_export_LayoutFcn, ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before guivolviewih_export is made visible.
function guivolviewih_export_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to guivolviewih_export (see VARARGIN)

% Choose default command line output for guivolviewih_export
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

clear global output
clear global pointlist
clear global spherelist
global output
global spherelist
global pointlist
global alp
global nPix0
global zspp

%spherelist=struct('sx',{},'sy',{},'sz',{},'sr',{},'si',{});
%spherelist=struct('x',{},'y',{},'z',{});
spherelist.sx(1)=0.0;
spherelist.sy(1)=0.0;
spherelist.sz(1)=0.0;
spherelist.sr(1)=0.0;
spherelist.si(1)=0.0;
pointlist.x(1)=0.0;
pointlist.y(1)=0.0;
pointlist.z(1)=0.0;

%zspp = 0.8032;        % 0.93961 (set 6); % 0.93955 (set 7); % (1.405/1.45)*82.50/335/(65.00/256)
zspp = input('Please enter the aspect ratio (nSnell*(boxZ/PixZ)/(boxX/256)) \n');

handles.zsp= zspp; 
guidata(hObject,handles);

% Decide not to autoload an image and droplet file 
% This sets up the initial plot - only do when we are invisible
% so window can get raised using guivolviewih_export.
% global output
% load sphereoverlayih.mat
% %Downsize the image if it is 512x512
% if size(output,1)==512
%     for i=1:size(output,3)
%         B(:,:,i)=imresize(output(:,:,i),[256 256]);
%     end %%for i
%     clear output
%     global output
%     output=B;
%     clear B
% end %end if
% handles.zsp=1.33;
% guidata(hObject,handles);
% 
% %Set the handles right
% nFrames=size(output,3);
% slider_step(1)=1/nFrames;
% slider_step(2)=10/nFrames;
% set(handles.slider2,'sliderstep',slider_step);
% 
% %This loads up the sphere list
% global spherelist
% droplets_file = 'droplets451_688.txt';
% droplets=load(droplets_file);
% spherelist.sx=droplets(:,1);
% spherelist.sy=droplets(:,2);
% spherelist.sz=droplets(:,3);
% %For some reason the radius is in microns, whilst everything else is in
% %pixels....So convert to pixels:
% spherelist.sr=droplets(:,4)/0.3;
% spherelist.si=droplets(:,5);
% 
% %Plot image
% if strcmp(get(hObject,'Visible'),'off')
%     imoverlay(20,1,handles);
% end

% UIWAIT makes guivolviewih_export wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = guivolviewih_export_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global output
global nPix0
global alp
global zspp
[file, path] = uigetfile({'*.mat', 'MAT files (*.mat)';'*.tif','Tif stack files (*.tif)'});
if ~isequal(file, 0)
    % Matlab 6.5 supports an 'index' for uigetfile, but 6.1 doesn't hence the
    % next line
    %[tmp1, tmp2, ext, tmp3]=fileparts(file);
    [tmp1, tmp2, ext]=fileparts(file);
    %clear tmp1,tmp2,tmp3;
    clear tmp1,tmp2;
    file=strcat(path,file);
    if ext=='.mat' %This code loads a MAT file, moderately robust in that it gets the fieldname
        tmp=load(file);
        %The curly brackets in the names{1,1} construct are necessary
        %because its a cell array
        names=fieldnames(tmp);
        output=getfield(tmp,names{1}); % 6.1/6.5 Compatibility - output=tmp.(names{1,1});
        nPix=size(output,1);
        nPix0=nPix;        
        alp= nPix0/256;
        nFrames=size(output,3);
    elseif ext=='.tif' %This code loads a TIFF file
        inf = imfinfo(file);
        nPix = inf(1).Width;
        nPix0=nPix;
        alp= nPix0/256;
        nFrames = size(inf,2); clear inf;
        output=uint8(zeros(nPix, nPix, nFrames));
        for i=1:nFrames
            [output(:,:,i)]=imread(file,i);
        end
    end
% Downsize 512x512 images to 256x256 7/12/02
    if nPix~=256        
        nPix=256;
        for i=1:size(output,3)
            B(:,:,i)=imresize(output(:,:,i),[256 256]);
        end %%for i
        clear output
        global output
        output=B;
        clear B
    end %end if
% Get an aspect ratio 1.33 for older data, 2.16 for newer data
    answer=inputdlg('Input zsp, do not forget downsizing:','Aspect ratio for image file', 1, {num2str(zspp)});
    if ~isempty(answer)
        handles.zsp=str2num(answer{1});
        guidata(hObject,handles);
        %We need to do some stuff here to set the slider up right
        set(handles.slider2,'Max',nFrames);
        slider_step(1)=1/nFrames;
        slider_step(2)=10/nFrames;
        set(handles.slider2,'sliderstep',slider_step);
        set(handles.slider2,'Value',round(nFrames/2));
        set(handles.slicenumber,'String',...
        sprintf('%6s   %4d','Slice:',round(nFrames/2)));
        nslice=round(nFrames/2);
        imoverlay(1,handles);
    end
end

function OpenDropletMenuItem_callback(hObject, eventdata, handles)
%This will load up a user selected droplet file
%Modify to combine successive droplet files together
[file, path] = uigetfile({'*.txt','txt files (*.txt)';'*.mat', 'MAT files (*.mat)'});
global spherelist
global alp
if ~isequal(file, 0)
    % Matlab 6.5 supports an 'index' for uigetfile, but 6.1 doesn't hence the
    % next line
    [tmp1, tmp2, ext]=fileparts(file);
    clear tmp1,tmp2;
    file=strcat(path,file);
    if spherelist.sx(1)==0
        start=1;
    else
        start=length(spherelist.sx)+1;
    end
    syms boxX;
    boxX = input('Please enter the box size in x or y direction in microns. \n');
    if ext=='.mat' %This loads up *.mat file data
        tmp=load(file);
        names=fieldnames(tmp);
        droplets=getfield(tmp,names{1}); % droplets=tmp.(names{1,1}); 6.5 / 6.1 compatibility
        len=length(droplets(:,1))-1;
        spherelist.sx(start:start+len,1)=droplets(:,1);
        spherelist.sy(start:start+len,1)=droplets(:,2);
        spherelist.sz(start:start+len,1)=droplets(:,3);
        %For some reason the radius is in microns, whilst everything else is in
        %pixels....So convert to pixels:
        spherelist.sr(start:start+len,1)=droplets(:,4)/(boxX/256.0); % x or y divided by 256
        spherelist.si(start:start+len,1)=droplets(:,5);
    elseif ext=='.txt' %This loads up *.txt file data        
        droplets=load(file);
        len=length(droplets(:,1))-1;
        spherelist.sx(start:start+len,1)=droplets(:,1)/alp;
        spherelist.sy(start:start+len,1)=droplets(:,2)/alp;
        spherelist.sz(start:start+len,1)=droplets(:,3);
        %For some reason the radius is in microns, whilst everything else is in
        %pixels....So convert to pixels:
        spherelist.sr(start:start+len,1)=droplets(:,4)/(boxX/256.0); % x or y divided by 256
        spherelist.si(start:start+len,1)=droplets(:,5);
    end
    imoverlay(1,handles);
end %~isequal(file,0)
clear spherelist

function SaveDropletMenuItem_Callback(hObject, eventdata, handles)
% This function saves the current spherelist IH 14/12/02
% Converts the radius, but nothing else, to microns.
global spherelist
global alp
[file, path] = uiputfile({'*.txt','txt files (*.txt)'});
syms boxX;
    boxX = input('Please enter the box size in x or y direction in microns. \n');
if ~isequal(file,0)    
    file=strcat(path,file);
    droplets(:,1)=spherelist.sx*alp;
    droplets(:,2)=spherelist.sy*alp;
    droplets(:,3)=spherelist.sz;
    droplets(:,4)=spherelist.sr*(boxX/256.0);
    droplets(:,5)=spherelist.si;
    droplets(droplets(:,1)==0,:) = [];
    save(file,'droplets','-ascii')
end

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

% 

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global output;
nslice=round(get(handles.slider2,'Value'));
imoverlay(1,handles);

set(handles.slicenumber,'String',...
    sprintf('%6s   %4d','Slice:',nslice));

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This is the callback that does all the work!

%Extract position of cursor and image
thing1=get(hObject,'CurrentPoint');
thing2=get(handles.axes1,'Position');
xmo=thing1(1);
ymo=thing1(2);
xorg=thing2(1);
yorg=thing2(2);
width=thing2(3);
height=thing2(4);

%Calculate pixel position in terms of image pixels
%Should fiddle this to convert to micron positions
%And to handle images of different sizes

global output
global pointlist

xsize=size(output,1);
ysize=size(output,2);
x=((xmo-xorg)/width)*xsize;
y=((height-(ymo-yorg))/height)*ysize;

%Don't act on mouse clicks outside the image location
%This is hard coded for 256x256 images, oh not it isn't
%This size refers to the canvas not the image


if (xmo-xorg)>width | (xmo-xorg)<1 
    return;
end

if (height-(ymo-yorg))>height | (height-(ymo-yorg))<1
    return;
end

%Get the z position from the slider
z=round(get(handles.slider2,'Value'));

%Set coordinate positions
set(handles.cursorx,'String',...
    sprintf('%6s   %4d','x:    ',x));
set(handles.cursory,'String',...
    sprintf('%6s   %4d','y:    ',y));
set(handles.cursorz,'String',...
    sprintf('%6s   %4d','z:    ',z));
%

if pointlist.x(1)==0 & pointlist.y(1)==0 & pointlist.z(1)==0
    pointlist.x(1)=x;
    pointlist.y(1)=y;
    pointlist.z(1)=z;
else
    i=length(pointlist.x);
    pointlist.x(i+1)=x;
    pointlist.y(i+1)=y;
    pointlist.z(i+1)=z;
end
%Now find out where the centre of the selected sphere was
%This is all done in pixels, the final number in the call is the xz aspect
%ratio. Behaves when no sphere is
%selected.
global spherelist;
index=cirwhich(spherelist.sx, spherelist.sy,spherelist.sz,spherelist.sr,...
                x,y,z,handles.zsp);
%Takes the first sphere if two spheres overlap
if length(index)>1
    disp('index    x      y      z     r     i')
    for j=1:length(index)
        k=index(j);
        disp(sprintf('%5d %6.2f %6.2f %6.2f %6.2f %6.4f',k,spherelist.sx(k),spherelist.sy(k),spherelist.sz(k),spherelist.sr(k),spherelist.si(k)))    
    end
end
index=index(1);
if index~=0
    sx=spherelist.sx(index);
    sy=spherelist.sy(index);
    sz=spherelist.sz(index);
    sr=spherelist.sr(index);
    si=spherelist.si(index)/max(max(spherelist.si));
    set(handles.spi,'String',...
        sprintf('%6s   %d','Index: ',index));
    set(handles.spx,'String',...
        sprintf('%6s   %6.1f','x:    ',sx));
    set(handles.spy,'String',...
        sprintf('%6s   %6.1f','y:    ',sy));
    set(handles.spz,'String',...
        sprintf('%6s   %6.1f','z:    ',sz));
    set(handles.spr,'String',...
        sprintf('%6s   %6.2f','r:    ',sr));
    set(handles.spint,'String',...
        sprintf('%6s   %6.4f','Intensity:    ',si));
%This is optional code to fill-in a selected sphere
%Need to reverse the x,y coordinates again!
    %output=cirpaint2(sy,sx,sz,sr,output,1,1.33,1);
    %imshow(output(:,:,z));
else %This executes if no sphere is found
    set(handles.spi,'String',...
        'No sphere found');
    set(handles.spx,'String',...
        ' ');
    set(handles.spy,'String',...
        ' ');
    set(handles.spz,'String',...
        ' ');
    set(handles.spr,'String',...
        ' ');
    set(handles.spint,'String',...
        ' ');
end

%This goes and gets the image data and puts a marker in
%Need to swap the coordinates around
%output(y-1:y+1,x-1:x+1,z)=255;
imoverlay(1,handles);

clear spherelist;
clear output;
clear pointlist;

% --- Executes on mouse press over figure background.
function figure1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


function imoverlay(on,handles)
% Function to plot an image with the sphere overlay
% IH 7/12/02
global output
global spherelist
global pointlist

nslice=round(get(handles.slider2,'Value'));
mode=get(handles.overlay_mode,'Value')-1;
zsp=handles.zsp; %2.16 for 4nov02.tif stack and 1.33 for previous data
if on
    if spherelist.sr(1)~=0
        %Temporary code to put outlines in black replace final flag with 0
        %to get back white outlines
        overlayed=cirpaint2(spherelist.sy,spherelist.sx,spherelist.sz,spherelist.sr,output(:,:,nslice),nslice,zsp,mode);
    else
        overlayed=output(:,:,nslice);
    end
    if pointlist.x(1)~=0
        for i=1:length(pointlist.x)
            %disp('index    x    y     z')
            if pointlist.z(i)==nslice
                %disp(sprintf('%d     %d     %d     %d',i,pointlist.x(i),pointlist.y(i),pointlist.z(i)))    
                overlayed(pointlist.y(i)-1:pointlist.y(i)+1,pointlist.x(i)-1:pointlist.x(i)+1)=255;
            end
        end
    end
else
    overlayed=output(:,:,nslice);
end
imshow(overlayed);



% --- Executes on button press in clear_pointlist.
function clear_pointlist_Callback(hObject, eventdata, handles)
% hObject    handle to clear_pointlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This reinitialises the pointlist IH 8/12/02
clear global pointlist;
global pointlist;
% pointlist=struct('x',{},'y',{},'z',{});
pointlist.x(1)=0.0;
pointlist.y(1)=0.0;
pointlist.z(1)=0.0;
imoverlay(1,handles);
clear pointlist;


% --- Executes on button press in fit_sphere.
function fit_sphere_Callback(hObject, eventdata, handles)
% hObject    handle to fit_sphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global pointlist;
global spherelist;

[x,y,z,r]=spherefit(pointlist.x,pointlist.y,pointlist.z*handles.zsp);

if spherelist.sr(1)==0
    i=0
else
    i=length(spherelist.sx);
end
spherelist.sx(i+1,1)=x;
spherelist.sy(i+1,1)=y;
spherelist.sz(i+1,1)=z/handles.zsp;
spherelist.sr(i+1,1)=r;
spherelist.si(i+1,1)=1;
disp('index    x      y      z     r')
disp(sprintf('%5d %6.2f %6.2f %6.2f %6.2f (%6.2f)',i+1,x,y,z,r,r/0.3))

%Replot overlay
imoverlay(1,handles);

clear pointlist;
clear spherelist;

% --- Executes on button press in delete_sphere.
function delete_sphere_Callback(hObject, eventdata, handles)
% hObject    handle to delete_sphere (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Added IH 14/12/02
answer=inputdlg('Index of sphere to delete:','Delete sphere', 1, {'1'});
if ~isempty(answer)
    index=round(str2num(answer{1}));
    global spherelist
    if index>0 & index<=length(spherelist.sx)
        spherelist.sx(index)=0;
        spherelist.sy(index)=0;
        spherelist.sz(index)=0;
        spherelist.sr(index)=0;
        spherelist.si(index)=0;
    end
    imoverlay(1,handles);
end %~isempty

% --- Executes on button press in clear_spherelist.
function clear_spherelist_Callback(hObject, eventdata, handles)
% hObject    handle to clear_spherelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This is the clear spherelist callback
clear global spherelist
global spherelist
spherelist.sx(1)=0.0;
spherelist.sy(1)=0.0;
spherelist.sz(1)=0.0;
spherelist.sr(1)=0.0;
spherelist.si(1)=0.0;
imoverlay(1,handles);
clear spherelist


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in tidy_spherelist.
function tidy_spherelist_Callback(hObject, eventdata, handles)
% hObject    handle to tidy_spherelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% The spherelist tidying function
% IH 15/12/02

%Only do any of this if the spherelist is not empty
global spherelist
disp('=========')
disp(sprintf('Initial spherelist length: %4d',length(spherelist.sr)))
if ~(spherelist.sr(1)==0 & length(spherelist.sr)==1)    
    %First average/merge spheres with close centres
    for i=1:length(spherelist.sx)
        %First of all we pop the current sphere centre into temporary
        %variables
        thisx=spherelist.sx(i);
        thisy=spherelist.sy(i);
        thisz=spherelist.sz(i);
        % Found which spheres the sphere centre falls within
        index=cirwhich(spherelist.sx, spherelist.sy,spherelist.sz,spherelist.sr,...
                thisx,thisy,thisz,handles.zsp);
        % Every sphere centre should fall in at least one sphere (its own)
        %If index is only 1 in length then do nothing
        if length(index)>1
            disp('Trimming spheres')
            disp('index    x      y      z     r     i')
            for m=1:length(index)
                k=index(m);
                disp(sprintf('%5d %6.2f %6.2f %6.2f %6.2f %6.4f',k,...
                    spherelist.sx(k),spherelist.sy(k),spherelist.sz(k),spherelist.sr(k),spherelist.si(k)))    
            end
            %Zero the sums
            sumx=0;sumy=0;sumz=0;sumr=0;sumi=0;
            %Loop over the list of found spheres
            for j=1:length(index)
                %Sum the sphere parameters
                %At this point
                sumx=sumx+spherelist.sx(index(j));
                sumy=sumy+spherelist.sy(index(j));
                sumz=sumz+spherelist.sz(index(j));
                sumr=sumr+spherelist.sr(index(j));
                sumi=sumi+spherelist.si(index(j));
                %Zero all but the first sphere in the index list
                if j>1
                    spherelist.sx(index(j))=0.0;
                    spherelist.sy(index(j))=0.0;
                    spherelist.sz(index(j))=0.0;
                    spherelist.sr(index(j))=0.0;
                    spherelist.si(index(j))=0.0;
                end %j>1
            end %for j
            %Save the averaged sphere back into the original slot
            spherelist.sx(i)=sumx/length(index);
            spherelist.sy(i)=sumy/length(index);
            spherelist.sz(i)=sumz/length(index);
            spherelist.sr(i)=sumr/length(index);
            spherelist.si(i)=sumi/length(index);
        end % length
    end % for i
%Then remove entries with sr=0
%First duplicate the spherelist
    temp=spherelist;
%Then clear it, nicely
    clear global spherelist
    global spherelist
    spherelist.sx(1)=0.0;
    spherelist.sy(1)=0.0;
    spherelist.sz(1)=0.0;
    spherelist.sr(1)=0.0;
    spherelist.si(1)=0.0;
%Now cycle through temp and only copy over non-zero elements
    j=1;
    for i=1:length(temp.sx)
        if ~(temp.sr(i)==0)
            spherelist.sx(j,1)=temp.sx(i);
            spherelist.sy(j,1)=temp.sy(i);
            spherelist.sz(j,1)=temp.sz(i);
            spherelist.sr(j,1)=temp.sr(i);
            spherelist.si(j,1)=temp.si(i);
            j=j+1;
        end 
    end % length(temp.sx)
end %empty check
disp(sprintf('Final spherelist length: %4d',length(spherelist.sr)))
imoverlay(1,handles);
clear spherelist


% --- Executes on button press in pushbutton7.
function find_neighbours_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% The find nieghbours function
% IH 20/1/03

global spherelist
answer=inputdlg('Index of parent sphere:','Find Neighbours', 1, {'1'});
if ~isempty(answer)
    org=round(str2num(answer{1}));
    xo=spherelist.sx(org);
    yo=spherelist.sy(org);
    zo=spherelist.sz(org);
    ro=spherelist.sr(org);
    j=1;
    neighbours(j)=org;
    j=j+1;
    disp('=========')
    disp('Origin sphere');
    disp('index    x      y      z     r')
    disp(sprintf('%5d %6.2f %6.2f %6.2f %6.2f',org,xo,yo,zo,ro))
    disp('Neighbour spheres');
    disp('index    x      y      z     r')
    for i=1:length(spherelist.sx)
        thisx=spherelist.sx(i);
        thisy=spherelist.sy(i);
        thisz=spherelist.sz(i);
        thisr=spherelist.sr(i);
        dist=sqrt((xo-thisx)^2+(yo-thisy)^2+(handles.zsp*(zo-thisz))^2);
        if dist<(ro+thisr) & i~=org
            neighbours(j)=i;
            j=j+1;
            disp(sprintf('%5d %6.2f %6.2f %6.2f %6.2f',i,thisx,thisy,thisz,thisr))
        end %dist<(ro+thisr)
    end % for i
    disp(sprintf('Co-ordination number: %5d', length(neighbours)-1));
end % ~isempty


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function overlay_mode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to overlay_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in overlay_mode.
function overlay_mode_Callback(hObject, eventdata, handles)
% hObject    handle to overlay_mode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns overlay_mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from overlay_mode
imoverlay(1,handles);



% --- Creates and returns a handle to the GUI figure. 
function h1 = guivolviewih_export_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end

h1 = figure(...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',[0.831372549019608 0.815686274509804 0.784313725490196],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','guivolviewih',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[644 292 492 415],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'WindowButtonMotionFcn','guivolviewih_export(''figure1_WindowButtonMotionFcn'',gcbo,[],guidata(gcbo))',...
'WindowButtonUpFcn','guivolviewih_export(''figure1_WindowButtonUpFcn'',gcbo,[],guidata(gcbo))',...
'ButtonDownFcn','guivolviewih_export(''figure1_ButtonDownFcn'',gcbo,[],guidata(gcbo))',...
'HandleVisibility','callback',...
'Tag','figure1',...
'UserData',zeros(1,0));

setappdata(h1, 'GUIDEOptions', struct(...
'active_h', 1.340007e+002, ...
'taginfo', struct(...
'figure', 2, ...
'axes', 8, ...
'pushbutton', 8, ...
'popupmenu', 2, ...
'slider', 6, ...
'text', 19, ...
'edit', 11, ...
'listbox', 3), ...
'override', 0, ...
'release', 13, ...
'resize', 'none', ...
'accessibility', 'callback', ...
'mfile', 1, ...
'callbacks', 1, ...
'singleton', 1, ...
'syscolorfig', 1, ...
'lastSavedFile', 'C:\Work\UMIST\Matlab\guivolviewih.m'));


h2 = axes(...
'Parent',h1,...
'Units','pixels',...
'ALim',get(0,'defaultaxesALim'),...
'ALimMode','manual',...
'CameraPosition',[0.5 0.5 9.16025403784439],...
'CameraPositionMode','manual',...
'CameraTarget',[0.5 0.5 0.5],...
'CameraTargetMode','manual',...
'CameraUpVector',[0 1 0],...
'CameraUpVectorMode','manual',...
'CameraViewAngle',6.60861036031192,...
'CameraViewAngleMode','manual',...
'CLim',get(0,'defaultaxesCLim'),...
'CLimMode','manual',...
'Color',get(0,'defaultaxesColor'),...
'ColorOrder',get(0,'defaultaxesColorOrder'),...
'DataAspectRatio',get(0,'defaultaxesDataAspectRatio'),...
'DataAspectRatioMode','manual',...
'PlotBoxAspectRatio',get(0,'defaultaxesPlotBoxAspectRatio'),...
'PlotBoxAspectRatioMode','manual',...
'Position',[50 136 256 256],...
'TickDir',get(0,'defaultaxesTickDir'),...
'TickDirMode','manual',...
'XColor',get(0,'defaultaxesXColor'),...
'XLim',get(0,'defaultaxesXLim'),...
'XLimMode','manual',...
'XTick',[0 0.2 0.4 0.6 0.8 1],...
'XTickLabel',{ '0  ' '0.2' '0.4' '0.6' '0.8' '1  ' },...
'XTickLabelMode','manual',...
'XTickMode','manual',...
'YColor',get(0,'defaultaxesYColor'),...
'YLim',get(0,'defaultaxesYLim'),...
'YLimMode','manual',...
'YTick',[0 0.2 0.4 0.6 0.8 1],...
'YTickLabel',{ '0  ' '0.2' '0.4' '0.6' '0.8' '1  ' },...
'YTickLabelMode','manual',...
'YTickMode','manual',...
'ZColor',get(0,'defaultaxesZColor'),...
'ZLim',get(0,'defaultaxesZLim'),...
'ZLimMode','manual',...
'ZTick',[0 0.5 1],...
'ZTickLabel','',...
'ZTickLabelMode','manual',...
'ZTickMode','manual',...
'ButtonDownFcn','guivolviewih_export(''axes1_ButtonDownFcn'',gcbo,[],guidata(gcbo))',...
'CreateFcn','guivolviewih_export(''axes1_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','axes1',...
'UserData',zeros(1,0));


h3 = get(h2,'title');

set(h3,...
'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[0.498046875 1.025390625 1.00005459937205],...
'VerticalAlignment','bottom',...
'HandleVisibility','off');

h4 = get(h2,'xlabel');

set(h4,...
'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[0.498046875 -0.091796875 1.00005459937205],...
'VerticalAlignment','cap',...
'HandleVisibility','off');

h5 = get(h2,'ylabel');

set(h5,...
'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','center',...
'Position',[-0.111328125 0.494140625 1.00005459937205],...
'Rotation',90,...
'VerticalAlignment','bottom',...
'HandleVisibility','off');

h6 = get(h2,'zlabel');

set(h6,...
'Parent',h2,...
'Color',[0 0 0],...
'HorizontalAlignment','right',...
'Position',[-0.193359375 1.087890625 1.00005459937205],...
'HandleVisibility','off',...
'Visible','off');

h7 = uimenu(...
'Parent',h1,...
'Callback','guivolviewih_export(''FileMenu_Callback'',gcbo,[],guidata(gcbo))',...
'Label','File',...
'Tag','FileMenu');

h8 = uimenu(...
'Parent',h7,...
'Callback','guivolviewih_export(''OpenMenuItem_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Open Image...',...
'Tag','OpenMenuItem');

h9 = uimenu(...
'Parent',h7,...
'Callback','guivolviewih_export(''OpenDropletMenuItem_callback'',gcbo,[],guidata(gcbo))',...
'Label','Open Droplet file...',...
'Tag','OpenDropletMenuItem');

h10 = uimenu(...
'Parent',h7,...
'Callback','guivolviewih_export(''SaveDropletMenuItem_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Save Droplet file...',...
'Tag','SaveDropletMenuItem');

h11 = uimenu(...
'Parent',h7,...
'Callback','guivolviewih_export(''CloseMenuItem_Callback'',gcbo,[],guidata(gcbo))',...
'Label','Close',...
'Separator','on',...
'Tag','CloseMenuItem');

h12 = uicontrol(...
'Parent',h1,...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback','guivolviewih_export(''slider2_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Max',59,...
'Min',1,...
'Position',[331 135 20 256],...
'String',{ '' },...
'Style','slider',...
'SliderStep',[0.017 0.085],...
'Value',20,...
'ButtonDownFcn','guivolviewih_export(''slider2_ButtonDownFcn'',gcbo,[],guidata(gcbo))',...
'CreateFcn','guivolviewih_export(''slider2_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','slider2');


h13 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 29.5384615384615 15.8 1.38461538461538],...
'String','Slice:',...
'Style','text',...
'Tag','slicenumber');


h14 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 28.2307692307692 21.4 1.15384615384615],...
'String','Cursor Position',...
'Style','text',...
'Tag','text2');


h15 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 26.3846153846154 15.8 1.46153846153846],...
'String','x:',...
'Style','text',...
'Tag','cursorx');


h16 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 24.6923076923077 15.8 1.61538461538462],...
'String','y:',...
'Style','text',...
'Tag','cursory');


h17 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 23.0769230769231 15.8 1.69230769230769],...
'String','z:',...
'Style','text',...
'Tag','cursorz');


h18 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 22 19.2 1.23076923076923],...
'String','Sphere parameters',...
'Style','text',...
'Tag','text7');


h19 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 18.6923076923077 21.2 1.46153846153846],...
'String','x:',...
'Style','text',...
'Tag','spx');


h20 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 17.1538461538462 21.4 1.46153846153846],...
'String','y:',...
'Style','text',...
'Tag','spy');


h21 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 15.3846153846154 20.6 1.69230769230769],...
'String','z:',...
'Style','text',...
'Tag','spz');


h22 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 20.2307692307692 21.4 1.46153846153846],...
'String','Index:',...
'Style','text',...
'Tag','spi');


h23 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 14 20.4 1.53846153846154],...
'String','r:',...
'Style','text',...
'Tag','spr');


h24 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'HorizontalAlignment','left',...
'ListboxTop',0,...
'Position',[64.2 12.5384615384615 20.2 1.15384615384615],...
'String','Intensity:',...
'Style','text',...
'Tag','spint');


h25 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','guivolviewih_export(''clear_pointlist_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[9.8 7 15.8 1.76923076923077],...
'String','Clear Pointlist',...
'Tag','clear_pointlist');


h26 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','guivolviewih_export(''fit_sphere_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[9.8 4.61538461538462 16 1.76923076923077],...
'String','Fit sphere',...
'Tag','fit_sphere');


h27 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','guivolviewih_export(''delete_sphere_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[9.8 2.23076923076923 16 1.76923076923077],...
'String','Delete Sphere',...
'Tag','delete_sphere');


h28 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','guivolviewih_export(''clear_spherelist_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[27.2 7 15.8 1.76923076923077],...
'String','Clear Spherelist',...
'Tag','clear_spherelist');


h29 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','guivolviewih_export(''tidy_spherelist_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[27.2 4.61538461538462 15.8 1.76923076923077],...
'String','Tidy Spherelist',...
'Tag','tidy_spherelist');


h30 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','guivolviewih_export(''find_neighbours_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[27.2 2.23076923076923 15.8 1.76923076923077],...
'String','Find neighbours',...
'Tag','pushbutton7');


h31 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'ListboxTop',0,...
'Position',[72.8 10.4615384615385 17.2 1.38461538461538],...
'String','Overlay Mode:',...
'Style','text',...
'Tag','text18');


h32 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','guivolviewih_export(''overlay_mode_Callback'',gcbo,[],guidata(gcbo))',...
'Max',0,...
'Min',1,...
'Position',[64.2 4.84615384615385 18 5.46153846153846],...
'String',{ 'White outline' 'White filled' 'Black outline' 'Black filled' 'None' },...
'Style','listbox',...
'Value',1,...
'CreateFcn','guivolviewih_export(''overlay_mode_CreateFcn'',gcbo,[],guidata(gcbo))',...
'Tag','overlay_mode');



hsingleton = h1;


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      GUIVOLVIEWIH_EXPORT, by itself, creates a new GUIVOLVIEWIH_EXPORT or raises the existing
%      singleton*.
%
%      H = GUIVOLVIEWIH_EXPORT returns the handle to a new GUIVOLVIEWIH_EXPORT or the handle to
%      the existing singleton*.
%
%      GUIVOLVIEWIH_EXPORT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUIVOLVIEWIH_EXPORT.M with the given input arguments.
%
%      GUIVOLVIEWIH_EXPORT('Property','Value',...) creates a new GUIVOLVIEWIH_EXPORT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.4 $ $Date: 2002/05/31 21:44:31 $

gui_StateFields =  {'gui_Name'
                    'gui_Singleton'
                    'gui_OpeningFcn'
                    'gui_OutputFcn'
                    'gui_LayoutFcn'
                    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);        
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [getfield(gui_State, gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % GUIVOLVIEWIH_EXPORT
    % create the GUI
    gui_Create = 1;
elseif numargin > 3 & ischar(varargin{1}) & ishandle(varargin{2})
    % GUIVOLVIEWIH_EXPORT('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % GUIVOLVIEWIH_EXPORT(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = 1;
end

if gui_Create == 0
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else
        feval(varargin{:});
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.
    
    % Do feval on layout code in m-file if it exists
    if ~isempty(gui_State.gui_LayoutFcn)
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        end
    end
    
    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    
    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig 
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA
        guidata(gui_hFigure, guihandles(gui_hFigure));
    end
    
    % If user specified 'Visible','off' in p/v pairs, don't make the figure
    % visible.
    gui_MakeVisible = 1;
    for ind=1:2:length(varargin)
        if length(varargin) == ind
            break;
        end
        len1 = min(length('visible'),length(varargin{ind}));
        len2 = min(length('off'),length(varargin{ind+1}));
        if ischar(varargin{ind}) & ischar(varargin{ind+1}) & ...
                strncmpi(varargin{ind},'visible',len1) & len2 > 1
            if strncmpi(varargin{ind+1},'off',len2)
                gui_MakeVisible = 0;
            elseif strncmpi(varargin{ind+1},'on',len2)
                gui_MakeVisible = 1;
            end
        end
    end
    
    % Check for figure param value pairs
    for index=1:2:length(varargin)
        if length(varargin) == index
            break;
        end
        try, set(gui_hFigure, varargin{index}, varargin{index+1}), catch, break, end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end
    
    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});
    
    if ishandle(gui_hFigure)
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
        
        % Make figure visible
        if gui_MakeVisible
            set(gui_hFigure, 'Visible', 'on')
            if gui_Options.singleton 
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        rmappdata(gui_hFigure,'InGUIInitialization');
    end
    
    % If handle visibility is set to 'callback', turn it on until finished with
    % OutputFcn
    if ishandle(gui_hFigure)
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end
    
    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end
    
    if ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end    

function gui_hFigure = local_openfig(name, singleton)
if nargin('openfig') == 3 
    gui_hFigure = openfig(name, singleton, 'auto');
else
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
end

