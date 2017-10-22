function diami = DropletDetector(img,lowDiam,highDiam,varargin)

%Optional:
%   - Low pass diam (default = 10)
%   - Canny param (Threshold, Gamma) (default = [0.4,1.])
%   - dilation diam (default = 1)
%   - radius correction factor (default = 1.075)


lpDiam = 10;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'LPDiam')
        lpDiam = varargin{n+1};
    end
end

CannyParam = [0.4,1.];
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Canny')
        CannyParam = varargin{n+1};
    end
end
CannyThr = CannyParam(1);
CannyGamma = CannyParam(2);

dilDiam = 1;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Dil')
        dilDiam = varargin{n+1};
    end
end

radCorr = 1.075;
for n = 1:2:length(varargin)
    if strcmpi(varargin{n},'Rad')
        radCorr = varargin{n+1};
    end
end


%Contrast enhancement
imgJ=adapthisteq(img);

%Low-pass filtering and high-pass filtering
h=fspecial('disk',lpDiam);
imgJ2=imfilter(imgJ,h);
imgJ3 = imsubtract(imgJ2,imgJ);
imgJ4=adapthisteq(imgJ3);

%Canny-edge detection & dilation
bw=edge(imgJ4,'canny',CannyThr,CannyGamma);
se=strel('disk',dilDiam);
bw2=imdilate(bw,se);


bw2=not(bw2);
bw3=ones(size(bw2));

%Connecting components & filling holes
obj=bwconncomp(bw2);
LM = labelmatrix(obj);
stats = regionprops(LM, 'MajorAxisLength','MinorAxisLength','Area');
[~,idx]=find((lowDiam<=[stats.MinorAxisLength])& ([stats.MajorAxisLength]<=highDiam));
bw3 = ismember(LM,idx);
bw4=imfill(bw3,'holes');

%Fit with circles
obj=bwconncomp(bw4);
props = regionprops(obj,'ConvexHull');
ND_temp = obj.NumObjects;
diami=zeros(ND_temp);

for j=1:ND_temp
    if(size(props(j).ConvexHull,1)<=20) %Check minimum numbers of points for a circle fit
        continue
    end
    x=props(j).ConvexHull(:,1);
    y=props(j).ConvexHull(:,2);
    k=convhull(x,y);
    XY=[x(k),y(k)];
    [~,RC]=fitcircle(XY);
    RC = radCorr*RC;
    if(2*RC < lowDiam || 2*RC > highDiam || ... %Check on size
            any(XY(:)<10) || ... %Check on low border
            any(XY(:,1) > size(bw4,2)-10) || ... %Check on high border, coord1
            any(XY(:,2) > size(bw4,1)-10))      %Check on high border, coord2
        continue
    end
    diami(j)=2*RC;    
end

diami = diami(diami~=0);




    
        

            










