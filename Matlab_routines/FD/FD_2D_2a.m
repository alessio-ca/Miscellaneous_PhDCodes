%DB (16/03/22): Conversion of FD code to run on 2D videos
%DB (16/03/28): From prelim tests seems FD works better than Hough
%DB (16/03/30): refactoring code 

function grandresult = FD_2D_2a(img,radii,RfConst)

stepR = .3;
maxjj = ceil((radii(2)-radii(1))/stepR);
%RfConst = 10^8;
maxNoise = .1;
if nargin < 4
    maxnumObjects = 50;
    maxRegionsPerRadius = 10;
elseif nargin < 5
    maxnumObjects = vargin{1};
    maxRegionsPerRadius = 10;
else
    maxnumObjects = vargin{1};
    maxRegionsPerRadius = vargin{2};
end
minRadiiNum = 4;
ThreshMaxI = .25;
Thresholding = .75*ThreshMaxI;

nPix = size(img,1);
R = radii(1) + stepR*[0:maxjj];
    
[mg sg] = ndgrid((1:nPix)-0.5,(1:nPix)-0.5);
ipDist = sqrt(((mg-nPix/2)).^2+((sg-nPix/2)).^2);
fii = fftshift(fftn(fftshift(img)));

for jj=1:maxjj
    ip = randn(nPix,nPix);
    ip = abs(ip)*maxNoise*0.01;
    iii = find(ipDist <= R(jj));
    ip(iii) = 100+ip(iii);
    
    fip = fftshift(fftn(fftshift(ip)));
        
    Rf = fip.*conj(fip);
    Rf = Rf./(Rf+RfConst);
        
    Rf = Rf./fip;
    Rf = fii.*Rf;
        
    nii = abs(fftshift(ifftn(fftshift(Rf))));
    
    c = max(nii,[],2);
    c = max(c,[],1);
        
    c(jj) = c;
    l = 1;
    numobjects = maxnumObjects + 1;
    while numobjects > maxnumObjects
        bw = nii > c(jj)*Thresholding*(0.9+l*0.1);
        bnii = bw.*nii;
        numpixels = nnz(bw);
        [L,numobjects] = bwlabel(bw,8);
        L = uint16(L);
        l=l+1;
    end
        
    stats = regionprops(L,'Area','PixelList');
    xcmplot = zeros(numobjects,2);tmpresult = [];
    indexM = 1;
    for l=1:numobjects
        ball = stats(l).PixelList;
        ar = stats(l).Area;
        intensities = bnii(sub2ind(size(bnii),ball(:,2),ball(:,1)));
                    
        if max(intensities) > c*ThreshMaxI
            xcmplot(l,:) = sum(ball.*repmat(intensities,1,2))./ ...
                sum(intensities);
            tmpresult(indexM,:) = [xcmplot(l,1:2)-.5 max(intensities)];
            indexM = indexM + 1;
        end
    end
    %tmpresult holds cm coordinates and max intensities of objects
        
    tmpresult = sortrows(tmpresult,3);
    endresult(1:maxRegionsPerRadius, size(tmpresult,2), jj) = 0;
    if indexM > maxRegionsPerRadius
        endresult(1:maxRegionsPerRadius,:,jj) = tmpresult(indexM-maxRegionsPerRadius:indexM-1,:);
    else
        endresult(1:indexM-1,:,jj) = tmpresult;
    end
end
%endresult now has (x,y,I) and its shape is (numobjects,[x y I],radii)
    
grandresult = zeros(1,4);
l=1;
for jj = 1:maxjj
    for m = 1:maxRegionsPerRadius
        if endresult(m,1,jj)
            grandresult(l,1:2) = endresult(m,1:2,jj); ...
                grandresult(l,3) = R(jj);
            grandresult(l,4) = endresult(m,3,jj);
            l=l+1;
        end
    end
end
sortedX = sortrows(grandresult,1);
totD = size(sortedX,1);
gres = zeros(1,4);
combined = zeros(1,totD);
deleted = zeros(1,totD);
l=1;
for m=1:totD
    if ~deleted(m)
        for n=m+1:totD
            if ~deleted(n)
                dsq = sqrt((sortedX(m,1) - sortedX(n,1))^2 + (sortedX(m,2) - sortedX(n,2))^2);
                temp = max(sortedX(m,3),sortedX(n,3)); 
                if dsq < temp
                    if sortedX(m,4) < sortedX(n,4)
                        deleted(m) = 1;
                        combined(n) = combined(n) + combined(m) + 1;
                        combined(m) = 0;
                        break;
                    else
                        combined(m) = combined(m) + combined(n) + 1;
                        deleted(n) = 1;
                        combined(n) = 0;
                    end
                end
            end
        end
        if combined(m) >= minRadiiNum - 1
            gres(l,:) = sortedX(m,:);
            l = l+1;
        end
    end
end
multS = l-1;
grandresult = sortrows(gres,4);