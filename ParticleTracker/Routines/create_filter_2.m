function [filterFin]=create_filter_2(rMax,rMid,rMin,int1,int2,int3)
% USAGE:    create_filter_2(rMax,rMid,rMin,int1,int2,int3)
% PURPOSE:  Create a filter mask composed of three disks with different
%           relative intensities.
%
%
% ARGUMENT MANAGEMENT:  rMax,rMid and rMin are mandatory and must be positive
%                       
%                       int1,int2 and int3 are optional
%                       Background intensity is always 1
%                       Outer rim is int1+1
%                       Mid rim is 1+int1-int2
%                       Inner rim is 1+int1-int2+int4
%                       
%                       Default is: 1 (outer rim = 2) 
%                                   2 (mid rim = 0) 
%                                   4 (inner rim = 4)
%
%                       If you wish an inverted image (black center and bright corona), simply
%                       type negative intensities (ex: -1 -2 -4 produces an inverted filter of
%                       the default one)
%
% 
if nargin <= 3
   int1 = 1;
   int2 = 2;
   int3 = 4;
elseif nargin >3 && nargin < 6
    error('Error in filter declaration.')
end

%Security control for incorrect user input
r=[rMax,rMid,rMin];
d=diff(r);
if any(r<0) || any(d>0)
    disp('Error in filter declaration.')
    filterFin=0;
    return 
end
if any(mod(r,2) == 1)
    disp('Diameter values must be even.')
    
    filterFin=0;
    return
end
filter_1=1+int1*mat2gray(fspecial('disk',rMax/2));

%If rMid or rMin are set to 0, set correspondent filter to 0
switch rMin
    case 0
        filter_3 = zeros(size(filter_1,1),size(filter_1,2));
    otherwise
        filter_3=int3*mat2gray(padarray(fspecial('disk',rMin/2),[(rMax-rMin)/2 (rMax-rMin)/2],'both'));
end
switch rMid
    case 0
        filter_2 = zeros(size(filter_1,1),size(filter_1,2));
    otherwise
        filter_2=int2*mat2gray(padarray(fspecial('disk',rMid/2),[(rMax-rMid)/2 (rMax-rMid)/2],'both'));
end   

filterFin=mat2gray(filter_1-filter_2+filter_3);


