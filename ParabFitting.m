function [out]=ParabFitting(A)
% 
% USAGE:    ParabFitting (A)
% PURPOSE:  Fit a 2D matrix of values to an ellipsoidal paraboloid centered in the matrix center by independent
%           parametrization of x and y direction to a 2nd order polynomial
%
% CREATED: Alessio Caciagli, University of Cambridge, 05/12/2015
out=zeros(size(A,3),2);
for i=1:size(A,3)
    x=(1:length(A(:,:,i)))'; 
    y=A(ceil(length(A(:,:,i))/2),x)'; %sweep on x axis
    fit=[ones(length(A(:,:,i)),1),x,x.^2]; %2nd order polynomial fit
    a=(fit\y)'; %Fitting parameters
    x_fit=(1:0.01:length(A(:,:,i)));
    y_fit=a(1) + a(2).*x_fit + a(3).*x_fit.^2;
    [~,xm]=max(y_fit); %x-position of max

    y=A(x,ceil(length(A(:,:,i))/2)); %sweep on x axis
    a=(fit\y)'; %Fitting parameters
    y_fit=a(1) + a(2).*x_fit + a(3).*x_fit.^2;
    [~,ym]=max(y_fit); %y-position of max
    out(i,:)=[x_fit(xm),x_fit(ym)];
end


