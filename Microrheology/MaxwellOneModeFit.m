function F = MaxwellOneModeFit(c,xdata,ydata)
% c is a two-element vector whose components are the Maxwell parameters n
% and G respectively.
% xdata is the angular frequency (type: real vector)
% ydata is the complex modulus G* (type: complex vector)

errors = ((c(1)*xdata).^2/c(2))./(1+(c(1).*xdata/c(2)).^2) ...
    + 1i * (c(1)*xdata)./(1+(c(1).*xdata/c(2)).^2)...
    - ydata;
F = [real(errors(:)),imag(errors(:))];

end