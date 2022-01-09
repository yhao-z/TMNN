function X = fft2c(x)

% X = fft2c(x)
% 
% orthonormal uncentered 2D fft for tensor
%
%@yhao

X = 1/sqrt(  size(x,1)*size(x,2)   )*fft2(x);

