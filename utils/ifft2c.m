function x = ifft2c(X)
%
%
% x = ifft2c(X)
% 
% orthonormal uncentered 2D ifft for tensor
%
% @yhao

x = sqrt(  size(X,1)*size(X,2)   )*ifft2(X);
