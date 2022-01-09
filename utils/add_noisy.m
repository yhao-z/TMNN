function b_noised = add_noisy(b, SNR)
% add complex gaussian noisy to the k-space sampling b in MRI
% reconstruction. Notes that, b can be a vector or a matrix.
%
% write by yinghao ZHANG, HIT

b_noised = b;
ind = find(b~=0);
b_vec = b(ind);
stdsignal = sqrt(sum(abs(b_vec).^2)/length(b_vec));
stdnoise = 0.1^(SNR/20)*stdsignal;
b_vec = b_vec + stdnoise*(randn(size(b_vec)) + 1i*randn(size(b_vec)));%*sqrt(0.5);
b_noised(ind) = b_vec;
end