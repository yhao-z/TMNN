function snr = SNR(Xfull,Xrecover)

MSE = norm(Xfull(:)-Xrecover(:))^2;
snr = 10*log10(norm(Xfull(:))^2/MSE);