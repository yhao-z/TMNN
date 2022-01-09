function tnn = calc_tnn(X)
X = fft(X,[],3);
tnn = 0;
for i = 1:size(X,3)
    s = svd(X(:,:,i),'econ');
    tnn = tnn+sum(s);
end
end