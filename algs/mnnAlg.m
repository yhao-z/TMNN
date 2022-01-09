function [x, cost, SNR_iter] = mnnAlg( Stb,StS,x_real,param,verbose )

number_of_iterations = 100;

x = ifft2c(Stb);
res = size(x);
[n1,n2,n3] = size(x);
L = zeros(res);

lambda = param.lambda_mnn;
mu = param.mu_mnn;

cost = NaN([1 number_of_iterations+1]);
SNR_iter = NaN([1 number_of_iterations+1]);
if verbose==1
    diff = StS.*fft2c(x)-Stb;
    cost(1) = 1/2.*norm( diff(:) ).^2+lambda.*sum(svd(reshape(x,n1*n2,n3),'econ'));
    SNR_iter(1) = 20*log10(norm(x_real(:))/norm(x(:)-x_real(:)));

    % plot/print curves
    figure(11);subplot(121);plot(cost);xlabel('iteration');ylabel('cost');drawnow;
    figure(11);subplot(122);plot(SNR_iter);xlabel('iteration');ylabel('SNR');drawnow;
    fprintf('Initial----> SNR: %6.4f \n',SNR_iter(1));
end

for iter=1:number_of_iterations
    
    Z = prox_nuclear(x+L,lambda/mu);
    
    x_prev = x;
    x = ( Stb+mu.*fft2c(Z-L) )./(StS+mu);
    x = ifft2c(x);
    
    L = L+(x-Z);
    
    if verbose == 1
        diff = StS.*fft2c(x)-Stb;
        cost(iter+1) = 1/2.*norm( diff(:) ).^2+lambda.*sum(svd(reshape(x,n1*n2,n3),'econ'));
        SNR_iter(iter+1) = 20*log10(norm(x_real(:))/norm(x(:)-x_real(:)));

        % plot/print curves
        figure(11);subplot(121);plot(cost);xlabel('iteration');ylabel('cost');drawnow;
        figure(11);subplot(122);plot(SNR_iter);xlabel('iteration');ylabel('SNR');drawnow;
        fprintf('iter: %d----> SNR: %6.4f \n',iter,SNR_iter(iter+1));
    end
    
    % iteration stopping condition
    stepsize_r = norm(x(:)-x_prev(:))/norm(x_prev(:));
    if(stepsize_r<1e-4)
        break;
    end
end
end
    