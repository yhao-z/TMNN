function [x, cost, SNR_iter, times_cost] = tmnnAlg( Stb,StS,x_real,param,verbose )

number_of_iterations = 1000;

x = ifft2c(Stb);
res = size(x);
[n1,n2,n3] = size(x);
L1 = zeros(res);
L2 = zeros(res);

lambda1 = param.lambda1;
lambda2 = param.lambda2;
mu = param.mu;

cost = NaN([1 number_of_iterations+1]);
SNR_iter = NaN([1 number_of_iterations+1]);
if verbose == 1
    diff = StS.*fft2c(x)-Stb;
    cost(1) = 1/2.*norm( diff(:) ).^2+lambda1.*calc_tnn(x)+lambda2.*sum(svd(reshape(x,n1*n2,n3),'econ'));
    SNR_iter(1) = 20*log10(norm(x_real(:))/norm(x(:)-x_real(:)));

    % plot/print curves
    figure(10);subplot(121);plot(cost);xlabel('iteration');ylabel('cost');drawnow;
    figure(10);subplot(122);plot(SNR_iter);xlabel('iteration');ylabel('SNR');drawnow;
    fprintf('Initial----> SNR: %6.4f \n',SNR_iter(1));
end

times_cost = NaN(1,number_of_iterations+1);
times_cost(1) = 0;
tic;
for iter=1:number_of_iterations %&& diff_X>error;
    
    A = prox_tnn(x+L1,lambda1/mu);
    
    B = prox_nuclear(x+L2,lambda2/mu);
    
    x_prev = x;
    x = ( Stb+mu.*fft2c(A+B-L1-L2) )./(StS+2.*mu);
    x = ifft2c(x);
   
    L1 = L1+(x-A);
    L2 = L2+(x-B);
    
    times_cost(iter) = toc;
    
    if verbose == 1
        diff = StS.*fft2c(x)-Stb;
        cost(iter+1) = 1/2.*norm( diff(:) ).^2+lambda1.*calc_tnn(x)+lambda2.*sum(svd(reshape(x,n1*n2,n3),'econ'));
        SNR_iter(iter+1) = 20*log10(norm(x_real(:))/norm(x(:)-x_real(:)));

        % plot/print curves
        figure(10);subplot(121);plot(cost);xlabel('iteration');ylabel('cost');drawnow;
        figure(10);subplot(122);plot(SNR_iter);xlabel('iteration');ylabel('SNR');drawnow;
        fprintf('iter: %d----> SNR: %6.4f \n',iter,SNR_iter(iter+1));
    end
    
    stepsize_r = norm(x(:)-x_prev(:))/norm(x_prev(:));
    if(stepsize_r<1e-4)
        break;
    end
end
end
    