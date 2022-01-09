function [X, cost, SNR_iter, times_cost] = tmnnAlg_fast( Stb,StS,X_real,param,verbose )

number_of_iterations = 1000;
lambda1 = param.lambda1;
lambda2 = param.lambda2;
mu = param.mu;
X = Stb;
res = size(X);
[n1,n2,n3] = size(X);
L1 = zeros(res);
L2 = zeros(res);

cost = NaN([1 number_of_iterations+1]);
SNR_iter = NaN([1 number_of_iterations+1]);
if verbose == 1
    diff = StS.*X-Stb;
    cost(1) = 1/2.*norm( diff(:) ).^2+lambda1.*calc_tnn(X)+lambda2.*sum(svd(reshape(X,n1*n2,n3),'econ'));
    SNR_iter(1) = 20*log10(norm(X_real(:))/norm(X(:)-X_real(:)));

    % plot/print curves
    figure(11);subplot(121);plot(cost);xlabel('iteration');ylabel('cost');drawnow;
    figure(11);subplot(122);plot(SNR_iter);xlabel('iteration');ylabel('SNR');drawnow;
    fprintf('Initial----> SNR: %6.4f \n',SNR_iter(1));
end

times_cost = NaN(1,number_of_iterations+1);
times_cost(1) = 0;
tic;
for iter=1:number_of_iterations
    
    A = prox_tnn(X+L1,lambda1/mu);
    
    B = prox_nuclear(X+L2,lambda2/mu);
    
    X_prev = X;
    X = ( Stb+mu.*(A+B-L1-L2) )./(StS+2.*mu);
   
    L1 = L1+(X-A);
    L2 = L2+(X-B);
    
    times_cost(iter) = toc;
    
    if verbose == 1
        diff = StS.*X-Stb;
        cost(iter+1) = 1/2.*norm( diff(:) ).^2+lambda1.*calc_tnn(X)+lambda2.*sum(svd(reshape(X,n1*n2,n3),'econ'));
        SNR_iter(iter+1) = 20*log10(norm(X_real(:))/norm(X(:)-X_real(:)));
    
        % plot/print curves
        figure(11);subplot(121);plot(cost);xlabel('iteration');ylabel('cost');drawnow;
        figure(11);subplot(122);plot(SNR_iter);xlabel('iteration');ylabel('SNR');drawnow;
        fprintf('iter: %d----> SNR: %6.4f \n',iter,SNR_iter(iter+1));
    end
    
    % iteration stopping condition
    stepsize_r = norm(X(:)-X_prev(:))/norm(X_prev(:));
    if(stepsize_r<1e-4)
        break;
    end
end
end
    