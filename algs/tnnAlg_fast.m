function [X, cost, SNR_iter] = tnnAlg_fast( Stb,StS,X_real,param,verbose )

number_of_iterations = 100;

X = Stb;
res = size(X);
L = zeros(res);

lambda = param.lambda_tnn;
mu = param.mu_tnn;

cost = NaN([1 number_of_iterations+1]);
SNR_iter = NaN([1 number_of_iterations+1]);
if verbose == 1
    diff = StS.*X-Stb;
    cost(1) = 1/2.*norm( diff(:) ).^2+lambda.*calc_tnn(X);
    SNR_iter(1) = 20*log10(norm(X_real(:))/norm(X(:)-X_real(:)));

    % plot/print curves
    figure(11);subplot(121);plot(cost);drawnow;
    figure(11);subplot(122);plot(SNR_iter);drawnow;
    fprintf('Initial----> SNR: %6.4f \n',SNR_iter(1));
end

for iter=1:number_of_iterations
    
    Z = prox_tnn(X+L,lambda/mu);
    
    X_prev = X;
    X = ( Stb+mu.*(Z-L) )./(StS+mu);
    
    L = L+(X-Z);
    
    if verbose == 1
        diff = StS.*X-Stb;
        cost(iter+1) = 1/2.*norm( diff(:) ).^2+lambda.*calc_tnn(X);
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

    