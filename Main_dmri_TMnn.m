%% Initialize:
clear;clc;
addpath(genpath(pwd));
%% 4 datasets
% --------------------------------------------------
% load('E:\Yhao\DATA\dmri\aperiodic_pincat.mat','new');
% X = fft2(new);
% --------------------------------------------------
% load('E:\Yhao\DATA\dmri\invivo_perfusion.mat')
% X = fft2c(x);
% --------------------------------------------------
load('E:\Yhao\DATA\dmri\data_tmi.mat','seq');
X = fft2c(seq);
% --------------------------------------------------
% load('E:\Yhao\DATA\dmri\breathing.mat','Data_xy_t');
% X = fft2(Data_xy_t);

[n1,n2,n3] = size(X);
%% normalize
maxX = max(abs(X(:)));
X = X./maxX;
%% 4 sampling mask
uds_ratio_or_lines =16;
disp('*************************************');
disp(['uds_ratio = ',num2str(uds_ratio_or_lines)]);
disp('*************************************');

% ------------------variable density random 2d sampling------------------
% sampling_mask = genrate_binary_sampling_map(n1,n2,uds_ratio_or_lines,n3); 

% ------------------variable density randome x sampling------------------
% sampling_mask = genrate_ylines_sampling_map(n1,n2,uds_ratio_or_lines,n3); 

% ------------------uniform density random 2d sampling------------------
% omega = find(rand(n1*n2*n3,1)<uds_ratio_or_lines);
% sampling_mask = zeros(n1,n2,n3);
% sampling_mask(omega) = 1;

% ------------------uniform density randome x sampling------------------
% raws = round(n1*uds_ratio_or_lines);
% ind_sample = randi(n1,raws,n3);
% sampling_mask = zeros(n1,n2,n3);
% for i = 1:n3
%     sampling_mask(ind_sample(:,i),:,i) = 1;
% end

% ---------------------------radio sampling-----------------------------------
line = uds_ratio_or_lines;
[T3D] = strucrand(n1,n2,n3,line);
sampling_mask = fftshift(T3D);
undersampling_ratio = sum(sampling_mask(:))./(n1*n2*n3);
%% obeserve data
b = sampling_mask.*X;
b = add_noisy(b,20);
%% TNN+MNN in admm
disp('============================');
disp('Recon using TNN+MNN');
disp('------------------------------------------------');
param.lambda1 = 0.1; % with normalization, recommand 0.1, or, without, recommand 2
param.lambda2 = 0.1; 
param.mu = 0.1;
X_tnnmnn = tmnnAlg_fast( b,sampling_mask,X,param,0 );
snr1 = SNR(X,X_tnnmnn);
disp(['lambda = ',num2str(param.lambda1),',',num2str(param.lambda2),',    ------>    SNR = ',num2str(snr1)]);
%% TNN fast in admm
disp('============================');
disp('Recon using TNN');
disp('------------------------------------------------');
param.lambda_tnn = 0.3;
param.mu_tnn = 0.1;
X_tnn = tnnAlg_fast( b,sampling_mask,X,param,0 );
snr2 = SNR(X,X_tnn);
disp(['lambda = ',num2str(param.lambda_tnn),',    ------>    SNR = ',num2str(snr2)]);
%% MNN fast in admm
disp('============================');
disp('Recon using MNN');
disp('------------------------------------------------');
param.lambda_mnn = 0.3;
param.mu_mnn = 0.1;
X_mnn = mnnAlg_fast( b,sampling_mask,X,param,0 );
snr3 = SNR(X,X_mnn);
disp(['lambda = ',num2str(param.lambda_mnn),',    ------>    SNR = ',num2str(snr3)]);