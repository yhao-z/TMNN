% NOTE: in this code, we compare the fast algorithm in frequency domain
% with the classical algorithm in spatiotemporal domain and show the
% difference in cost times. For precision, we recommand you to annotate the
% cost calculation and the plot codes. 
%% Initialize:
clear;clc;
addpath(genpath(pwd));
%% 4 datasets
% --------------------------------------------------
% load('E:\Yhao\DATA\dmri\aperiodic_pincat.mat','new');
% X = fft2(new);
% --------------------------------------------------
% load('E:\Yhao\DATA\dmri\invivo_perfusion.mat')
% x = x./max(abs(x(:)));
% X = fft2c(x);
% --------------------------------------------------
load('E:\Yhao\DATA\dmri\data_tmi.mat','seq');
x = seq;
X = fft2c(x); 
% --------------------------------------------------
% load('E:\Yhao\DATA\dmri\breathing.mat','Data_xy_t');
% X = fft2(Data_xy_t);

[n1,n2,n3] = size(X);
%% normalize
maxX = max(abs(X(:)));
X = X./maxX;
x = x./maxX;
%% 4 sampling mask
uds_ratio_or_lines =30;
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
%% TMNN fast
disp('============================');
disp('Recon using TMNN_fast');
disp('------------------------------------------------');
param.lambda1 = 0.1;
param.lambda2 = 0.1;
param.mu = 0.1;
tic;
[X_tmnn_fast, ~, ~, ~] = tmnnAlg_fast( b,sampling_mask,X,param,1 );
times_tmnn_fast = toc;
snr2 = SNR(X,X_tmnn_fast);
disp(['lambda = ',num2str(param.lambda1),', ',num2str(param.lambda2),',    ------>    SNR = ',num2str(snr2)]);
disp(['Fast TMNN cost ',num2str(times_tmnn_fast),' s']);
%% TMNN
disp('============================');
disp('Recon using TMNN');
disp('------------------------------------------------');
param.lambda1 = 0.1;
param.lambda2 = 0.1;
param.mu = 0.1;
tic;
[x_tmnn, ~, ~, ~] = tmnnAlg( b,sampling_mask,x,param,1 );
times_tmnn=toc;
snr1 = SNR(x,x_tmnn);
disp(['lambda = ',num2str(param.lambda1),', ',num2str(param.lambda2),',    ------>    SNR = ',num2str(snr1)]);
disp(['Classical TMNN cost ',num2str(times_tmnn),' s']);
