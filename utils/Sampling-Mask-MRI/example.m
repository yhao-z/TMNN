%% initial
% set undersampling_ratio a scalar between [0 1] except the radio sampling
% mode, which set undersampling_ratio a integer, like 16, 30, 60 and etc. .
% [n1,n2,n3] is the dimension of the data.
undersampling_ratio=30;
n1 = 120; n2=240; n3=50;
%% variable density random 2d sampling
sampling_mask = genrate_binary_sampling_map(n1,n2,undersampling_ratio,n3); 
%% uniform density random 2d sampling
omega = find(rand(n1*n2*n3,1)<undersampling_ratio);
sampling_mask = zeros(n1,n2,n3);
sampling_mask(omega) = 1;
%% variable density randome x sampling
sampling_mask = genrate_ylines_sampling_map(n1,n2,undersampling_ratio,n3); 
%% uniform density randome x sampling
raws = round(n1*undersampling_ratio);
ind_sample = randi(n1,raws,n3);
sampling_mask = zeros(n1,n2,n3);
for i = 1:n3
    sampling_mask(ind_sample(:,i),:,i) = 1;
end
%% radio sampling
line = undersampling_ratio;
[T3D] = strucrand(n1,n2,n3,line); % Generate the (kx,ky)-t sampling pattern; Each frame has uniformly spaced radial rays, with random rotations across frames
sampling_mask = fftshift(T3D);
undersampling_ratio = sum(sampling_mask(:))./(n1*n2*n3);