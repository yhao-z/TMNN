function sampling_mat=genrate_ylines_sampling_map(n1,n2,undersampling_ratio,n3)
% This function generates the ylines under-sampling  masks for each image, low
% frequencies have higher probabilty to be chosen.
% sampling_mat is 3D(n1,n2,n3)
sampling_mat= false(n1,n2,n3);
pdf_vardens_cut=genPDF([n1,1],9,undersampling_ratio,2,0,0);
pdf_vardens_cut = pdf_vardens_cut(:);
for i=1:n3
    r_mat=rand(n1,1);  
    pdf_vardens2=(r_mat.*pdf_vardens_cut);
    [~,b]=sort(pdf_vardens2);
    b=flipud(b);
    threshold_for_sampling=pdf_vardens2(b(round(undersampling_ratio*length(b))));
    pdf_vardens3=zeros(n1,1);
    pdf_vardens3(pdf_vardens2>=threshold_for_sampling)=1;
    sampling_mat(:,:,i)=repmat(logical(ifftshift(pdf_vardens3)) , [1 n2]);
end
