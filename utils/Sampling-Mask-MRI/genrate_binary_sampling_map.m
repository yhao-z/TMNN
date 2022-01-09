function sampling_mat=genrate_binary_sampling_map(n1,n2,undersampling_ratio,n3)
% This function generates the binary under-sampling  masks for each image, low
% frequencies have higher probabilty to be chosen.
% sampling_mat is 3D(N,N,L)
sampling_mat= false(n1,n2,n3);
pdf_vardens_cut=genPDF([n1,n2],9,undersampling_ratio,2,0,0);
for i=1:n3
    r_mat=rand(n1,n2);  
    pdf_vardens2=(r_mat.*pdf_vardens_cut);
    pdf_vardens3=pdf_vardens2(:);
    [~,b]=sort(pdf_vardens3);
    b=flipud(b);
    threshold_for_sampling=pdf_vardens3(b(round(undersampling_ratio*length(b))));
    pdf_vardens4=zeros(n1,n2);
    pdf_vardens4(pdf_vardens2>=threshold_for_sampling)=1;
    sampling_mat(1:n1,1:n2,i)=logical(ifftshift(pdf_vardens4));
end
