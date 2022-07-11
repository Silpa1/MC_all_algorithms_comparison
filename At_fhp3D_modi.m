function At = At_fhp3D_modi(z)
global q S2 n1 n2 mk m coil_sens nc kk

z=double(z); S2 = double(S2); coils = double(coil_sens);
data_img  = zeros(n1,n2,nc);
T = zeros(n1,n2,nc);
T(S2) = z ;
for p = 1: nc,
    %for t = 1: n3,
    data_img (:,:,p) =(sqrt (n1*n2))*ifft2(T(:,:,p)).*(conj(coils(:,:,1,p))); 
    %end
end

%At = data_img(:,:,:,1)+data_img(:,:,:,2)+data_img(:,:,:,3)+data_img(:,:,:,4);
At=sum(data_img,3);