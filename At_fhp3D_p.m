function At = At_fhp3D_p(z)
global q S1 n1 n2 mk m coil_sens nc

z=double(z); S1 = double(S1); coils = double(coil_sens);
data_img  = zeros(n1,n2,q,nc);
T = zeros(n1,n2,q,nc);
T(S1) = z ;
for p = 1: nc,
    %for t = 1: n3,
    data_img (:,:,:,p) =(sqrt (n1*n2))*ifft2(T(:,:,:,p)).*(conj(coils(:,:,:,p))); 
    %end
end

%At = data_img(:,:,:,1)+data_img(:,:,:,2)+data_img(:,:,:,3)+data_img(:,:,:,4);
At=sum(data_img,4);