function A = A_fhp3D_p(z)
global n1 n2 q nc S1 coil_sens
z=double(z); S = double(S1); coils = double(coil_sens);
Q= zeros(n1,n2,q,nc);
if length(z(1,1,:))==1
    z=repmat(z,[1,1,q]);
    for p =1 : nc,
        Q(:,:,:,p)= (1/sqrt(n1*n2))*fft2(z.*squeeze(coils(:,:,:,p)));
    end
else
    for p =1 : nc,
        Q(:,:,:,p)= (1/sqrt(n1*n2))*fft2(z.*squeeze(coils(:,:,:,p)));
    end
end
A= Q(S);
