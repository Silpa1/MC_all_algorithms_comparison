function A = A_fhp3D_modi(z)
global n1 n2 q nc S2 coil_sens
z=double(z); S = double(S2); coils = double(coil_sens);
Q= zeros(n1,n2,1,nc);

    for p =1 : nc,
        Q(:,:,p)= (1/sqrt(n1*n2))*fft2(z.*squeeze(coils(:,:,1,p)));
    end
A= Q(S);
