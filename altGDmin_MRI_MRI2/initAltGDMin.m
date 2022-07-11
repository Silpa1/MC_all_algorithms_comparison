function [Uhat,r] = initAltGDMin(n1,n2,nc,y11,coil_sens, S1)
global mk S n1 n2  nc S1 coil_sens
Cy = 6;
[m,q] = size(S);
n=n1*n2;

threshold1 = ( Cy * norm(y11) ) / sqrt(nc*sum(mk));
Y_trunc1 = y11;
Y_trunc1(abs(y11) > threshold1) = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X2=At_fhp3D_p(Y_trunc1);
X3=reshape(X2,[n1*n2,q]);
% tic;
% for k=1:1:q
%     X01(:,k)=X3(:,k)/sqrt(mk(k)*m);
% end
% Time=toc;
% tic;
% X0=X3./(sqrt(mk));
% X0=X0/sqrt(m);
% Time1=toc;
% tic;
X01=X3./(sqrt(mk *m));
% Time2=toc;
%%%%%%%%%%% Automated r %%%%%%%%%%%%%%%%%%%
r_big=floor(min([n/10,q/10,m*nc/10]));
[Un,Sigman,Vn]=svds(double(X01),r_big);
SS=diag(Sigman);
E=sum(SS.^2);
Esum=0;
for i=1:1:r_big
    Esum=Esum+((SS(i))^2);
    if Esum >(E*0.85)
        break
    end
end
r=i+1;
r=min(r,r_big);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Uhat=Un(:,1:r);
rrrr=1;
