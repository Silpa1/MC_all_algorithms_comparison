function [U_init, eta,Uhat, Bhat, ErrorU, ErrorUF, ErrorX, ExeTime] = AltGDMin(S,n1,n2,nc,  T,X,y11,coil_sens,S1)
global   mk  q S n1 n2 m nc S1 coil_sens
ErrorU = [];
ErrorUF = [];
ErrorX = [];
ExeTime = [];
%tic;
[U_init,r] = initAltGDMin( n1,n2,nc,y11,coil_sens,S1);
n=n1*n2;
% [m, q] = size(S);
Uhat = U_init;
% [Ustar, Sstar, Vstar]=svds(X,r);
% ErrorU = [abs(sin(subspace(Ustar, Uhat))), ErrorU];
% ErrorUF(1) =( norm((Uhat - Ustar*(Ustar'*Uhat)), 'fro')/sqrt(r))

sm=sum(mk);
for t = 1 : T
   Uhat_mat=reshape(Uhat,[n1,n2,r]);
% AU=A_Uhat(Uhat_mat,coil_sens,n1,n2,r,nc);
% tic;
%     coil2=coil_sens(:,:,1:r,:);
% 
%     for k =  1 : q
%         tmp1=[];
%         y=[];
%         for i=1:1:nc
%             tmp=zeros(n,r);
%           %%  tmp =(1/ sqrt(n1*n2))* reshape( fft2(reshape(Uhat.*coil(:,i),[n1 n2,r])), [n,r]); %out = tmp(mask_k)
%             tmp =(1/ sqrt(n1*n2))* reshape( fft2(reshape(Uhat,[n1 n2,r]).*coil_sens(:,:,1,i)), [n,r]); %out = tmp(mask_k)
%             tmp1=[tmp1;tmp(S(1:mk(k),k),:)];
%             y=[y;y11(((i-1)*sm)+sum(mk(1:(k-1)))+1:((i-1)*sm)+sum(mk(1:(k-1)))+mk(k))];
%            
%         end
%         Bhat(:, k) = tmp1\y;
%         
%     end
%   Time1=toc;

    for j=1 : r
        Utemp=repmat(Uhat_mat(:,:,j),[1,1,q]);
        AU(:,j)=A_fhp3D_p(Utemp);
    end
    for k=1:q
        AUk=[];
        y=[];
    for i=1:nc
        AUk=[AUk;AU(((i-1)*sm)+sum(mk(1:(k-1)))+1:((i-1)*sm)+sum(mk(1:(k-1)))+mk(k),:)];
        y=[y;y11(((i-1)*sm)+sum(mk(1:(k-1)))+1:((i-1)*sm)+sum(mk(1:(k-1)))+mk(k))];
    end
    Bhat(:, k) = AUk\y;
    end
    Xhat = Uhat * Bhat;
    
    
    Grad_U = zeros(n, r);
    
    Inter=zeros(nc*m,q);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Xhat_mat=reshape(Xhat,[n1,n2,q]);
    yinter=A_fhp3D_p(Xhat_mat);
    ydif=yinter-y11;
    Grad_U= (reshape(At_fhp3D_p(ydif),[n1*n2,q]))*(Bhat');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    if t==1
        eta=1/(7*norm(Grad_U));
    end
    Uhat_t0=Uhat;
    Uhat = Uhat - eta * Grad_U;
    [Qu,~]  =  qr(Uhat,0);
    Uhat  =  Qu(:, 1:r);
    Uhat_t1=Uhat;
    Subspace_d= ( norm((Uhat_t0 - Uhat_t1*(Uhat_t1'*Uhat_t0)), 'fro')/sqrt(r));
    if  (Subspace_d <=.01)
        break;
    end
%     ErrorU = [ ErrorU, abs(sin(subspace(Ustar, Uhat)))];
%     ErrorUF(t+1) = ( norm((Uhat - Ustar*(Ustar'*Uhat)), 'fro')/sqrt(r));

end

for j=1 : r
        Utemp=repmat(Uhat_mat(:,:,j),[1,1,q]);
        AU(:,j)=A_fhp3D_p(Utemp);
    end
    for k=1:q
        AUk=[];
        y=[];
    for i=1:nc
        AUk=[AUk;AU(((i-1)*sm)+sum(mk(1:(k-1)))+1:((i-1)*sm)+sum(mk(1:(k-1)))+mk(k),:)];
        y=[y;y11(((i-1)*sm)+sum(mk(1:(k-1)))+1:((i-1)*sm)+sum(mk(1:(k-1)))+mk(k))];
    end
    Bhat(:, k) = AUk\y;

end