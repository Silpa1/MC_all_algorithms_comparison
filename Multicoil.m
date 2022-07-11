clc;
clear all;
close all;
%filenames=  {'Pincat.mat','brain_T2_T1.mat','speech_seq.mat','Cardiac_ocmr_data.mat','lowres_speech.mat','FB_ungated.mat'};


[fid,msg] = fopen('Comparison_error.txt','wt');
[fid2,msg] = fopen('Comparison_sim.txt','wt');
fprintf(fid, '%s(%s) & %s & %s & %s & %s  &  %s  \n','Dataset','Radial','ktslr','L+S-Otazo','L+S-Lin','altGDmin-MRI','altGDmin-MRI2');
fprintf(fid2, '%s(%s) & %s & %s & %s & %s  &  %s   \n','Dataset','Radial','ktslr','L+S-Otazo','L+S-Lin','altGDmin-MRI','altGDmin-MRI2');




global mk m S1 coil_sens n1 n2 q nc
% % [fid,msg] = fopen('Comparison_error.txt','wt');
% % [fid2,msg] = fopen('Comparison_sim.txt','wt');
% % fprintf(fid, '%s(%s) & %s & %s &  %s &  %s &  %s  &  %s  \n','Dataset','Radial','ktslr','L+S-Otazo','L+S-Lin','altGDmin','altGDmin-MRI','altGDmin-MRI2');
% % fprintf(fid2, '%s(%s) & %s & %s &  %s &  %s &  %s &  %s   \n','Dataset','Radial','ktslr','L+S-Otazo','L+S-Lin','altGDmin','altGDmin-MRI','altGDmin-MRI2');

% load('cardiac_perf_R8.mat')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% load multicoil_ungated_cmr_perf.mat
%  filenames={ 'multicoil_ungated_cmr_perf.mat'}
%  load(filenames{1})
%  [~,name,~] = fileparts(filenames{1});
% no_comp=8;% NUMBR OF virtual coils
% [kSpace]= coil_compress_withpca(kSpace,no_comp);% use this for coil compression; else comment it out.
% [nx,ny,nc,q]=size(kSpace);
% n=nx*ny;
% im = sqrt(n)*ifft2(kSpace);
% xm = mean(im,4);
% csm = ismrm_estimate_csm_walsh_modified(xm);
% x_coil_combined_groundTruth = squeeze(sum(im.*repmat(conj(csm),[1 1 1 q]),3));
% for i =1:1:nc
%     kdata(:,:,:,i)=kSpace(:,:,i,:);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load multi_coil_lowres_speech.mat
% no_comp=8;% NUMBR OF virtual coils
% [k]= coil_compress_withpca(k,no_comp);% use this for coil compression; else comment it out.
% [nx,ny,nc,nt]=size(k);
% n=nx*ny;
% q=nt;
% im = sqrt(nx*ny)*ifft2(k);
% xm = mean(im,4);
% csm = ismrm_estimate_csm_walsh_modified(xm);
% x_coil_combined_groundTruth = squeeze(sum(im.*repmat(conj(csm),[1 1 1 nt]),3));
% for i =1:1:nc
%     kdata(:,:,:,i)=k(:,:,i,:);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filenames={ 'brain_T2_T1rho.mat'}
load(filenames{1})
[~,name,~] = fileparts(filenames{1});


[nx,ny,nc,nt]=size(img1);
n=nx*ny;
q=nt;
k= (1/sqrt(nx*ny))*fft2(img1);
no_comp=8;% NUMBR OF virtual coils
[k]= coil_compress_withpca(k,no_comp);% use this for coil compression; else comment it out.
[nx,ny,nc,nt]=size(k);

im = sqrt(nx*ny)*ifft2(k);
xm = mean(im,4);

csm = ismrm_estimate_csm_walsh_modified(xm);

x_coil_combined_groundTruth = squeeze(sum(im.*repmat(conj(csm),[1 1 1 nt]),3));


for i =1:1:nc
    kdata(:,:,:,i)=k(:,:,i,:);
end

% % % % % % % % % % % % % % % % % % % % % % % %







radial=[4,8,16]
for ii=1:1:length(radial)
    n=nx*ny;
x=reshape(x_coil_combined_groundTruth ,[n,q]);
Xtrue=x_coil_combined_groundTruth;
    [n1,n2,q,nc]=size(kdata);
    n=n1*n2;
    samp = goldencart(n1,n2,q,radial(ii));
    % im(samp)
    % load('Xinf.mat')
    %  Xtrue=double(Xinf.perf);
    % Xtrue=double(Image1);
    mask = fftshift(fftshift(samp,1),2);
    % mask=samp;
    mask3=reshape(mask,[n, q]);
    
    mk=[];
    for i=1:1:q
        mk(i)=length(find(logical(mask3(:,i))));
        S2(1:mk(i),i)=double(find(logical(mask3(:,i))));
    end
    mask1=repmat(mask,[1,1,1,nc]);
    S1=find(mask1~=0);
    kdata1=kdata.*mask1;
    y=kdata(S1);
    % b1=ones(n1,n2,nc);
    %  tmp = sqrt(sum(abs((b1)).^2,3));
    % b1c = div0(b1,tmp);
    for i=1:1:q
        coil_sens(:,:,i,:)=csm;
    end
%     masknew=logical(kdata1(:,:,:,1)~=0);
%     rrr=2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % X=reshape(Xtrue,[n,q]);
    % m=max(mk);
    % Y=zeros(m,q);
    %         for k=1:1:q
    %             ksc = reshape((1/sqrt(n1*n2))* fft2( reshape(X(:,k), [n1 n2]) ), [n,1]) ;
    %             Y2(:,k)=ksc;
    %             Y(1:mk(k),k)=double(ksc(S2(1:mk(k),k)));
    %         end
    %        kdata=Y2;
    % y=A_fhp3D_p(Xtrue);
    % z=Xtrue;
    %     for p =1 : nc,
    %         Q(:,:,:,p)= (1/sqrt(n1*n2))*fft2(z.*squeeze(coil_sens(:,:,:,p)));
    %     end
    %     kdata1=Q.*mask1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Ktslr %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     A = @(z)A_fhp3D_p(z);
%     
%     At=@(z)At_fhp3D_p(z);
%     b=y;
%     addpath('ktslr/');
%     tic;
%     
%     
%      x_init = At(b);
%      mu =20e-6; % Regularization parameter
%      opts.mu = mu;
%      opts.betarate = 35;
%      opts.outer_iter =12;
%      opts.inner_iter = 50;
%      opts.beta=1e3; %
%     
%     U=x_init; [m n d] = size(U);
%     Lam= zeros(m,n,d);
%     
%     o=0;
%     earray=zeros(1,opts.outer_iter*opts.inner_iter); cost=zeros(1,opts.outer_iter*opts.inner_iter);
%     
%     fac=length((find(b~=0)));
%     U = double(U); b=double(b); fac=double(fac);
%     for out = single(1:opts.outer_iter),
%         o=o+1;
%         for in = single(1:opts.inner_iter)
%     
%        z = temp_FT(U);
%     
%     
%     z1 = (abs(z)-1/opts.beta); z1 = z1.*(z1>0);
%     z2 = abs(z) + (abs(z)<1e-12);
%     z = z1.*z./z2;
%     z = itemp_FT(z);
%     
%     
%     
%      [U,earray1] = xupdateCG(b,A, At,z,opts,U,Lam, 1e-7,10);
%     
%      e = A(U) - b;
%     
%      e_ft = temp_FT(U);
%     
%     
%      cost = [cost, sum(abs(e(:)).^2)  +  sum(abs(e_ft(:)))*opts.mu];
%     
%             if in>1
%             if abs(cost(end) - cost(end-1))/abs(cost(end)) < 1e-3
%                 break;
%             end
%             end
%     
%     
%     
%         end
%         opts.beta=opts.beta*opts.betarate;
%     end
%     Zhat_ktslr=U;
%     Time_ktslr=toc;
%      similarity_index=[];
%             for i =1:1:q
%                 mssim=ssim(abs(Zhat_ktslr(:,:,i)/max(max(Zhat_ktslr(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
%                 similarity_index(i)=mssim;
%             end
%             sim_ktslr=min(similarity_index);
    
%     NMSE_ktslr=RMSE_modi(Zhat_ktslr,Xtrue)
    % %%%%%%%%%%%%%%%%%%%%   LplusS Otazo (cardiac perfusion parameters) modified  %%%%%%%%%%%%%
%              A = @(z)A_fhp3D_p(z);
%              At = @(z) At_fhp3D_p(z);
%             param.A=A;
%             param.At=At;
%             param.d = y;
%             param.T=TempFFT(3);
%             tic;
%             param.lambda_L=0.01;
%             param.lambda_S=0.01;
%             param.nite=50;
%             param.tol=0.0025;
%             M=At(param.d);
%             M=reshape(M,[n1*n2,q]);
%             Lpre=M;
%             S=zeros(n1*n2,q);
%             ite=0;
%             while(1)
%                 ite=ite+1;
%                 M0=M;
%                 [Ut,St,Vt]=svd(M-S,0);
%                 St=diag(SoftThresh(diag(St),St(1)*param.lambda_L));
%                 L=Ut*St*Vt';
%                 S=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[n1,n2,q]),param.lambda_S)),[n1*n2,q]);
%                 resk=param.A(reshape(L+S,[n1,n2,q]))-param.d;
%                 M=L+S-reshape(param.At(resk),[n1*n2,q]);
%                 Lpre=L;
%                 tmp2=param.T*reshape(S,[n1,n2,q]);
%                 if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
%             end
%             Xhat_LpS1=L+S;
%             Xhat_LpS=reshape(Xhat_LpS1,[n1,n2,q]);
%             Time_Otazo= toc;
%             NMSE_Otazo=RMSE_modi(Xhat_LpS,Xtrue);
%             similarity_index=[];
%             for i =1:1:q
%                 mssim=ssim(abs(Xhat_LpS(:,:,i)/max(max(Xhat_LpS(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
%                 similarity_index(i)=mssim;
%             end
%             sim_Otazo=min(similarity_index)
    
    %%%%%%%%%%%%%%%%%%%%   LplusS Otazo (cardiac perfusion parameters %%%%%%%%%%%%%
    % tic;
    % L+S reconstruction ******************************************************
    % kdata1=fftshift(kdata1);
    % for i=1 : nc
    %     kdata2(:,:,:,i)=fftshift(kdata1(:,:,:,i));
    % end
    % param.E=Emat_xyt(kdata1(:,:,:,1)~=0,squeeze(coil_sens(:,:,1,:)));
    % % % param.E=Emat_xyt(mask,squeeze(coil_sens(:,:,1,:)));
    % % % param.d=param.E*Xtrue;
    % % % param.T=TempFFT(3);
    % % % param.lambda_L=0.01;
    % % % param.lambda_S=0.01;
    % % % param.nite=50;
    % % % param.tol=0.0025;
    % % % [L,S] = lps_ist(param);
    % % %
    % % % Zhat_LplusS_otazo=double(L+S);
    % % %
    % % % similarity_index=[];
    % % % for i =1:1:q
    % % %     mssim=ssim(abs(Zhat_LplusS_otazo(:,:,i)/max(max(Zhat_LplusS_otazo(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
    % % %     similarity_index(i)=mssim;
    % % % end
    % % % sim_LplusS_lin=min(similarity_index);
    % % % Time_LplusS_otazo=toc;
    % % % NMSE_LplusS_otazo=RMSE_modi(Zhat_LplusS_otazo,Xtrue)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LplusS jeff (cardiac perfusion parameters) %%%%%%%%%%%
    
    
%     
% %     tmp = sqrt(sum(abs((csm)).^2,3));
% %     b1c = div0(csm,tmp);
% %     
% % %     param.E=getE(csm,q,'samp',samp);
% % %        kdata2=param.E*Xtrue;
% addpath('LplusS_lin\');
%       A = @(z)A_fhp3D_p(z);
%              At = @(z) At_fhp3D_p(z);
%             param.A=A;
%             param.At=At;
%             param.d = y;
%     tic;
%  
%     param.T = getT(n1,n2,q);
%     param.nite=10;
%     param.scaleL = 130/1.2775;
%     param.scaleS = 1/1.2775;
%     param.lambda_L=0.01;
%     param.lambda_S=0.01*param.scaleS;
%      param.Xinf = zeros(n1*n2,q);
%     
%     [L_pogm,S_pogm,xdiff_pogm,cost_pogm,time_pogm,rankL_pogm] = PGM(param,'pogmS',1,'pogmL',1);
%     L = L_pogm;S = S_pogm;
%     Zhat_LplusS_lin=double(L+S);
%     Time_lin=toc;
%     similarity_index=[];
%     for i =1:1:q
%         mssim=ssim(abs(Zhat_LplusS_lin(:,:,i)/max(max(Zhat_LplusS_lin(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
%         similarity_index(i)=mssim;
%     end
%     sim_lin=min(similarity_index);
%     
%     NMSE_lin=RMSE_modi(Zhat_LplusS_lin,Xtrue)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%altGDMin %%%%%%%%%%%%%
    
    % tic;
    % T=70;
    % [U_init, eta,Uhat, Bhat, ErrorU, ErrorUF, ErrorX, ExeTime] = AltGDMin(S2,n1,n2,nc,  T,x,y,coil_sens,S1);
    % Xhat=Uhat*Bhat;
    % Zhat_GD=double(reshape(Xhat,[n1,n2,q]));
    % Time_GD=toc;
    % similarity_index=[];
    % for i =1:1:q
    %     mssim=ssim(abs(Zhat_GD(:,:,i)/max(max(Zhat_GD(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
    %     similarity_index(i)=mssim;
    % end
    % sim_GD=min(similarity_index);
    %
    % NMSE_GD=RMSE_modi(Zhat_GD,Xtrue)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%altGDMin_MRI%%%%%%%%%%%%%
     tic;
     addpath('altGDmin_MRI_MRI2\');
    [zbar_hat,flag,resNE,iter] = cgls(@A_fhp3D_p,@At_fhp3D_p, y,0,1e-36,10);
    
    Ytemp=A_fhp3D_p(zbar_hat);
    
    Yinter=y-Ytemp;
    
   
    T=70;
    [U_init, eta,Uhat, Bhat, ErrorU, ErrorUF, ErrorX, ExeTime] = AltGDMin(S2,n1,n2,nc,  T,x,Yinter,coil_sens,S1);
    Xhat=Uhat*Bhat;
    Xhat_mat=reshape(Xhat,[n1,n2,q]);
    Yhat_hat=y-A_fhp3D_p(Xhat_mat+zbar_hat);
    Ehat=[];
    sm=sum(mk);
    global kk S2
    
    for kk=1:1:q
        tmp1=[];
        y1=[];
        S2=[];
        for i=1:1:nc
            S2=find(mask1(:,:,kk,:)~=0);
            
            y1=[y1;Yhat_hat(((i-1)*sm)+sum(mk(1:(kk-1)))+1:((i-1)*sm)+sum(mk(1:(kk-1)))+mk(kk))];
            
        end
        Ehat(:,:,kk)=cgls_modi(@A_fhp3D_modi,@At_fhp3D_modi, y1 ,0,1e-36,3);
    end
    Zhat_MRI=double(reshape(Xhat,[n1,n2,q])+repmat(zbar_hat,[1,1,q])+Ehat);
    Time_MRI=toc;
    similarity_index=[];
    for i =1:1:q
        mssim=ssim(abs(Zhat_MRI(:,:,i)/max(max(Zhat_MRI(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
        similarity_index(i)=mssim;
    end
    sim_MRI=min(similarity_index);
    NMSE_MRI=RMSE_modi(Zhat_MRI,Xtrue)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% altGDMin_MRI2%%%%%%%%%%%%%%%%%
    [zbar_hat,flag,resNE,iter] = cgls(@A_fhp3D_p,@At_fhp3D_p, y,0,1e-36,10);
    
    Ytemp=A_fhp3D_p(zbar_hat);
    
    Yinter=y-Ytemp;
    L=[];
    tic;
    T=70;
    [U_init, eta,Uhat, Bhat, ErrorU, ErrorUF, ErrorX, ExeTime] = AltGDMin(S2,n1,n2,nc,  T,x,Yinter,coil_sens,S1);
    Xhat=Uhat*Bhat;
    Xhat_mat=reshape(Xhat,[n1,n2,q]);
    Yhat_hat=y-A_fhp3D_p(Xhat_mat+zbar_hat);
    A = @(z)A_fhp3D_p(z);
    At = @(z) At_fhp3D_p(z);
    param.A=A;
    param.At=At;
    param.d = y;
    param.T=TempFFT(3);
    param.lambda_L=0.01;
    param.nite=10;
    param.tol=0.0025;
    M=At(param.d);
    M=reshape(M,[n1*n2,q]);
    Lpre=M;
    Ehat=zeros(n1*n2,q);
    L(:,1:q)=reshape(Xhat_mat+zbar_hat,[n1*n2,q]);
    param.lambda_S=0.001*max(max(abs(M-L)));
    ite=0;
    while(1)
        ite=ite+1;
        M0=M;
        Ehat=reshape(param.T'*(SoftThresh(param.T*reshape(M-Lpre,[n1,n2,q]),param.lambda_S)),[n1*n2,q]);
        resk=param.A(reshape(L+Ehat,[n1,n2,q]))-param.d;
        M=L+Ehat-reshape(param.At(resk),[n1*n2,q]);
        Lpre=L;
        tmp2=param.T*reshape(Ehat,[n1,n2,q]);
        if (ite > param.nite) || (norm(M(:)-M0(:))<param.tol*norm(M0(:))), break;end
    end
    
    Zhat=L+Ehat;
    Zhat_MRI2=reshape(Zhat,n1,n2,q);
    Time_MRI2=  toc;
    NMSE_MRI2=RMSE_modi(Zhat_MRI2,Xtrue);
    similarity_index=[];
    for i =1:1:q
        mssim=ssim(abs(Zhat_MRI2(:,:,i)/max(max(Zhat_MRI2(:,:,i)))),abs(Xtrue(:,:,i)/max(max(Xtrue(:,:,i)))));
        similarity_index(i)=mssim;
    end
    sim_MRI2=min(similarity_index);
    fprintf(fid, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f) & %8.4f (%5.2f) \n', name, radial(ii),NMSE_ktslr,Time_ktslr,NMSE_Otazo,Time_Otazo,NMSE_lin, Time_lin,NMSE_MRI,Time_MRI,NMSE_MRI2,Time_MRI2);
    fprintf(fid2, '%s(%d) & %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)& %8.4f (%5.2f)  \n', name, radial(ii),sim_ktslr,Time_ktslr,sim_Otazo,Time_Otazo,sim_lin,Time_lin,sim_MRI,Time_MRI,sim_MRI2,Time_MRI2);
        
end
fclose(fid);
fclose(fid2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(1,2,1);imagesc(abs(X_image(:,:,1)));title('Original');colormap(gray);
% subplot(1,2,2);imagesc(abs(Zhat_MRI(:,:,1)));title('AltGDmin');colormap(gray);
% figure(2);plot(ErrorUF)


function y=SoftThresh(x,p)
y=(abs(x)-p).*x./abs(x).*(abs(x)>p);
y(isnan(y))=0;
end
