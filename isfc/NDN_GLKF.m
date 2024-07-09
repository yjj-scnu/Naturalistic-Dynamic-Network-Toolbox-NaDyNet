clear all;clc;
rootdir='E:\yjj\scnu_work\matlab_APP\data\sfc\data\ROI_mat\raw';
%（nt *  nr * nsub）
data=read_2Dmat_2_3DmatrixROITC(rootdir);
N_sub = size(data, 3);
N_time = size(data, 1);
N_roi = size(data, 2);

Nwin = N_time-6;
TR=1;
pKF=6;
ucKF=0.03;

%% dynamic FC
dFC_result=[];
for s=1:N_sub

    %is_dcc
    %[tmp_dFC]=pp_ReHo_dALFF_dFC_gift(subtc,method,TR,wsize);%trme * ROI paris, 2D. r*(r-1)/2
    %             [~,Ct2,~,~] = DCC_X(subtc2,allpair, parallel);
    %             Ct2 = yuan_DynamicBC_fls_FC(subtc2Z,mu);
    subtc=squeeze(data(:,:,s));%time * ROI
    subtcZ=zscore(subtc);%time * ROI
    fprintf('GLKF for sub %s\n', num2str(s));
    YKF(1,:,:)=subtcZ';

    FKF = dynet_SSM_KF(YKF,pKF,ucKF);
    for i=1:N_time
        FKFR(:,:,i)=icatb_corrcov(squeeze(FKF.R(:,:,i)));
    end
    % moving average DCC with window length extract the upper right ISDCC values
    atmp=zeros(size(FKFR,1),size(FKFR,1));
    tmp_dFC_DCCX=zeros(Nwin,length(mat2vec(atmp)));
    for iw=1:Nwin
        tmpr=FKFR(:,:,iw+6);
        tmp_dFC_DCCX(iw,:)=mat2vec(squeeze(tmpr));
    end
    tmp_dFC=tmp_dFC_DCCX;
    DEV = std(tmp_dFC, [], 2);%STD OF NODE
    [xmax, imax, xmin, imin] = icatb_extrema(DEV);%local maxima in FC variance
    pIND = sort(imax);%?
    k1_peaks(s) = length(pIND);%?
    SP{s,1} = tmp_dFC(pIND, :);%Subsampling
    dFC_result=[dFC_result;tmp_dFC];
end%s
cd(resultdir)
save('SP.mat','SP','-v7.3')

save('dFC_result.mat','dFC_result','-v7.3')

