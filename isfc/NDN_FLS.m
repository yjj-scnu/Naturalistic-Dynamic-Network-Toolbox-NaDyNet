clear all;clc;
rootdir='E:\yjj\scnu_work\matlab_APP\data\sfc\data\ROI_mat\raw';

% method='L1';
TR=1;
% wsize=1;
data=read_2Dmat_2_3DmatrixROITC(rootdir);
N_sub = size(data, 3);
N_time = size(data, 1);
N_roi = size(data, 2);
Nwin = N_time;
% atlastype={'DMN_3mm_11nodes.nii'};
mu=100;

%% dynamic FC
resultdir=[];

dFC_result=[];
for s=1:N_sub

    %is_dcc
    subtc=squeeze(data(:,:,s));%time * ROI
    subtcZ=zscore(subtc);%time * ROI

    fprintf('FLS for sub %s\n', num2str(s));
    Ct2 = yuan_DynamicBC_fls_FC(subtcZ,mu);
    % extract the upper right ISDCC values
    % moving average DCC with window length
    atmp=zeros(size(Ct2,1),size(Ct2,1));
    tmp_dFC_DCCX=zeros(Nwin,length(mat2vec(atmp)));
    for iw=1:Nwin
        tmpr=Ct2(:,:,iw);
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


