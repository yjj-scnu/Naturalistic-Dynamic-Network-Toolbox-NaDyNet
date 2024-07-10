clear all;
clc;
rootdir='E:\yjj\scnu_work\matlab_APP\data\sfc\data\ROI_mat\raw';

methodType = 'ISDCC';% Valid values are "SWFC" or "SWISFC".
% ISA_type (str): The type of intersubject analysis. Valid values are "LOO_regress" or "LOO".
ISA_type = 'LOO';
%         load data （nt *  nr * nsub）
data=read_2Dmat_2_3DmatrixROITC(rootdir);
N_sub = size(data, 3);
N_time = size(data, 1);
N_roi = size(data, 2);
Nwin = N_time;
allpair = 0; parallel = 0;

%% dynamic FC
savedDirdir=[];
dFC_result=[];

for s=1:N_sub
    subtc=squeeze(data(:,:,s));%time * ROI

    if isequal(methodType, 'ISDCC') % isxxx

        if ~exist("ISA_type","var")
            error("Argument ISA_type should not be empty, it must be assgind ")
        end
        if isequal(ISA_type, "LOO")
            subR=data;
            subR(:,:,s)=[];
            subtc2=[subtc,squeeze(mean(subR,3))];
            subtc2Z=zscore(subtc2);
            fprintf('ISDCC for sub %s\n', num2str(s));
            [~,Ct2,~,~] = DCC_X(subtc2Z,allpair, parallel);
            % extract the upper right ISDCC values
            ISCt2=Ct2(1:N_roi, N_roi + 1 : N_roi * 2, :);

            % moving average DCC with window length
            tmp_dFC_DCCX=zeros(Nwin,N_roi*N_roi);
            for iw=1:Nwin
                tmpr=ISCt2(:,:,iw);
                tmp_dFC_DCCX(iw,:)=mat2vec_Asym(tmpr);
            end

        else
            error("暂时没有开发出来")
        end
    else% xxx
        subtcZ=zscore(subtc);%time * ROI
        fprintf('DCC for sub %s\n', num2str(s));
        [~,Ct2,~,~] = DCC_X(subtcZ,allpair, parallel);
        % extract the upper right ISDCC values
        % moving average DCC with window length
        atmp=zeros(size(Ct2,1),size(Ct2,1));
        tmp_dFC_DCCX=zeros(Nwin,length(mat2vec(atmp)));
        for iw=1:Nwin
            tmpr=Ct2(:,:,iw);
            tmp_dFC_DCCX(iw,:)=mat2vec(squeeze(tmpr));
        end

    end

    tmp_dFC=tmp_dFC_DCCX;
    DEV = std(tmp_dFC, [], 2);%STD OF NODE
    [xmax, imax, xmin, imin] = icatb_extrema(DEV);%local maxima in FC variance
    pIND = sort(imax);%?
    k1_peaks(s) = length(pIND);%?
    SP{s,1} = tmp_dFC(pIND, :);%Subsampling
    dFC_result=[dFC_result; tmp_dFC];
end%s

cd(savedDirdir)
save('SP.mat','SP','-v7.3')
save('dFC_result.mat','dFC_result','-v7.3')


