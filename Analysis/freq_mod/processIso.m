function [data] = processIso(params, data)
%Process Dual Color Photometry
%
%   [data] = processDual(data,params)
%
%   Description: This function is designed to process photometry
%   experiments with an isosbestic control. The code takes the photometry
%   channel and asks the user for the modulation frequencies of the
%   excitation 470 and the isosbestic control 405
%
%   Input:
%   - data - A data structure specific to the Tritsch Lab. Created using
%   the convertH5_FP script
%   - params - A structure created from a variant of the processParams
%   script
%
%   Output:
%   - data - Updated data structure containing processed data
%
%   Author: Pratik Mistry 2019

%Pull parameters required for this analysis
nAcq = length(data.acq);
%Filter Properties
lpCut = params.FP.lpCut; filtOrder = params.FP.filtOrder;
%General downsampling parameter
dsRate = params.dsRate;
%Baselining parameters
interpType = params.FP.interpType;
fitType = params.FP.fitType; winSize = params.FP.winSize;
winOv = params.FP.winOv;
basePrc = params.FP.basePrc;
%Demodulation-specific property
modFreq = params.FP.modFreq;
sigEdge = params.FP.sigEdge;
%Downsampling property
rawFs = data.acq.Fs;
Fs = rawFs/dsRate;
data.final.Fs = Fs;
%This outer for-loop goes performs the analysis on each sweep acquired
%during the experiment
for x = 1:nAcq
    Ls = length(data.acq(x).time);
    L = 1:Ls;
    nFP = length(data.acq(x).FP);
    refSig = data.acq(x).refSig;
    data.final(x).FP = cell(nFP,1); data.final(x).nbFP = cell(nFP,1); data.final(x).FPbaseline = cell(nFP,1);
    %The for loop will go through all FP traces assuming all of them were
    %recorded using an isosbestic control.
    for y = 1:nFP
        rawFP = data.acq(x).FP{y,1};
        isoFreq = modFreq(1); excFreq = modFreq(2);
        isoRef = findRef(isoFreq,refSig,rawFs); excRef = findRef(excFreq,refSig,rawFs);
        isoDemod = digitalLIA(rawFP,isoRef,isoFreq,rawFs,lpCut,filtOrder);
        excDemod = digitalLIA(rawFP,excRef,excFreq,rawFs,lpCut,filtOrder);
        if sigEdge ~= 0
            isoDemod = isoDemod((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
            excDemod = excDemod((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
        end
        excDemod = downsample(excDemod, dsRate);
        isoDemod = downsample(isoDemod, dsRate);
        data.final(x).nbFP{y} = excDemod;
        data.final(x).iso{y} = isoDemod;
        %Use isosbestic signal to baseline the signal using linear regression
        [FP,baseline] = linregFP(isoDemod,excDemod,basePrc);
%         [FP,baseline] = baselineFP(excDemod,interpType,fitType,basePrc,winSize,winOv,Fs);
        data.final(x).FP{y} = FP;
        data.final(x).FPbaseline{y} = baseline;
    end
    %Create the time vector based on the new length of the photometry
    %signal and store new photometry signal
    if sigEdge ~= 0
        L = L((sigEdge*rawFs)+1:end-(sigEdge*rawFs));
        Ls = length(L);
    end
    Ls = length(1:dsRate:Ls); timeVec = [1:Ls]/Fs;
    data.final(x).time = timeVec';
end
end
