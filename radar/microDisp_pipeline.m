%% Pipeline
function [distances, GTErrors1, GTErrors2] = microDisp_pipeline(start_and_end_delay, experiments)
%start_and_end_delay are number of frames that the object is held still at
%the beginning/end of a frame. 
%experiments should be the number of the experiment to analyze

%Will return:
%distances(fileNo, tagNo, sampleNo) of estimated distances moved. To get a single
%estimate, call mean(distances, 2)
%GTErrors1 is distanes - errors
%GTErrors2 is distances + errors, use either depending on motion stage direction
%Assumes ground truth files are measured in cms
%------------------------------------------------ --------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
c0 = 1/sqrt(8.85e-12*4.*pi.*1e-7);
CurPath = pwd();
addpath([CurPath,'/DemoRadBasic/DemoRadUsb']);
addpath([CurPath,'/DemoRadBasic/Class']);
switchPer=[625*10^-6]; %#ok<*NBRAK2> %625, 500
% minRD=1000;
% maxRD=-1000;
% minR=1000;
% maxR=-1000;
% minAng=1000;
% maxAng=-1000;

bkndSubtract=0;

filePre='./Data/';
c=3e8;
constOffset = 0;

%1, 3, 4, 7, 9, 20, 21, 22, 23, 25, 27, 28
    Brd         =   TinyRad();
    skip=0;
    loadGT = 1;
    if experiments==1
        maxRang=15;
        foldern='blockages-diatomite-1';
        switchPer=[625]*10^-6; %#ok<*NBRAK1> 
    elseif experiments==2
        maxRang=15;
        foldern='blockages-metal-1';
        switchPer=[625]*10^-6;
    elseif experiments==3
        maxRang=15;
        foldern='blockages-cardboard-1';
        switchPer=[625]*10^-6;
    elseif experiments==4
        maxRang=15;
        foldern='blockages-plaster-1';
        switchPer=[625]*10^-6;
    elseif experiments==5
        maxRang=15;
        foldern='blockages-metal-1';
    elseif experiments==6
        maxRang=15;
        foldern='blockages-posterboard-1';
    elseif experiments==7
        maxRang=15;
        foldern='blockages-posterboard-2';
    elseif experiments==8
        maxRang=15;
        foldern='blockages-posterboard-3';
    elseif experiments==9
        maxRang=15;
        foldern='blockages-wood-1';
    elseif experiments==10
        maxRang=15;
        foldern='vspeed-10-1';
    elseif experiments==12
        maxRang=15;
        foldern='vspeed-10-3';
    elseif experiments==13
        maxRang=26;
        foldern='vspeed-20-1';
    elseif experiments==14
        maxRang=26;
        foldern='vspeed-20-2';
    elseif experiments==15
        maxRang=26;
        foldern='vspeed-20-3';
    elseif experiments==17
        maxRang=36;
        foldern='vspeed-30-1';
    elseif experiments==18
        maxRang=36;
        foldern='vspeed-30-2';
    elseif experiments==19
        maxRang=36;
        foldern='vspeed-30-3';
    elseif experiments==20
        maxRang=15;
        foldern='refl-diff-bin';
    elseif experiments==21
        maxRang=15;
        foldern='refl-same-bin';
    elseif experiments==22
        maxRang=15;
        foldern='angles-0-1';
    elseif experiments==23
        maxRang=15;
        foldern='angles-10-1';
    elseif experiments==24
        maxRang=15;
        foldern='angles-10-2';
    elseif experiments==25
        maxRang=15;
        foldern='angles-20-1';
    elseif experiments==26
        maxRang=15;
        foldern='angles-20-2';
    elseif experiments==27
        maxRang=15;
        foldern='angles-30-1';
    elseif experiments==28
        maxRang=15;
        foldern='angles-40-1';
    elseif experiments==29
        maxRang=15;
        foldern='varyspeed-10';
    elseif experiments==30
        maxRang=25;
        foldern='varyspeed-20';
    elseif experiments==31
        maxRang=35;
        foldern='varyspeed-30';
    elseif experiments==33
        maxRang=10;
        foldern='radar1-2';
    elseif experiments==34
        maxRang=10;
        foldern='radar2-2';
    elseif experiments==35
        maxRang=10;
        foldern='radar1-3';
    elseif experiments==36
        maxRang=10;
        foldern='radar2-3';
    elseif experiments==37
        maxRang=10;
        foldern='radar1-4';
    elseif experiments==38
        maxRang=10;
        foldern='radar2-4';
    elseif experiments==39
        maxRang=10;
        foldern='radar2-5';
    elseif experiments==41
        maxRang=10;
        foldern='radar1-5';
    elseif experiments==42
        maxRang=12;
        foldern='doubleradar';
        loadGT=0;
    elseif experiments==43
        maxRang=10;
        foldern='bandf';
        loadGT=0;
    elseif experiments==44
        maxRang=15;
        foldern='test2-uiuc';
        loadGT=0;
    elseif experiments==45
        maxRang=10;
        foldern='220628-154158test-uiuc-bridge';
        loadGT=0;
    else
        skip=1;
    end
    %ADD FILENAME RANGE LOADGT MODF HERE
    %Create new if statement for
    fileName=dir([filePre foldern '/*.mat']);

    startFrame=1;
    modDuty=50;
    switchPer=unique(switchPer);

    selRAll = []; rangeAll = []; velocityAll= []; corrAll = []; angleAll=[];
    
    if(loadGT)
        fileMax = size(fileName,1);
        GTAll = load([filePre foldern '/' fileName(fileMax).name]);
    end
    
    unwrappedDistance2=[];

    for jj=1:size(fileName,1)-1%-8:-1:1 %Change file index
        realFFT = []; imagFFT = []; absFFT = []; angleFFT = []; allFFT=[];
        tagImagAll = []; tagRealAll = []; tagMagAll = []; tagAngAll = []; tagAllAll = [];
        %disp("Running new file: ");
        %disp(jj)
        if(1) %Basic parameter retrieval and creation (Keep collapsed unless modifying)
            fileInd=jj;
            fileName(fileInd).name;
            RadData=load([filePre foldern '/' fileName(fileInd).name]);
            DataAll=RadData.Data;
            RadData.dtime.TimeZone='America/New_York';
            numFrame=size(DataAll,2);%For jj%size(DataAll,2);
            frameTime=RadData.dtime;
            %timeAll(jj, :) = frameTime;

            Cfg=RadData.Cfg;
            if(length(Cfg.Seq) == 1)
                num_tx = 1;
            else
                num_tx = 2;
            end

            N=floor(RadData.N); %number of samples
            fs=RadData.fs; %sampling freq
            num_rx = RadData.NrChn;
            NrChn = num_rx * num_tx; %number of receivers
            CalDat=RadData.CalDat;%ones(4,1);%
            DataAll(1:N:(size(DataAll,1)), :, :)=[]; %remove the chirp number
            antenna_dist=0.006217;%0.00612;

            %sincFunBased
            modF=1./(2*[switchPer]);%1./(2*[3.5e-3]);%1./(2*[1.5e-3]);%1./(2*[500e-6:200e-6:2e-3]);

            % Processing of range profile
            Win2D           =   repmat(Brd.hanning(N-1),1,Cfg.FrmMeasSiz,NrChn);
            ScaWin          =   sum(Win2D(:,1,1));
            NFFT            =   2^12;%RadData.N-1;% 10

            NFFTVel         =   2^8;%size(RadData.Data,1)/RadData.N;%

            kf              =   (Cfg.fStop - Cfg.fStrt)/Cfg.TRampUp;%Cfg.Perd;%
            fc              =   (Cfg.fStop + Cfg.fStrt)/2;
            lambda          =   c/fc;
            vRange          =   [0:NFFT-1].'./NFFT.*fs.*c0/(2.*kf);

            RMin            =   1.5;%0.1;
            RMax            =   min(maxRang,((Cfg.N/Cfg.TRampUp)*c0)/(4*(250e6/Cfg.TRampUp)));

            [~, RMinIdx]   =   min(abs(vRange - RMin));
            [~, RMaxIdx]   =   min(abs(vRange - RMax));
            vRangeExt       =   vRange(RMinIdx:RMaxIdx);
            rOffset         =   max(length(vRangeExt(vRangeExt<RMin)), 1); %Range offset

            doppSize        =   Cfg.FrmMeasSiz;%min(128,Cfg.FrmMeasSiz);
            WinVel          =   Brd.hanning(doppSize);
            ScaWinVel       =   sum(WinVel);
            WinVel2D        =   repmat(WinVel.',numel(vRangeExt),1);

            %vFreqVel        =   [-NFFTVel./2:NFFTVel./2-1].'./NFFTVel.*(1/Cfg.Perd);
            %vVel            =   vFreqVel*c0/(2.*fc);

            % Window function for receive channelsf
            NFFTAnt         =   256;
            WinAnt          =   Brd.hanning(NrChn);
            ScaWinAnt       =   sum(WinAnt);
            WinAnt2D        =   permute(repmat(WinAnt,1,Cfg.FrmMeasSiz,numel(vRangeExt)),[3,2,1]);
            %WinAnt2D        =   repmat(WinAnt.',numel(vRangeExt)*Cfg.FrmMeasSiz,1);
            vAngDeg         =   asin(2*[-NFFTAnt/2:NFFTAnt/2-1].'./NFFTAnt)./pi*180;


            % Calibration data
            mCalData        =   permute(repmat(CalDat(1:NrChn),1,Cfg.FrmMeasSiz,N-1),[3 2 1]);

            % Positions for MeasIdxpolar plot of cfost function
            % vU              =   linspace(-1,1,NFFTAnt);
            % [mRange , mU]   =   ndgrid(vRangeExt,vU);
            % mX              =   mRange.*mU;
            % mY              =   mRange.*cos(asin(mU));

            tagSideWidth = 1.0; %Length of one side of a tag's signature in a range domain in m(approximate for now)
            vRangeWidth  = round(tagSideWidth.*NFFT./fs./c0*(2.*kf) * 125/N); %Convert above measurement into a number of bins
        end

        RPExtAll=zeros(numFrame,numel(vRangeExt),Cfg.FrmMeasSiz,NrChn);
        RDAll=zeros(numFrame,numel(vRangeExt),NFFTVel,NrChn);
    
        %% Initial Processing and FFTs
        for MeasIdx = 1:1:numFrame-startFrame+1
            Data        =  reshape(squeeze(DataAll(:,MeasIdx,:)),N-1,[], NrChn);

            if(bkndSubtract)
                Data=Data-Data(:,1,:); %#ok<UNRCH> 
            end

            % Calculate range profile including calibration
            RP          =   2*fft(Data.*Win2D.*mCalData,NFFT,1).*Brd.FuSca/ScaWin; %.*Win2D
            RPExt       =   RP(RMinIdx:RMaxIdx,:,:);
            RD          =   fft(RPExt.*WinVel2D, NFFTVel, 2)./ScaWinVel;
            RDA=RD;
            RPExtAll(MeasIdx-startFrame+1,:,:,:)=RPExt;
            RDAll(MeasIdx-startFrame+1,:,:,:)=RDA;
        end

        %% Find Range Bins/Do Matched Filtering       
        for MeasIdx = 1:numFrame-1+startFrame
            [selR, range, velocity, corr] = matched_filtering(Cfg,modF,modDuty,NFFTVel,vRangeExt, fc, c,fs, numFrame, startFrame, squeeze(RDAll(MeasIdx,:,:,:)));
            selRAll(jj, MeasIdx-startFrame+1, :) = selR; %#ok<*AGROW> 
            rangeAll(jj, MeasIdx-startFrame+1, :) = range;
            velocityAll(jj, MeasIdx-startFrame+1, :) = velocity;
            corrAll(jj, MeasIdx-startFrame+1, :) = corr;
        end

        selROverall(jj, :) = round(median(selRAll(jj, :, :), 2));
        rangeOverall(jj, :) = median(rangeAll(jj, :, :), 2);
        velocityOverall(jj, :) = mean(velocityAll(jj, :, :), 2);
        %% Find Angle from Radar
        % function [aoa_f, aoa_ind]=calculate_angle(Cfg, modF(tagNum), modDuty, NFFTVel, vRangeExt, fc, c, fs, velocities(tagNum))
        for MeasIdx = 1:1:numFrame-startFrame+1
            %disp(MeasIdx);
            [aoa_f, aoa_ind] = estimate_aoa(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc, c, fs, velocity, WinAnt2D, NFFTAnt, ScaWinAnt, vAngDeg, squeeze(RDAll(MeasIdx, selROverall(jj, :), :, :)));

            angleAll(jj,MeasIdx-startFrame+1,:)=aoa_f;
            angleIndAll(jj,MeasIdx-startFrame+1,:)=aoa_ind;
        end

        angleOverall(jj, :) = mean(angleAll(jj, :, :), 2);
        selAngleOverall(jj, :) = round(mean(angleIndAll(jj, :, :), 2));

        %% Calculate Tag FFTs for Phase Extraction
        for tagNum = 1:length(modF)
            expDopFinals=[];

            expDopFinals(:,:)=calculateFilter_wMobility(Cfg, modF(tagNum), modDuty, NFFTVel, vRangeExt, fc(1), c, fs, velocityOverall(jj, tagNum));

            [~, indxAll] = max(expDopFinals(1, 1:floor(NFFTVel/2), :), [], 2);
            indx = squeeze(indxAll);
            %             [~, indxAll2] = max(expDopFinals(1, (floor(NFFTVel/2):end), :), [], 2);
            %             indxAll2 = squeeze(indxAll2);
            for MeasIdx = 1:1:numFrame-startFrame+1
                %disp(MeasIdx);
                indsAll = selROverall(jj,tagNum) + (-vRangeWidth:vRangeWidth);
                indsUse = indsAll(indsAll >= 1 & indsAll <= size(RPExtAll, 2));
                RPExt = squeeze(RPExtAll(MeasIdx-startFrame+1, indsUse, :,  :));
                RPExt = permute(RPExt, [2 1 3]);

                realFFT(:, :, :) = fft(real(RPExt), NFFTVel, 1);
                imagFFT(:, :, :) = fft(imag(RPExt), NFFTVel, 1);
                absFFT(:, :, :) = fft(abs(RPExt), NFFTVel, 1);
                angleFFT(:, :, :) = fft(unwrap(angle(RPExt)), NFFTVel, 1);
                allFFT(:, :, :) = fft((RPExt), NFFTVel, 1);


                %max(imagFFT(max(indx-3, 1):min(indx+3, end), :)
                tagImagAll(jj, MeasIdx-startFrame+1, tagNum, :, :) = squeeze(max(imagFFT(max(indx-5, 1):min(indx+5, end), :, :), [], 1));
                tagRealAll(jj, MeasIdx-startFrame+1, tagNum, :, :) = squeeze(max(realFFT(max(indx-5, 1):min(indx+5, end), :, :), [], 1));
                tagMagAll(jj, MeasIdx-startFrame+1, tagNum, :, :)  = squeeze(max(absFFT (max(indx-5, 1):min(indx+5, end), :, :), [], 1));
                %                         tagImagAll(jj, MeasIdx-startFrame+1, sel, indx, tagNum, :) = squeeze(imagFFT(indx, :));
                %                         tagRealAll(jj, MeasIdx-startFrame+1, sel, indx, tagNum, :) = squeeze(realFFT(indx, :));
                %                         tagMagAll(jj, MeasIdx-startFrame+1, sel, indx, tagNum, :) = squeeze(absFFT(indx, :));
                tagAngAll(jj, MeasIdx-startFrame+1, tagNum, :, :) = squeeze(max(angleFFT(max(indx-5, 1):min(indx+5, end), :, :), [], 1));
                tagAllAll(jj, MeasIdx-startFrame+1, tagNum, :, :) = squeeze(max(allFFT(max(indx-5, 1):min(indx+5, end), :, :), [], 1));
                %                 all1(tagNum, MeasIdx-startFrame+1, :) = squeeze(max(allFFT( max(indx-3, 1): min(indx+3, end), :), [], 1));
                %                 all2(tagNum, MeasIdx-startFrame+1, :) = squeeze(max(allFFT(max(indx2-3, 1):min(indx2+3, end), :), [], 1));
            end
        end

        %% Extract Tag Phase from FFTs - METHOD 2
        for(tagNum = 1:length(modF))
            measIndices = 1:numFrame-startFrame+1;
            realFFTAll = squeeze((tagRealAll(jj, measIndices, tagNum, 1:size(tagRealAll, 4), :)));
            imagFFTAll = squeeze((tagImagAll(jj, measIndices, tagNum, 1:size(tagRealAll, 4), :)));
            magFFTAll  = squeeze((tagAllAll(jj, measIndices, tagNum, 1:size(tagRealAll, 4), :)));

            scaledReal = abs(realFFTAll)./abs(magFFTAll);
            scaledImag = abs(imagFFTAll)./abs(magFFTAll);

            phi=[0:pi/1000:pi]; %43 = size of distance btwn two peaks
            widths= pi*[1, 1/10, 1/100, 1/1000, 1/10000];
            x=1:size(magFFTAll, 2);

            %Quick and dirty frequency extraction
            freqIndex=[];

            indind = 1;
            for(MeasIdx=1:(numFrame-startFrame+1)/10:numFrame-startFrame+1) %select 10 represenative samples
                for(antNumber=2:3) %across the two central antennas
                    [~, I1] = findpeaks(scaledReal(MeasIdx,:,antNumber), 'MinPeakProminence', 0.5);
                    [~, I2] = findpeaks(scaledImag(MeasIdx,:,antNumber), 'MinPeakProminence', 0.5);
                    per1temp = diff(I1);
                    per2temp = diff(I2);
                    l1=length(per1temp); l2=length(per2temp);
                    freqIndex = [freqIndex, per1temp, per2temp];
                end
            end
            period = trimmean(freqIndex, 10);

            tri1 = repmat(triang(length(x))', [size(scaledReal, 1)  1 4]);
            %for(MeasIdx = 1:250)%numFrame-startFrame+1)
            autoCorr2a = [];
            autoCorr2b = [];
            for(phiInd=1:length(phi))
                sinW = repmat(abs(sin(x*pi/period + phi(phiInd))), [size(scaledReal, 1) 1 4]);
                cosW = repmat(abs(cos(x*pi/period + phi(phiInd))), [size(scaledReal, 1) 1 4]);

                autoCorr2a(:, phiInd, :) = (sum(((scaledReal(:, :, :))) .* sinW .* tri1, 2));
                autoCorr2b(:, phiInd, :) = (sum(((scaledImag(:, :, :))) .* cosW .* tri1, 2));
            end
            %end
            autocorr = autoCorr2a(:, :, :)+autoCorr2b(:, :, :);
            [~, I]   = max(autocorr, [], 2);
            angles   = squeeze(phi(I));
            %Do some comparison of angle as well?
            unwrappedAngle2(jj, tagNum, 1:size(angles, 1), :) = unwrapPi(angles);
            unwrappedDistance2(jj, tagNum, :, :) = 1*lambda/(4*pi)*unwrappedAngle2(jj, tagNum, :, :);
        end

        %% Compare displacements to GT
        coherenceFrames=[15];
        numSamples = 50;

        if(loadGT)
            for tagNum = 1:length(modF)
                for sample = 1:numSamples
                    randindsStart = randperm(start_and_end_delay,round(start_and_end_delay/2))+5;
                    randindsEnd = numFrame + 1 - randperm(start_and_end_delay,round(start_and_end_delay/2));

                    pp    = (squeeze(unwrappedDistance2(jj, tagNum, :, :)));
                    dist2 = squeeze(mean(trimmean(pp(randindsEnd, :), 10) - mean(pp(randindsStart, :), 1), 2));

                    %GTErrors1 or 2 depends on the direction of your motion stage - I'm not sure how your setup works.
                    GTErrors1(jj, tagNum, sample) = dist2 + 1/100*GTAll.posn(jj); %May need to change scaling based on input in gtdists files
                    GTErrors2(jj, tagNum, sample) = dist2 - 1/100*GTAll.posn(jj);
                    distances(jj, tagNum, sample) = dist2; %*3.7778; %Some cnst scale
                    gtdists1(jj, tagNum, sample)    = GTAll.posn(jj);
                end
                %end
            end
        
        end

        if(0)
            figure(experiments) %#ok<UNRCH> 
            if(jj==1)
                title(["Experiment ", ""+experiments])
            end
            subplot(5,4,jj)
            plot(squeeze(unwrappedDistance2(jj, 1, :, 1:3)))
            title(["jj =  "+jj+" : GT "+(1/100*GTAll.posn(jj))])
        end
    end
    
    %plots unwrapped distances to verify correctness

    if(0)
        for(ind=1:jj)
            figure(1) %#ok<UNRCH> 
            if(jj==1)
                title(["Experiment "+experiments+" diagnostics unwrap plots"])
            end
            subplot(5,4,jj)
            plot(squeeze(unwrappedDistance2(jj, 1, :, 1:3)))
            title(["jj =  "+jj+" : GT "+(1/100*GTAll.posn(jj))])
        end
    end
clear Brd;
end

%% Filter Function Generator
function expDopFinal=calculateFilter_wMobility(Cfg,modF,modDuty,NFFTVel,vRangeExt,fc,c,fs,vel)
t=[0:Cfg.Perd:Cfg.FrmMeasSiz*Cfg.Perd];%oversampling to get rid of the noise
expDopFinal=zeros(length(vRangeExt),NFFTVel);

doppShift=fft(cos(2*pi*(2*vel*fc/c).*t), NFFTVel,2);
doppShift=doppShift(1:NFFTVel/2) ;
doppShift = doppShift/max(doppShift);
doppShift(doppShift < 0.6) = 0;

undopp_square=fft(square(2*pi*(modF.').*t(1:end-1),modDuty),NFFTVel,2);
expdopt=abs(undopp_square);
expdopt=expdopt./max(expdopt);
expdopt(expdopt<0.1)=0;
expdopt=expdopt./sum(expdopt)*100; %Normalize such that each filter will have the same 'amount' of autocorrelation

if(vel>0)
    doppShiftPos=[doppShift zeros(1,NFFTVel/2)];
    sq_wav=conv(doppShiftPos, expdopt);
    expdopt=abs(sq_wav(1:NFFTVel));
elseif(vel<0)
    doppShiftNeg=flip([doppShift  zeros(1,NFFTVel/2)]);
    sq_wav=conv(doppShiftNeg, expdopt);
    expdopt=abs(sq_wav(NFFTVel:end));
else
    expdopt=abs(expdopt);
end

expdopt=expdopt./max(expdopt);
expdopt(expdopt<0.1)=0;

expDopFinal=repmat(expdopt, length(vRangeExt), 1);
end

%{
Calculates the polar median of complex data. Instead of using purely real
and imaginary data, it converts to magnitude and phase, calculates the
median in the polar domain, and then reconverts to complex numbers. 
%}

%% Polar Median
function med=polarMedian(input, dimension)
if(nargin<2)
    dimension=1;
end
a=abs(input);
b=angle(input);
a2=median(a, dimension, 'omitnan');
b2=median(b, dimension, 'omitnan');
[realPart, imagPart]=pol2cart(b2, a2);
med=squeeze(complex(realPart, imagPart));
end
%% AoA Estimate
function aoa_estimate = temp_f() %(solves optimization problem for angles)
phi=angdiff(repmat(angle(tagMeans(end,1)), 1, size(tagMeans, 2)), angle(tagMeans(end,:)));
d=vRangeExt(selROverall(jj, 1)-rOffset);
phi=phi*lambda/(4*pi);
A=[2*phi(2), -lambda; 2*phi(3), -2*lambda; 2*phi(4), -3*lambda];
b=[lambda^2/4 - 2*d*phi(2) - phi(2)^2; lambda^2 - 2*d*phi(3) - phi(3)^2; 9*lambda^2/4 - 2*d*phi(4) - phi(4)^2];
[Q, R]=qr(A, 0);
soln=inv(R)*(Q' * b);
soln(1) = soln(1) + d; %Check if equations are wrong, seems like this correction is needed
vertD=soln(2);
unknE = soln(1);
angleEst = asin( repmat([vertD], 4, 1)./(d + unknE + [0; phi(2); phi(3); phi(4)]) );
aoa_estimate = median(angleEst);
end

%State2 is imaginary
%% Corrector
function correctAngle = corrector(state1, state2, angleIn)
angleOut=[];
for i=1:length(state1)
    if(state1(i) && state2(i))
        angleOut=pi/2 - angleIn;
    elseif(state1(i))
        angleOut=pi/2 + angleIn;
    elseif (state2(i))
        angleOut=-pi/2 - angleIn;
    else
        angleOut=-pi/2 + angleIn;
    end
end
end

%% Find GT2
function idx=findGT2(frameTime,GT_time)
[h1, m1, s1] = hms(frameTime);
[h2, m2, s2] = hms(GT_time);
[val,idx]=min(abs(m1*60 + s1 - m2*60 - s2));
if val>=1
    idx=-1;
end
end

%% Find GT
function idx=findGT(frameTime,GT_time)
[val,idx]=min(abs(frameTime-GT_time));
if val>=seconds(1)
    idx=-1;
end
end
% function ud = correctPiJumps(ud2)
%     ud= [];
%     for(i = 1:size(ud2, 1))
%         udx = squeeze(ud2(i, :, :));
%         for(kk=1:size(udx,  3))
%            for(jj = 2:size(udx, 2))
%
%            end
%         end
%     end
% end

%% Extract Angle
function angles_all = extract_angle(unwrappedAngle, rangeDetected, antenna_dist, lambda)
angles_all=[];
for(MeasIdx=1:size(unwrappedAngle,1))
    vrange = round(rangeDetected/lambda)*lambda;

    phases = mod(squeeze(unwrappedAngle(MeasIdx, :)), 2*pi);
    x_r= [0 antenna_dist 2*antenna_dist 3*antenna_dist];
    y_r= [0 0 0 0];
    centers = [x_r', y_r'];
    dist_offset = phases/(4*pi)*lambda;
    radii = dist_offset + rangeDetected;

    angles_temp=[];
    for(i=1:4)
        for(j=1:4)
            if(i ~= j)
                [x_out, y_out] = circcirc(x_r(i), y_r(i), radii(i), x_r(j), y_r(j), radii(j));
                angle = atan(y_out(1)/x_out(1));
                angles_temp = cat(1, angles_all, angle);
            end
        end
    end
    angles_all = cat(1, angles_all, median(angles_temp));
end
end

%% Angle from Unwrap
function angleFromPhase = angFromUnwrap(unwrappedAngle, jj, tagNum)
a=diff(squeeze(unwrappedAngle(jj, tagNum, :, :))')';
angleFromPhase = asin(unwrapPi(wrapToPi2(a(3:end,:)))/pi)*180/pi;
end

%% Differential Angle
function new_theta = differentialAngle(init_ang, dPhase)
old_diff = pi * sin(init_ang*pi/180);
new_theta = asin((old_diff + dPhase) / pi) *180/pi;
end

%% Phase Match Filter
function phase_match_filter = output(rang_mag, IQ_period_est)
phaseOuterBins = 1:180;
for(i = 1:length(phaseOuterBins))
    T=14;
    filter_f = abs(sin(2*pi / T * 1:length(rang_mag) + phaseOuterBins(i)));
    correlation = rang_mag * zscore(filter_f);
    phaseOuterEst(phaseOuter) = correlation;
end
[~, ind] = max(phaseOuterEst);
phaseOuter = phaseOuterBins(ind);

phaseInnerBins = -1.5:0.01:1.5;
for(i = 1:length(phaseOuterBins))
    T=14;
    filter_f = abs(sin(2*pi/T + phaseInnerBins(i) * 1:length(rang_mag) + phaseOuter + phaseInnerBins(i)));
    correlation = rang_mag * zscore(filter_f);
    phaseInnerEst(i) = correlation;
end
[~, ind] = max(phaseInnerEst);
phaseInner = phaseOuterBins(ind);
output = phaseOuter+phaseInner;
end

%% Plot RDA
function plotRDA(RPExtAll, MeasIdx, ind)
RDA1 = fft(imag(squeeze(RPExtAll(:, :, MeasIdx, 1))), [], 2); ind= 36;
RDA2 = fft(real(squeeze(RPExtAll(:, :, MeasIdx, 1))), [], 2);
RDA3 = fft((squeeze(RPExtAll(:, :, MeasIdx, 1))), [], 2);
hold off;
a1 = abs(((RDA1(:, ind))));
a2 = abs(((RDA2(:, ind))));
a3 = abs(((RDA3(:, ind))));
figure
plot(a1); hold on
plot(a2)
plot(a3)
p1 = angle(((RDA1(:, ind))));
p2 = angle(((RDA2(:, ind))));
p3 = angle(((RDA3(:, ind))));
figure
plot(p1); hold on;
plot(p2);
plot(p3);
end

%% Helper Function 1
function helper1(normSig)
[~,b]=max(normSig(:, 6:end-6), [], 'all', 'linear');
[maxR, maxC]   = ind2sub(size(normSig(:, 6:end-6)), b);
nPks = findpeaks(-normSig(:,6+maxC-1));
[~, inds] = findpeaks(-normSig(:,6+maxC-1));
allInds = inds;
i1=find(allInds > maxR, 1, 'first');
i2=find(allInds < maxR, 1, 'last');
if(isempty(i1))
    i1=0;
end
if(isempty(i2))
    i2=size(normSig, 1);
end
[maxR-allInds(i2), allInds(i1)-maxR]
end

%% Matched Filtering
function [selR, range, velocity, corr] = matched_filtering(Cfg,modF,modDuty,NFFTVel,vRangeExt, fc, c,fs, numFrame, startFrame, RDA)
% define the signal for correlation matching
sig=squeeze(abs(RDA(:,:,:)));%squeeze(max(RDA(:,:,:),[],3)); %get the max value across 4 channels

% normalize the signal
normA = sig - min(sig(:));
normSig = normA ./ max(normA(:));
normSig = sum(normSig(:,:,1:3), 3);

%Calculate background noise across frequency bins
doppNoise=sum(sum(normSig(:,:,:),3),1);
convW=17;
doppNoise=conv(doppNoise, triang(convW));
doppNoise=doppNoise(ceil(convW/2):end-floor(convW/2));
doppNoise=normalize(doppNoise,2,'range');

%default uses values 10 and -0.2
%mobi using values 400*(-0.03)
sigmoid(:) = 1./(1+exp(10*(-0.2+doppNoise)));

vel=-0.2:0.05:0.2; %for no mobility cases

normsumR2d=zeros(length(modF), length(vel), length(vRangeExt));

%Sweep each filter across results
for vind=1:length(vel)
    for k=1:length(modF)
        expDopFinal=calculateFilter_wMobility(Cfg,modF(k),modDuty,NFFTVel,vRangeExt, fc, c,fs,vel(vind));
        expDopFinal(isnan(expDopFinal)) = 0;
        expDopFinal = expDopFinal ./ (sum(expDopFinal(1, :))) *3; %Want to bias away from a matched filter with more peaks/just larger
        expDopSigmoid=expDopFinal.*sigmoid;%.*sigmoid(MeasIdx-startFrame+1, :);
        normexpDopFinal=expDopSigmoid;

        normR2=squeeze(sum(abs(normSig).*repmat(abs(normexpDopFinal),1,1,size(normSig,3)),2));
        normsumR2d(k, vind, :)=normR2;
    end
end

newFinal=zeros(length(modF), 1);
velocities=zeros(length(modF), 1);
newSelR=zeros(length(modF), 1);

%Find maximum value and record

for tagNum=1:length(modF)
    %rOffset=120; %long range data
    [peakInds, newFinal(tagNum)]=findPeaks2D(squeeze(normsumR2d(tagNum,:,1:end)), 1);
    velocities(tagNum)=vel(peakInds(1));
    newSelR(tagNum)=peakInds(2);
end

selR = newSelR;
range = vRangeExt(newSelR);
velocity = velocities;
corr = newFinal;
end

%% Estimate_AoA
function [aoa_f, aoa_ind]=estimate_aoa(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc, c, fs, velocity, WinAnt2D, NFFTAnt, ScaWinAnt, vAngDeg, RDABin)
for tagNum = 1:length(modF)
    expDopFinal=[];
    expDopFinal(:,:)=calculateFilter_wMobility(Cfg, modF(tagNum), modDuty, NFFTVel, vRangeExt, fc, c, fs, velocity(tagNum));

    %apply to fft across antennas
    RDAngle = fftshift(fft(RDABin.*squeeze(WinAnt2D(1,1,:))', NFFTAnt, 2)/ScaWinAnt,2);
    rangeBin(:,:)=squeeze(abs(RDAngle));
    normexpDopFinal=expDopFinal(1,:);

    %Find maximum correlation
    normR2angle=squeeze(sum(rangeBin.*normexpDopFinal', 1));
    [~, aoa_ind]=max(normR2angle);

    aoa_f=vAngDeg(aoa_ind);
end
end