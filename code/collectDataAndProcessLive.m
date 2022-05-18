% AN24_02 -- FMCW Basics 
function collectDataAndProcessLive()
maxRange = 10;
chirpDur=125; 
c = 3e8;
%***************
%filename='test';
% (1) Connect to TinyRad: Check if Brd exists: Problem with USB driver
% (3) Configure RX
% (4) Configure TX
% (5) Start Measurements
% (6) Configure calculation of range profile
% Configure script
%--------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
CurPath = pwd();
addpath([CurPath,'/DemoRadBasic/DemoRadUsb']);
addpath([CurPath,'/DemoRadBasic/Class']);

%--------------------------------------------------------------------------
% Setup Connectionr
%--------------------------------------------------------------------------
Brd         =   TinyRad();

% Reset prvious configuration if board was already configured
Brd.BrdRst();

%--------------------------------------------------------------------------
% Software Version
%--------------------------------------------------------------------------
Brd.BrdDispSwVers();

%--------------------------------------------------------------------------
% Configure Receiver
%--------------------------------------------------------------------------
Brd.RfRxEna();


%--------------------------------------------------------------------------
% Configure Transmitter (Antenna 0 - 2, Pwr 0 - 100)
%--------------------------------------------------------------------------
TxPwr           =   100;
Brd.RfTxEna(1, TxPwr);
CalDat          =   Brd.BrdGetCalDat();
CalDat          =   CalDat(1:4);

%--------------------------------------------------------------------------
% Configure Up-Chirp and timing for the measurements
%--------------------------------------------------------------------------
Cfg.fStrt       =   24.0e9;                    %   Start frequency   
Cfg.fStop       =   24.25e9;                    %   Stop frequency
Cfg.TRampUp     =   chirpDur*10^-6;                     %   UpChirp duration 
Cfg.Perd        =   Cfg.TRampUp+100e-6;                     %   Period between measurements
Cfg.N           =   100;  %   Number of samples taken at start of chirp 
Cfg.Seq         =   [1];                        %   Antenna transmit sequence
Cfg.CycSiz      =   2;                          %   Number of buffers in the acquisition framework 2 are required
Cfg.FrmSiz      =   256;                          %   Number of chirp sequences for one measurement cycle
Cfg.FrmMeasSiz  =   256;                          %   Number of chirps sequences for collecting IF data
Brd.RfGet('kf')

Brd.RfMeas(Cfg);

disp('1')
%--------------------------------------------------------------------------
% Read actual configuration
%--------------------------------------------------------------------------
N               =   Brd.Get('N');
NrChn           =   Brd.Get('NrChn');
fs              =   Brd.Get('fs');

%--------------------------------------------------------------------------
% Check TCP/IP data rate:
% 16 Bit * Number of Enabled Channels * Number of Samples are measureed in
% the interval TInt. If the data rate is too high, than frames can be losed
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Configure Signal Processing
%--------------------------------------------------------------------------
disp('Get Measurement data')
measRounds=1000;%700 for moving traces %300 for regular stationary ones   %floor(7/(Cfg.FrmMeasSiz*Cfg.Perd));

drawGraphs = 1;
drawHeatmap = 0;
connectAR = 0;

c0 = 1/sqrt(8.85e-12*4.*pi.*1e-7);
%--------------------------------------------------------------------------
% Include all nefdracessary directories
%--------------------------------------------------------------------------

%Brd         =   TinyRad();
minRD=1000;
maxRD=-1000;
minR=1000;
maxR=-1000;
minAng=1000;
maxAng=-1000;

bknSubtract=0;
filePre='./Data/';
c=3e8;

if connectAR == 1
    %Modify IP and port here
    %Create a matlab tcp server to send data
    clientAR = tcpip("127.0.0.1", 5311, 'NetworkRole', 'client');  %192.168.1.2
    fopen(clientAR);
end

maxRang=maxRange;

 %*******************
modDuty=50;
vel=[0];

rangeAll=[];angleAll=[];xAll=[];yAll=[];
rangeAllMUSIC=[];angleAllMUSIC=[];xAllMUSIC=[];yAllMUSIC=[];
rangeAllKAL=[];angleAllKAL=[];


tagHeatmap=[];

antenna_dist=0.00612;

%sincFunBased2
%TODO
switchPer=[625]*10^-6;%[500]*10^-6; %625 1000%1./(2*[3.5e-3]);%1./(2*[1.5e-3]);%1./(2*[500e-6:200e-6:2e-3]);
switchPer=unique(switchPer);
modF=[1./(2*[switchPer])]; %1000

% Processing of range profile
Win2D           =   repmat(Brd.hanning(N-1),1,Cfg.FrmMeasSiz,NrChn);
ScaWin          =   sum(Win2D(:,1,1));
NFFT            =   2^10;%RadData.N-1;%

NFFTVel         =   2^8;%size(RadData.Data,1)/RadData.N;%
kf              =   (Cfg.fStop - Cfg.fStrt)/Cfg.TRampUp;%Cfg.Perd;%
vRange          =   [0:NFFT-1].'./NFFT.*fs.*c0/(2.*kf);

fc              =   (Cfg.fStop + Cfg.fStrt)/2;
lambda= c/fc;

RMin            =   0.2;
RMax            =   min(maxRang,((Cfg.N/Cfg.TRampUp)*c0)/(4*(250e6/Cfg.TRampUp)));

[Val RMinIdx]   =   min(abs(vRange - RMin));
[Val RMaxIdx]   =   min(abs(vRange - RMax));
vRangeExt       =   vRange(RMinIdx:RMaxIdx);
rOffset         =   length(vRangeExt(vRangeExt<RMin)); %Range offset

doppSize        =   Cfg.FrmMeasSiz;%min(128,Cfg.FrmMeasSiz);
WinVel          =   Brd.hanning(doppSize);
ScaWinVel       =   sum(WinVel);
WinVel2D        =   repmat(WinVel.',numel(vRangeExt),1);

vFreqVel        =   [-NFFTVel./2:NFFTVel./2-1].'./NFFTVel.*(1/Cfg.Perd);
vVel            =   vFreqVel*c0/(2.*fc);

% Window function for receive channels
NFFTAnt         =   127;
WinAnt          =   Brd.hanning(NrChn);
ScaWinAnt       =   sum(WinAnt);
WinAnt2D        =   permute(repmat(WinAnt,1,Cfg.FrmMeasSiz,numel(vRangeExt)),[3,2,1]);
%WinAnt2D        =   repmat(WinAnt.',numel(vRangeExt)*Cfg.FrmMeasSiz,1);
vAngDeg         =   asin(2*[-NFFTAnt./2:NFFTAnt./2-1].'./NFFTAnt)./pi*180;


% Calibration data
mCalData        =   permute(repmat(CalDat,1,Cfg.FrmMeasSiz,N-1),[3 2 1]);

% Positions for polar plot of cost function
vU              =   linspace(-1,1,NFFTAnt);
[mRange , mU]   =   ndgrid(vRangeExt,vU);
mX              =   mRange.*mU;
mY              =   mRange.*cos(asin(mU));

JNormAll=[];
RPExtAll=zeros(numel(vRangeExt),Cfg.FrmMeasSiz,measRounds,NrChn);
RDAll=zeros(numel(vRangeExt),NFFTVel,measRounds     ,NrChn);       

dtime=[];
DataAll=[];
if(drawHeatmap==1)
        %Initialize heatmap
        tagHeatmap=heatmapSetup(mX, mY, length(modF), 'Live heatmap');
end

DataRate        =   16*NrChn*N.*Cfg.FrmMeasSiz./(Cfg.Perd*Cfg.FrmSiz);
disp(['DataRate: ', num2str(DataRate/1e6), ' MBit/s'])

%PrevState
KalmanXk_prev = [];
% Current state estimate
KalmanXk = [];

KalmanP=[];
for eyeInd = 1:length(switchPer)
    KalmanP(:, :, eyeInd) = eye(4);
end

% Q is the process noise covariance. 
% Should probably find out if these are correlated: TODO
KalmanQ = 1/2*[1 0 0.1 0; 0 1 0 0.1; 0 0 1 0; 0 0 0 1];

% M is the measurement matrix.

        KalmanM = eye(4);

% R is the measurement noise covariance. Varies between samples
%Higher -> more on prediction, lower -> more on measurement
KalmanR = 2*[0.5 0 0 0; 0 1 0 0; 0 0 2 0; 0 0 0 1].^2;
frameTime = 0;
MeasIdx = 0; 
 while 1
    MeasIdx= MeasIdx + 1;
    pause(0.1)
%     freqs = linspace(24.0e9, 24.25e9, 6)
%     Cfg.fStrt       =   freqs(mod(MeasIdx, 5) + 1);                    %   Start frequency   
%     Cfg.fStop       =   freqs(mod(MeasIdx, 5) + 2);                    %   Stop frequency
%     Brd.BrdRst();
%     Brd.RfMeas(Cfg);
    
    disp(['reading frame', num2str(MeasIdx)])
         %remove the chirp number
        % read data for configured measurement sequence
       DataAll(:,1,:)=Brd.BrdGetData();
       %dtime(MeasIdx,:)=datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS');
    tic
        
        
    %--------------------------------------------------------------------------
    % Measure and calculate range, doppler, DBF
    %--------------------------------------------------------------------------
    Data        =  DataAll(:,1,:);
    Data(1:N:size(Data,1),:,:)=[]; 
    Data = reshape(squeeze(Data(:, :, :)),N-1,[],size(Data,3));

    %background subtraction
    Data=Data-Data(:,1,:);

    % Calculate range profile including calibration
    RP          =   2*fft(Data.*Win2D.*mCalData,NFFT,1).*Brd.FuSca/ScaWin;
    RPExt       =   RP(RMinIdx:RMaxIdx,:,:);    
    %RPExtAll(:,:,MeasIdx,:)=RPExt;

    % calculate fourier transform over receive channels
%     JOpt        =   fftshift(fft(RPExt.*WinAnt2D,NFFTAnt,3)/ScaWinAnt,3);
%     % normalize cost function
%     JdB         =   20.*log10(abs(JOpt));
%     JMax        =   max(JdB(:));
%     JNorm       =   JdB - JMax;
%     JNorm(JNorm < -25)  =   -25;    
%     JNorm=(squeeze(max(JNorm,[],2)));
%     JNormAll(MeasIdx,:,:)=(squeeze(max(JdB,[],2)));

    if(drawHeatmap==1)
        tagHeatmap=heatmapBack(tagHeatmap, JNorm);
    end
    
    % calculate range-doppler profile
     RD          =   fft(RPExt.*WinVel2D, NFFTVel, 2)./ScaWinVel;
     RDA=RD;%cat(2,RDA,RD);

%      RDAll(:,:,MeasIdx,:)=RDA;

     % define the signal for correlation matching
      sig=squeeze(abs(RDA(:,:,1:4)));%squeeze(max(RDA(:,:,:),[],3)); %get the max value across 4 channels

      % normalize the signal
      normA = sig - min(sig(:));
      normSig = normA ./ max(normA(:));
    
    %template matching 
            indind=1;
            normsumR2d=[];
            doppNoise=sum(sum(normSig(:,:,:),3),1);
            doppNoise=normalize(doppNoise,2,'range');
            
            %Convolve with triangle to smooth
            convW=17;
            doppNoise=conv(doppNoise, triang(convW));
            doppNoise=doppNoise(ceil(convW/2):end-floor(convW/2));
            doppNoise=normalize(doppNoise,2,'range');
            
            %[~, noiseLoc] = findpeaks(doppNoise, 'NPeaks', 1, 'MinPeakHeight', 0.6);
            %Turns noise into attuning function. Can increase the multiplier
            %(15) to make it steeper or the shift (-0.05) to change the 0.5
            %gain point
            sigmoid=1./(1+exp(200*(-0.16+doppNoise)));
            %gtVel=2.5 <> shift of about 36 bins
%             if(~isempty(noiseLoc)) %Center velocity around dopp frequency
%                 velLoc=noiseLoc;
%                 if(velLoc>NFFTVel/2)
%                     velLoc=velLoc-NFFTVel;
%                 end
%                 velCenter=7.5*(velLoc/(NFFTVel/2)); %(c/fc*modF)? FIX THIS, 7.5 is a place holder, not relient on modF
%                 velUncertainty=1;
%                 %approx one bin per 5 cm/s
%                 vel=linspace(velCenter-velUncertainty, velCenter+velUncertainty, 60);
%             else
%                 vel=linspace(-1, 1, 60);
%             end
            
            vel=[-3.6:0.05:3.6];
            
            for vind=1:length(vel)
                for k=1:length(modF)
                    expDopFinal=calculateFilter_wMobility(Cfg,modF(k),modDuty,NFFTVel,vRangeExt, fc, c,fs,vel(vind));
                    expDopFinal(isnan(expDopFinal)) = 0;
                    expDopFinal = expDopFinal ./ (sum(expDopFinal(1, :))); %Want to bias away from a matched filter with more peaks/just larger
                    expDopSigmoid=expDopFinal.*sigmoid;
                    normexpDopFinal=expDopSigmoid;
                    
                    normR2=squeeze(sum(abs(normSig).*repmat(abs(normexpDopFinal),1,1,size(normSig,3)),2));
                    normsumR2d(k, vind, :)=sum(normR2, 2);
                end
                indind=indind+1;
            end
            
            newFinal=zeros(length(modF), 1);
            velocities=zeros(length(modF), 1);
            newSelR=zeros(length(modF), 1);
    for tagNum=1:length(modF)
        %rOffset=120; %long range data
        [peakInds, newFinal(tagNum) ]=findPeaks2D(squeeze(normsumR2d(tagNum,:,1:end)), 1);
        velocities(tagNum)=vel(peakInds(1));
        newSelR(tagNum)=peakInds(2);
    end
    finalVal=newFinal;
    selR=newSelR;
    % extract the tag parameters
     threshTagDet=.01;
     for tagNum=1:length(modF)
         if finalVal(tagNum)>threshTagDet
            range_f=vRangeExt(min(length(vRangeExt), selR(tagNum)+rOffset)) + 0.3; %TODO
              
            %Recreate filter
            expDopFinal = calculateFilter_wMobility(Cfg, modF(tagNum), modDuty, NFFTVel, vRangeExt, fc, c, fs, velocities(tagNum));

            %apply to fft across antennas
            RDASmall = squeeze(RDA(selR(tagNum),:,:));
            RDAngle = fftshift(fft(RDASmall.*squeeze(WinAnt2D(1,1,:))', NFFTAnt, 2)/ScaWinAnt,2);
            rangeBin(:,:)=squeeze(abs(RDAngle));
            normexpDopFinal=expDopFinal(1,:).*sigmoid;

            %Find maximum correlation
            normR2angle=squeeze(sum(rangeBin.*normexpDopFinal', 1));
            [aoa_val, aoa_ind]=max(normR2angle);

            aoa_f=vAngDeg(aoa_ind);

            angleAll(MeasIdx,tagNum)=aoa_f;
            angleIndAll(MeasIdx,tagNum)=aoa_ind;
                
             indices_New=nan;
             %[indices_New,Ahat1]=music_RAP_2D_Real(reshape(permute(Data.*Win2D.*mCalData,[1,3,2]),[],Cfg.FrmMeasSiz), 5, 0.02,1/fs,kf,fc,antenna_dist,NrChn,0,range_f-0.25,range_f+0.25,max(aoa_f-25,-50),min(aoa_f+25,50));%-50,50);%max(aoa_f-20,-50),min(aoa_f+20,50));      
             %[indices_New,Ahat1]=music_RAP_2D_Real(reshape(permute(Data.*Win2D.*mCalData,[1,3,2]),[],Cfg.FrmMeasSiz), 5, 0.1,1/fs,kf,fc,antenna_dist,NrChn,0,2,20,-50,50);%-50,50);%max(aoa_f-20,-50),min(aoa_f+20,50));      
             %[indices_New,Ahat1]=music_RAP_2D_Real(reshape(permute(Data,[1,3,2]),[],Cfg.FrmMeasSiz), 5, 0.1,1/fs,kf,fc,antenna_dist,NrChn,0,2,20,-50,50);%-50,50);%max(aoa_f-20,-50),min(aoa_f+20,50));      
             %[indices_New,Ahat1]=music_RAP_2D_Real(newData, 3, 0.05,1/fs,kf,fc,antenna_dist,NrChn,0);
              if isnan(indices_New)
                  indices_New=[aoa_f,range_f];
              end
             corrVal(MeasIdx, tagNum)=finalVal(tagNum,1);
            
             rangeAll(MeasIdx, tagNum)=range_f;
             angleAll(MeasIdx, tagNum)=aoa_f;
             velocityAll(MeasIdx, :)=velocities(:);
             [xtemp,ytemp]=pol2cart(deg2rad(aoa_f+90),range_f);
             
              if(drawHeatmap == 1)
                 tagHeatmap=heatmapPoint(tagHeatmap, [xtemp, ytemp], tagNum, 1);
             end
             
             corrVal(MeasIdx,tagNum)=finalVal(tagNum,1);
             
             rangeAllMUSIC(MeasIdx,tagNum)=indices_New(2);
             angleAllMUSIC(MeasIdx,tagNum)=indices_New(1);
             [xtemp,ytemp]=pol2cart(deg2rad(indices_New(1)+90),indices_New(2));
             xAllMUSIC(MeasIdx,tagNum)=-xtemp
             yAllMUSIC(MeasIdx,tagNum)=ytemp
             
             if(isempty(KalmanXk_prev) || isnan(KalmanXk_prev(1, tagNum)))
                KalmanXk_prev=NaN(4, length(modF));
                KalmanXk=NaN(4, length(modF));

                KalmanXk_prev(:, tagNum)=[rangeAllMUSIC(MeasIdx,tagNum); angleAllMUSIC(MeasIdx,tagNum); 0; 0];
                KalmanXk(:, tagNum)=[rangeAllMUSIC(MeasIdx,tagNum); angleAllMUSIC(MeasIdx,tagNum); 0; 0];
                KalmanCov(:,:,tagNum) = KalmanQ;
             else
                %TODO: Explain this
                dt=frameTime;
                newRange=rangeAllMUSIC(MeasIdx,tagNum);
                newAngle=angleAllMUSIC(MeasIdx,tagNum);
                %velOpt3 = velocityAll(MeasIdx,tagNum);
                Z = [rangeAll(MeasIdx,tagNum); angleAll(MeasIdx,tagNum); 
                        velocityAll(MeasIdx,tagNum); 0];
                %State Transition Matrix
                Phi = [1 0 dt 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

                % Prediction
                predictVal = Phi * KalmanXk(:, tagNum);
                predictCov = Phi * KalmanCov(:,:,tagNum) * Phi' + KalmanQ;

                % Kalman iteration
                KalmanError = Z - predictVal;
                measCov = KalmanM*predictCov*KalmanM' + KalmanR;

                mahalPenalty = 5;
                mahalToR = 10;
                mahalError = Z - KalmanXk(:, tagNum);
                mahalDist = sqrt(mahalError' * inv(measCov) * mahalError);
                weightedMahalDist = mahalToR ./ (1+ exp(0.5*(-mahalDist + mahalPenalty)));

                %Remembering parameter
                alphParam=0.7;

                %Scale covariance by current distance and old covariance - want some level
                %of permanence to covariance
                KalmanR = alphParam*KalmanR + (1-alphParam)*(mahalError*mahalError' + KalmanCov(:, :, tagNum));
                %mahalDist*mahalDist'
                weightedAll(MeasIdx)=mahalDist(1);

                %Gain
                K  = predictCov*KalmanM'*inv(measCov);
                allGain(MeasIdx)=det(K);
                if(mahalDist(1)>2 && KalmanCov(1,1,tagNum) < 3)
                    solnt(MeasIdx) = 1;
                    KalmanXk(:, tagNum)  = predictVal;
                else
                    solnt(MeasIdx) = 0;
                    KalmanXk(:, tagNum)  = predictVal + K*KalmanError;
                end

                %Update Estimates
                KalmanCov(:, :, tagNum)  = predictCov - K * measCov * K';
                KalmanCov(KalmanCov > 50) = 50;
                KalmanCov(KalmanCov < 0)  = 0;
                
                allCoswitchv(MeasIdx) = KalmanCov(1, 1, tagNum);
                allCov2(MeasIdx) = KalmanR(1,1);
                allCov3(MeasIdx) = measCov(1);
                [xtemp,wd]=pol2cart(deg2rad(KalmanXk(2)+90),KalmanXk(1));
            end
            rangeAllKAL(MeasIdx,tagNum) = KalmanXk(1,tagNum);
            angleAllKAL(MeasIdx,tagNum) = KalmanXk(2,tagNum);
             
             [xtemp,ytemp]=pol2cart(deg2rad(angleAllKAL(MeasIdx,tagNum)+90),rangeAllKAL(MeasIdx,tagNum));
             ytemp=-ytemp;
             
             if(drawHeatmap == 1)
                 tagHeatmap=heatmapPoint(tagHeatmap, [xtemp, ytemp], tagNum, 3);
             end
             pause(0.5)
             if connectAR == 1
                %Write single tag value to AR Here
                tagNum;
                xStr=num2str(-xtemp);
                yStr=num2str(ytemp);
                %posStr format:
                %Tag(TagNum),(xCoordinate),(yCoordinate),(zCoordinate)
                posStr=strcat('Tag', num2str(tagNum), ',', xStr, ',', yStr, ',0');
                fwrite(clientAR, strcat(posStr));
             end
         end
     end
     
     if drawGraphs
        figure(1)

        %title([num2str(MeasIdx) ' BknSub ' num2str(bknSubtract)]);
        axd(1)=subplot(1,2,1);   
        imagesc(1:Cfg.FrmMeasSiz,vRangeExt,abs(squeeze(RPExt(:,:,1))));
        %title(['Range FFT Amplitude RX' num2str(i) '']);
        xlabel('Chirps');
        ylabel('Range (m)');
        set(gca,'YDir','normal');
        ylim([vRangeExt(1) vRangeExt(end)]);
        colorbar
        %caxis([minR,maxR]);

        axd(2)=subplot(1,2,2);
        imagesc(normSig(1:end,5:end-4,1));
        %caxis([minRD,maxRD]);
        set(gca,'YDir','normal');
        
%         figure(2);
%         surf(mX,mY, JNorm); 
%         shading flat;
%         view(0,90);
%         axis equal
%         xlabel('x (m)');
%         ylabel('y (m)');
%         colormap('jet');
     end
    
   frameTime=toc
end
    
if(drawHeatmap == 1)
    tagHeatmap=heatmapIterate(tagHeatmap);
    heatmapSave(tagHeatmap)
end
    
xAll
yAll

if(connectAR==1)
    fclose(clientAR);    
end
%Data=reshape(tempD,measRounds,N,[],NrChn);
    
save(['../Data/' datestr(now, 'yymmdd-HHMMSS') filename],'Data','dtime','Brd','Cfg','N','NrChn','fs','measRounds','CalDat');
clear Brd;
end


function expDopFinal=calculateFilter_wMobility(Cfg,modF,modDuty,NFFTVel,vRangeExt,fc,c,fs,vel)
t=[0:Cfg.Perd:20*Cfg.Perd*10];%oversampling to get rid of the noise
tagFFT=round((Cfg.Perd*NFFTVel)*modF);%+128;
expDopFinal=zeros(length(vRangeExt),NFFTVel);
expdop=[];
sumR2d=[]; 
doppShift=fft(cos(2*pi*(2*vel*fc/c).*t), NFFTVel,2);


undopp_square=fft(square(2*pi*(modF.').*t(1:end-1),modDuty),NFFTVel,2);
%Time spent here?
if(vel>0)
    doppShiftPos=[doppShift(1:NFFTVel/2)   zeros(1,NFFTVel/2)];
    sq_wav=conv(doppShiftPos, undopp_square);
    expdopt=abs(sq_wav(1:NFFTVel));
else
    doppShiftNeg=flip([doppShift(1:NFFTVel/2)   zeros(1,NFFTVel/2)]);
    sq_wav=conv(doppShiftNeg, undopp_square);
    expdopt=abs(sq_wav(NFFTVel:end));
end
%sq_wav=conv(fft(square(2*pi*(modF.').*t(1:end-1),modDuty),NFFTVel,2),  fft(cos(2*pi*(2*vel*fc/c).*t),NFFTVel,2));
%sq_wav=conv(fft(square(2*pi*(modF.').*t(1:end-1),modDuty),NFFTVel,2),  fft(cos(2*pi*(2*vel*fc/c).*t),NFFTVel,2));
%expdopt=abs(sq_wav(1:NFFTVel));
%   sq_wav=square(2*pi*(modF.‘).*t,modDuty).*  cos(2*pi*(2*vel*fc/c)*t*Cfg.TRampUp);
%cos(2 * pi * kf* t * (2*vRangeExt(rk)) * (1/(c*fs))));
%sq_wav=square(2*pi*(modF.’).*t,modDuty);
%expdop=fftshift(fft(sq_wav,NFFTVel,2)); %reference signal
%expdopt=abs(fft(sq_wav,NFFTVel,2)); %reference signal
%expdop(rk,:)=expdopt;
expdopt=expdopt./max(expdopt);%normalize(expdopt,2,'range');
%[pks, locs]=findpeaks(expdopt,'MinPeakDistance',15,'MinPeakHeight',0.1);
%expdopt(~(expdopt>0.1 & expdopt > 1.6 * [0 expdopt(1:end-1)]))=0;
expdopt(expdopt<0.1)=0;
%expdopt(expdopt>0.5)=0.5;
%expDopFinal(1, locs)=expdopt(locs);
%expDopFinal(1, locs(1))=expDopFinal(1, locs(1))/2;
%expDopFinal(1, locs(end))=expDopFinal(1, locs(end))/2;
%expDopFinal(1,:)=expdopt; %expDopFinal(1,ffInd)=expdopt(ffInd); %
%expDopFinal(1,1:5)=0; %Zero out noise at beginning and end
%expDopFinal(1, end-5:end)=0;
expDopFinal=repmat(expdopt, length(vRangeExt), 1);
%figure;plot(normalize(expDopFinal(1,:),‘range’))

end


function template = calTemp_mob(v,tau,T,lambda,Tc,N,vrangeExtLeng)
% close all;
% lambda=0.0125;
% N=256; % number of chirps
% Tc=150e-6; % chirp time
% T=T_k*Tc; % on+off time
% tau=tau_k*Tc; % on time
bins=10000;

% creating DTFT
lhs=zeros(bins,1);
rhs=zeros(bins,1);
vel_ft=zeros(bins,1);

w1=linspace(-pi,pi,bins);

doppler_w = 2*pi*2*v/lambda*Tc;
res = find(w1>=doppler_w);
vel_ft(res(1)) = 1;
res = find(w1>=-doppler_w);
vel_ft(res(1)) = 1;

del=zeros(bins,1);
i=-floor(T/Tc):1:ceil(T/Tc);
for ii=1:length(i)
    if (2*pi*i(ii)*Tc/T <= pi) & (2*pi*i(ii)*Tc/T >= -pi)
        res = find(w1>=2*pi*i(ii)*Tc/T);
        del(res(1)) = 1;
    end
end
lhs=(sin(tau/Tc*w1(1:end)/2)./sin(w1(1:end)/2)).';
lhs=lhs.*del;

rhs=(sin(N*w1(1:end)/2)./sin(w1(1:end)/2)).';
% cont_ft = abs(2*pi*Tc/T*ifft(fft(lhs).*fft(rhs))); % convolution
cont_ft = conv(lhs,rhs,'same');

cont_ft_new = abs(conv(cont_ft,vel_ft,'same'));

if mod(N,2)==0
    ind = floor(1:(bins-1)/(N):bins); % sampling DTFT
    ind = ind(1:end-1);
else
    ind = floor(1:(bins-1)/(N-1):bins); % sampling DTFT
%     ind = ind(1:end-1);
end

template=cont_ft_new(ind); % Final DFT

% figure;
% plot(cont_ft);
% figure;
% plot(vel_ft);
% figure;
% plot(fftshift(template));
template=repmat(fftshift(template)',vrangeExtLeng,1);
end

function expDopFinal=calculateFilter(Cfg,modF,modDuty,NFFTVel,vRangeExt,fc,c,fs,vel)
    t=[0:Cfg.Perd:Cfg.FrmMeasSiz*Cfg.Perd*4];%oversampling to get rid of the noise
    tagFFT=round((Cfg.Perd*NFFTVel)*modF);%+128;
%       vFreqVel= ((vel  *(2.*fc))/c)*(Cfg.Perd*NFFTVel)

    expDopFinal=zeros(1,NFFTVel);

    %create the template (square wave)
%     sq_wav=conv(fft(square(2*pi*(modF.').*t,modDuty),NFFTVel,2),    fft(cos(2*pi*(2*vel*fc/c)*Cfg.TRampUp),NFFTVel,2));
%     expdopt=abs(sq_wav(1:NFFTVel));
%     sq_wav=square(2*pi*(modF.').*t,modDuty).*   cos(2*pi*(2*vel*fc/c)*t*Cfg.TRampUp);
    %cos(2 * pi * kf* t * (2*vRangeExt(rk)) * (1/(c*fs))));
    sq_wav=square(2*pi*(modF.').*t,modDuty);
    
    % apply FFT to capture sinc template
    %expdop=fftshift(fft(sq_wav,NFFTVel,2)); %reference signal 
    expdopt=abs(fft(sq_wav,NFFTVel,2)); %reference signal 

    % only keep the harmonics of the sync function 
    ind=1:round(NFFTVel/(tagFFT*2));
    fInd=tagFFT*ind+1;
    fIndL=NFFTVel-tagFFT*ind+1;
    ffInd=unique(sort([fInd,fIndL]));
    ffInd=ffInd(ffInd<NFFTVel & ffInd>1 );
    %ffInd=ffInd(ffInd>NFFTVel/2 & ffInd<NFFTVel);
    expDopFinal(ffInd)=expdopt(ffInd);

    %normalize the template
    %figure;plot(normalize(expDopFinal(1,:),'range'))
    expDopFinal=normalize(expDopFinal,2,'range');
end

%another template generator function 
function expDopFinal=calculateFilter1(Cfg,modF,modDuty,NFFTVel,vRangeExt)
%     modF=modF-3;
    t=linspace(0,Cfg.FrmMeasSiz*Cfg.Perd,20000);
    sq_wav=square(2*pi*(modF).*t,modDuty);
    sq_wav_n=zeros(1,Cfg.FrmMeasSiz);
    for i=1:Cfg.FrmMeasSiz
        sq_wav_n(i) = sq_wav(round((i-1)/Cfg.FrmMeasSiz*20000)+1);
    end
    sq_wav_n=(sq_wav_n+1)/2;
    expdopt=(fft(sq_wav_n,NFFTVel));
    expdopt=expdopt/norm(expdopt);
    expDopFinal=repmat(expdopt,length(vRangeExt),1,4);

%     [~,b,~,d]=findpeaks(abs(expdopt));
%     expdopt_t=zeros(size(expdopt));
%     expdopt_t(b(d>0.05))=expdopt(b(d>0.05));
%     expdopt_t=expdopt_t/norm(expdopt_t);
%     expDopFinal=repmat(expdopt_t,length(vRangeExt),1,4);
end

function idx=findGT(frameTime,GT_time)
 [val,idx]=min(abs(frameTime-GT_time));
    if val>=seconds(1)
        idx=-1;
    end
end
