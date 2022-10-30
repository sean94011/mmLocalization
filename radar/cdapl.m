% Based on AN24_06 -- Range-Doppler basics
% Evaluate the range-Doppler map for a single channel
function cdapl(process, filename)
arguments
    process (1,1) double {mustBeInteger} = 0;
    filename (1,1) string {mustBeNonzeroLengthText} = 'nope_no_filename_specified';
end

% Variables for MMtro
maxRange            =   10;
chirpDur            =   125; 
c0                  =   2.998e8;
switchPer           =   625*10^-6;
measRounds          =   10000;
modF                =   1 / (2 * switchPer);
modDuty             =   50;
tDelay              =   0.1;

% (1) Connect to DemoRad: Check if Brd exists: Problem with USB driver
% (3) Configure RX
% (4) Configure TX
% (5) Start Measurements
% (6) Configure calculation of range profile and range doppler map for
% channel 1

%--------------------------------------------------------------------------
% Include all necessary directories
%--------------------------------------------------------------------------
CurPath = pwd();
addpath([CurPath,'/DemoRadBasic/DemoRadUsb']);
addpath([CurPath,'/DemoRadBasic/Class']);

if process == 0
    %--------------------------------------------------------------------------
    % Setup Connection
    %--------------------------------------------------------------------------
    Brd         =   TinyRad();

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
    % TODO figure out what CalDat does

    %--------------------------------------------------------------------------
    % Configure Up-Chirp and timing for the measurements
    %--------------------------------------------------------------------------
    Cfg.fStrt       =   24.00e9;                    %   Start frequency   
    Cfg.fStop       =   24.25e9;                    %   Stop frequency
    Cfg.TRampUp     =   chirpDur*10^-6;             %   UpChirp duration
    Cfg.Perd        =   Cfg.TRampUp+100e-6;         %   Period between measurements
    Cfg.N           =   100;                        %   Number of samples taken at start of chirp 
    Cfg.Seq         =   [1];                        %   Antenna transmit sequence
    Cfg.CycSiz      =   2;                          %   Number of buffers in the acquisition framework 2
    Cfg.FrmSiz      =   256;                        %   Number of chirp sequences for one measurement cycle
    Cfg.FrmMeasSiz  =   256;                        %   Number of chirps sequences for collecting IF data

    Brd.RfMeas(Cfg);

    %--------------------------------------------------------------------------
    % Read actual configuration
    %--------------------------------------------------------------------------
    NrChn           =   Brd.Get('NrChn');
    N               =   Brd.Get('N');
    fs              =   Brd.Get('fs');

else
    RadData = load(CurPath + "/Data" + filename);
    Brd = TinyRad();
    
    Cfg = RadData.Cfg;
    N = double(RadData.N);
    fs = RadData.fs;
    NrChn = double(RadData.NrChn);
    CalDat = RadData.CalDat;
    measRounds = RadData.measRounds;
    Cfg.FrmMeasSiz = double(Cfg.FrmMeasSiz);
end

%--------------------------------------------------------------------------
% Configure Signal Processing
%--------------------------------------------------------------------------
% Processing of range profile
Win2D           =   Brd.hanning(N, Cfg.FrmMeasSiz); % Make this 3D for AoA
ScaWin          =   sum(Win2D(:,1));
NFFT            =   2^10;
NFFTVel         =   2^8;
kf              =   (Cfg.fStop - Cfg.fStrt)/Cfg.TRampUp;
vRange          =   [0:NFFT-1].'./NFFT.*fs.*c0/(2.*kf);
fc              =   (Cfg.fStop + Cfg.fStrt)/2;

RMin            =   0.2;
RMax            =   min(maxRange,((N/Cfg.TRampUp)*c0)/(4*(250e6/Cfg.TRampUp)));

[~, RMinIdx]    =   min(abs(vRange - RMin));
[~, RMaxIdx]    =   min(abs(vRange - RMax));
vRangeExt       =   vRange(RMinIdx:RMaxIdx);
[~, ROffset]    =   max(vRangeExt>=1.0);

WinVel          =   Brd.hanning(Cfg.FrmMeasSiz);
ScaWinVel       =   sum(WinVel);
WinVel2D        =   repmat(WinVel.',numel(vRangeExt),1);

vFreqVel        =   [-NFFTVel./2:NFFTVel./2-1].'./NFFTVel.*(1/Cfg.Perd);
vVel            =   vFreqVel*c0/(2.*fc); 

mCalData        =   permute(repmat(CalDat,1,Cfg.FrmMeasSiz,N-1),[3 2 1]);

% velArr          =   (0);
velArr          =   (-3 : 0.1 : 3);

%--------------------------------------------------------------------------
% Measure and calculate range-Doppler Map
%--------------------------------------------------------------------------

visualizeEvery  =   1;

for i           =   1:floor(measRounds / visualizeEvery)
    MeasIdx     =   i * visualizeEvery;
    tic;
    
    if process == 0
        % Get the measures of the current frame and only take 1 channel.
        % DataFrame in dimensions (chirps * samples) * (channels);
        % resized into (samples) * (chirps).
        DataFrame   =   Brd.BrdGetData();
        Data        =   reshape(DataFrame(:,1), N, []); 
    else
        Data        =   reshape(squeeze(RadData.Data(:,MeasIdx,1)), N, []);
    end
    
    % Background subtraction by first chirp of each frame
    Data        =   Data - Data(:,1);

    % Calculate range profile including calibration 
    % Perform 1st FFT (2D range FFT) w/ NFFT bins, then take range
    % bins [RMinIdx, RMaxIdx) which should contain the tag
    % This is done wrt the samples (dim=1). New data becomes (bins)
    % * (chirps in each frame).
    RP          =   2*fft(Data.*Win2D,NFFT,1).*Brd.FuSca/ScaWin;
    RPExt       =   RP(RMinIdx:RMaxIdx,:);
    % TODO Add calibration data

    % Perform Doppler-FFT to get velocity. 
    % This is done wrt the chirps to compare different phases
    % across each chirp (dim=2). New data becomes (rangeFFT bins) *
    % (dopplerFFT bins).
    RD          =   fft(RPExt.*WinVel2D, NFFTVel, 2)./ScaWinVel;
    
    % Match sinc function template. The template is of size (doppFFT bins),
    % thus repmat-ed (rangeFFT bins) times to match the size of RDNorm. 
    RDDiff      =   abs(RD) - min(abs(RD));
    RDNorm      =   RDDiff ./ max(RDDiff);
    
    amplitude   =   0;
    corrAll     =   zeros(length(velArr), length(vRangeExt));
    templateAll =   zeros(length(velArr), NFFTVel);
    for vIdx = 1:length(velArr)
        vel         =   velArr(vIdx);
        template    =   sincTempVel0(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc);
%         template    =   sincTemplate(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc, vel);
        templateAll(vIdx, :) = template;
        corr        =   abs(RDNorm).*repmat(abs(template), size(RDNorm, 1), 1);
        corrSum     =   sum(abs(corr), 2);
        [ampArr, ampIdxArr] = max(corrSum(ROffset:end));
        corrAll(vIdx, :) = corrSum;

%         if length(ampArr) > 1 && ampArr(1) > amplitude\
        if ampArr(1) > amplitude
            amplitude = ampArr(1);
            rangeIdx = ampIdxArr(1) + ROffset - 1;
            velIdx = vIdx;
        end
    end
    
    frameTime = toc;
    figure(1)
    axd(1) = subplot(1, 2, 1);
    imagesc(1:Cfg.FrmMeasSiz, vRangeExt, abs(RPExt));
    
%     imagesc(1:Cfg.FrmMeasSiz, vRange, abs(RP));
    xlabel('Chirps');
    ylabel('Range (m)');
    set(gca, 'YDir', 'normal');
    ylim([vRangeExt(1) vRangeExt(end)]);
%     title(['Range: ', num2str(vRangeExt(rangeIdx)), 'm; Vel: ', num2str(velArr(velIdx)), 'm/s']);
    colorbar

    axd(2) = subplot(1, 2, 2);
    plotRD = abs(RD(:, 15:end-15));
    plotRD(rangeIdx, :) = max(max(plotRD));
    imagesc(vVel(15:end-15), vRangeExt, plotRD);
    xlabel('Velocity Bins');
    ylabel('Range (m)');
    set(gca, 'YDir', 'normal');
    if process == 0
        title(['Current frametime: ', num2str(frameTime), ' seconds']);
    else
        title(['Data captured time: ', RadData.dtime(MeasIdx, :)]);
    end
%     colorbar
    
%     axd(3) = subplot(2, 2, 3);
%     imagesc(vRangeExt, velArr, corrAll);

%     axd(4) = subplot(2, 2, 4);
%     findpeaks(corrAll(velIdx, ROffset:end))

%     axd(4) = subplot(2, 2, 4);
%     imagesc(vVel(15:end-15), velArr, abs(templateAll(:, 15:end-15)));
    
    % Somehow matlab needs at least 0.05 seconds before next graph
    % rendering to finish proper rescaling based on graph's window
    % size. Less than this will make the graph fixed size.
    pause(tDelay);
end

clear Brd

end

function expdopt = sincTempVel0(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc)
    t = linspace(0,Cfg.FrmMeasSiz*Cfg.Perd,20000);
    sq_wav = square(2*pi*(modF).*t,modDuty);
    sq_wav_n = zeros(1,Cfg.FrmMeasSiz);
    for i=1:Cfg.FrmMeasSiz
        sq_wav_n(i) = sq_wav(round((i-1)/Cfg.FrmMeasSiz*20000)+1);
    end
    sq_wav_n=(sq_wav_n+1)/2;
    expdopt=fft(sq_wav_n,NFFTVel);
    expdopt=expdopt/norm(expdopt);
    % expDopFinal=repmat(expdopt,length(vRangeExt),1,4);
end

function expdopt = sincTemplate(Cfg, modF, modDuty, NFFTVel, vRangeExt, fc, vel)
    c0 = 2.998e8;

    t = (0 : Cfg.Perd : 20*Cfg.Perd*10);
    tagFFT = round((Cfg.Perd*NFFTVel)*modF);
    doppShift=fft(cos(2*pi*(2*vel*fc/c0).*t), NFFTVel, 2); % Doppler shift of sq. wave

    undopp_square=fft(square(2*pi*(modF.').*t(1:end-1),modDuty), NFFTVel, 2);
    if (vel>0)
        doppShiftPos = [doppShift(1:NFFTVel/2)   zeros(1,NFFTVel/2)];
        sq_wav = conv(doppShiftPos, undopp_square);
        expdopt = abs(sq_wav(1:NFFTVel));
    elseif (vel < 0)
        doppShiftNeg = flip([doppShift(1:NFFTVel/2)   zeros(1,NFFTVel/2)]);
        sq_wav = conv(doppShiftNeg, undopp_square);
        expdopt = abs(sq_wav(NFFTVel:end));
    else
        expdopt = abs(undopp_square);
    end

    expdopt = expdopt ./ max(expdopt);
    expdopt(expdopt<0.1) = 0;
end