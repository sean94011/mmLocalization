% AN24_02 -- FMCW Basics 
function collectDataAuto(filename, chirpDur)
filename='10m-together-moving-together-far';
chirpDur=150;
% (1) Connect to TinyRad: Check if Brd exists: Problem with USB driver
% (3) Configure RX
% (4) Configure TX
% (5) Start Measurements
% (6) Configure calculation of range profile
% Configure script
%--------------------------------------------------------------------------
% Include all necessary directories

%--------------------------------------------------------------------------
clientAR = tcpip("127.0.0.1", 5311, 'NetworkRole', 'client');  %192.168.1.2
fopen(clientAR);

pause(1)

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
CalDat          =   CalDat(1:4); %Note number of antennas, should be 1-8 or 1-4

%--------------------------------------------------------------------------
% Configure Up-Chirp and timing for the measurements
%--------------------------------------------------------------------------
Cfg.fStrt       =   24.0e9;                    %   Start frequency   
Cfg.fStop       =   24.25e9;                   %   Stop frequency
Cfg.TRampUp     =   chirpDur*10^-6;            %   UpChirp duration 
Cfg.Perd        =   Cfg.TRampUp+25e-6;         %   Period between measurements
Cfg.N           =   chirpDur/2;                %   Number of samples taken at start of chirp 
Cfg.Seq         =   [1];                       %   Antenna transmit sequence
Cfg.CycSiz      =   2;                         %   Number of buffers in the acquisition framework 2 are required
Cfg.FrmSiz      =   256;                       %   Number of chirp sequences for one measurement cycle
Cfg.FrmMeasSiz  =   256;                       %   Number of chirps sequences for collecting IF data

Brd.RfMeas(Cfg);


%--------------------------------------------------------------------------
% Read actual configuration
%--------------------------------------------------------------------------
N               =   Brd.Get('N');
NrChn           =   Brd.Get('NrChn');
fs              =   Brd.Get('fs');

%--------------------------------------------------------------------------
% Check TCP/IP data rate:
% 16 Bit * Number of Enabled Channels * Number of Samples are measured in
% the interval TInt. If the data rate is too high, than frames can be losed
%--------------------------------------------------------------------------
DataRate        =   16*NrChn*N.*Cfg.FrmMeasSiz./(Cfg.Perd*Cfg.FrmSiz);
disp(['DataRate: ', num2str(DataRate/1e6), ' MBit/s'])

%--------------------------------------------------------------------------
% Configure Signal Processing
%--------------------------------------------------------------------------
disp('Get Measurement data')
measRounds=4000;%10/(Cfg.FrmMeasSiz*Cfg.Perd);

freqList = [24266936299.2922, 24217961654.894, 24169184290.0302, 24072216649.9499, 24000000000];

%fcs = repelem(freqList, ceil(measRoundslength(freqList)/length(freqList)));
bw = 0.4e9;
%startFreqs = fcs - bw/2;
%stopFreqs = fcs + bw/2;
tic

pause(2)
posn = [];
time1 = datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS');
currentPos = [0 0];

for i=1:100
    Brd.BrdRst();
    while(length(currentPos) <= i || 0>currentPos(i+1) || 25<currentPos(i+1))
        posn(i) = round((rand-0.5)*2);
        currentPos(i+1) = posn(i)+ currentPos(i);
    end
    
    fwrite(clientAR, num2str(currentPos(i+1)));
    disp(currentPos(i+1));
    tic
    for Idx= 1:measRounds
    %     Cfg.fStrt       =   startFreqs(Idx);                    %   Start frequency   
    %     Cfg.fStop       =   stopFreqs(Idx);                    %   Stop frequency
    %     Brd.BrdRst();
    %     Brd.RfMeas(Cfg);
       if(Idx == 50)
        disp(1);
       end
       disp(Idx)
       Data(:,Idx,:)=Brd.BrdGetData();
       time2=datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS');
       dtime(Idx, :) = time2;
       pause(max(0, 0.1-seconds(time2-time1)));
       time1=datetime('now','Format','yyyy-MM-dd HH:mm:ss.SSS');
    end
    disp(i);
    save([[CurPath, '/Data/'] datestr(now, 'yymmdd-HHMMSS') filename i],'Data','dtime','Brd','Cfg','N','NrChn','fs','measRounds','CalDat');%, 'fcs', 'bw');
    Data=[]; dtime=[];
    toc
    pause(10);
end
%Data=reshape(tempD,measRounds,N,[],NrChn);

save([[CurPath, '/Data/'] datestr(now, 'yymmdd-HHMMSS') filename 'gt'],'posn');%, 'fcs', 'bw');
clear Brd;
end
