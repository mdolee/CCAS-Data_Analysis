clc; clear;
format long g 
tic
scriptFolder = pwd;
calDataFolder = 'CCAS_vibration_data_sample/';
dataFolder = "C:\Users\hl6726\Box\CCAS-Upgrade\CCAS_Folder\Analysis\aoa_video\2021-05-03\";
cx1File = "Inner_Tower - CX1_2612 - AC - 20210503_075432.csv";

% dataFile = fullfile(scriptFolder,dataFolder,'Inner_Tower - CX1_2612 - AC - 20201101_015108.csv');
% dataFile = fullfile(scriptFolder,dataFolder,'Outer_Tower - CX1_2613 - AC - 20201101_015058.csv');

% HL update 11/11/20
% [t,acc] = readCXfile(dataFile); % read CX file
% time = t - t(1);  % normalize to the start time
% g = 9.81; % gravitational constant m/s^2
% Accel_X = acc(:,1) * g; % scale to m/s^2
% Accel_Y = acc(:,2) * g; % scale to m/s^2
% Accel_Z = acc(:,3) * g; % scale to m/s^2

% Read In Calibration and Raw Data
%data = xlsread(dataFile);
calib=xlsread(fullfile(scriptFolder,calDataFolder,'CX1_1027 - CX1_1027 - A - 20190726_210802.csv'));
%data=xlsread('Inner_Tower - CX1_2612 - AC - 20201030_055108.csv');%Noisy Data
% data=xlsread('Inner_Tower - CX1_2612 - AC - 20210416_102111.csv');%Recent Data
% data=xlsread('Mars_Platform - CX1_2613 - AC - 20210416_102100.csv');%Platform Data
% data=xlsread('Inner_Tower - CX1_2612 - AC - 20201101_015108.csv');%Calm Data
data=xlsread(fullfile(dataFolder,cx1File));

%DA Update 12/21/20
%Raw Data to Clean Data, Aceel to Disp
[C_Ax,C_Ay,C_Az,C_t] = cleandata(calib);
[Accel_X,Accel_Y,Accel_Z,time] = cleandata(data);

[dispX,calib_dispX] = accel2disp(Accel_X,time,C_Ax,C_t);
[dispY,calib_dispY] = accel2disp(Accel_Y,time,C_Ay,C_t);
[dispZ,calib_dispZ] = accel2disp(Accel_Z,time,C_Az,C_t);
% Savitzky-Golay filtering to remove large scale trend which is likely
% intrinsic to the accelerometer rather than anything physical.
% order 3 and filter window size 51 samples. (Signal processing toolbox needed).

% DA update 12/8/20
% Conversion from Time Domain to Freq Domain
Fs=1/mean(diff(time)); %Sampling Frequency, Nyquist=Sampling/2
[PSDx, freqx]=PSD2(Fs, dispX/125e-6, calib_dispX/125e-6); %disp in arcsec
[PSDy, freqy]=PSD2(Fs, dispY/125e-6, calib_dispY/125e-6); %disp in arcsec
[PSDz, freqz]=PSD2(Fs, dispZ/125e-6, calib_dispZ/125e-6); %disp in arcsec
%Power of origional Displacement
fprintf("Standard Dev of Origional Disp x: %1.8f microns\n",sqrt(std(dispX).^2-std(calib_dispX).^2)*1e6)
fprintf("Standard Dev of Origional Disp y: %1.8f microns\n",sqrt(std(dispY).^2-std(calib_dispY).^2)*1e6)
fprintf("Standard Dev of Origional Disp z: %1.8f microns\n",sqrt(std(dispZ).^2-std(calib_dispZ).^2)*1e6)
%DA update 12/21/20
%Save PSD in each dimension to a matfile
save('PSDPlatform.mat','Fs','PSDx','PSDy','PSDz','freqx','freqy','freqz');
%%% Plotting Section %%%
% figure()
% hold on; grid on;
% plot(time, Accel_X)
% xlabel("Time (sec)")
% ylabel("Acceleration (m/s^2)")
% title("X Accleration vs Time")
% figure()
% hold on; grid on;
% plot(time, Accel_Y)
% xlabel("Time (sec)")
% ylabel("Acceleration (m/s^2)")
% title("Y Accleration vs Time")
% figure()
% hold on; grid on;
% plot(time, Accel_Z)
% xlabel("Time (sec)")
% ylabel("Acceleration (m/s^2)")
% title("Z Accleration vs Time")
% 
% figure()
% hold on; grid on;
% plot(time, dispX)
% xlabel("Time (sec)")
% ylabel("Position (m)")
% title("X Displacement vs Time")
% figure()
% hold on; grid on;
% plot(time, dispY)
% xlabel("Time (sec)")
% ylabel("Position (m)")
% title("Y Displacement vs Time")
% figure()
% hold on; grid on;
% plot(time, dispZ)
% xlabel("Time (sec)")
% ylabel("Position (m)")
% title("Z Displacement vs Time")
% 
fig = figure();
fig.Position = [100,100,900,900];
fig.Color = [1 1 1];

ax = subplot(1,1,1);
hold on;
ph = plot(ax,freqx,PSDx); 
ph(1).LineStyle = 'none'; ph(1).Marker='o'; ph(1).MarkerFaceColor='w';ph(1).LineWidth=1;ph(1).MarkerSize=3;
ph = plot(ax,freqy,PSDy);
ph(1).LineStyle = 'none'; ph(1).Marker='x'; ph(1).MarkerFaceColor='w';ph(1).LineWidth=1;ph(1).MarkerSize=3;
lh = yline(1); lh.LineWidth=2;
hold off;
ax.XScale = 'log';
ax.YScale = 'log';
ax.XLim=[0.1,30];
ax.YLim=[1e-19,1e-6]/0.125e-3^2;
ax.XTick=[0.1,0.2,0.3,0.4,0.5,1,2,3,4,5,10,20,30];
ax.XTickLabel = cellstr(string(round(ax.XTick,1))');
ax.GridLineStyle = ':';
legend(["x","y","1'' fwhm"]);
% th = title(ax,"(A) Pier data (as calibration)",'Units','Normalized', ...
%      'Position',[0.27 0.95 0],'FontName','Calibri','FontSize',16);
xlabel(ax,"Frequency [Hz]",'FontName','Calibri','FontSize',12);
ylabel(ax,"PSD [''^2/Hz]",'FontName','Calibri','FontSize',12);
grid on; box on;
pos = ax.Position;
pos(1)=0.075; pos(3) = 0.9; ax.Position = pos;
ax2 = axes('Position',pos,'XAxisLocation','top','YAxisLocation','left','Color','none');
ax2.XScale = ax.XScale;
ax2.YScale = ax.YScale;
ax2.XLim=ax.XLim;
ax2.YLim=ax.YLim;
ax2.YTick=ax.YTick;
ax2.XTick=ax.XTick;
ax2.XTickLabel = cellstr(string(round(1./ax2.XTick/2.*1000,0))');
ax2.XAxis.FontSize = 8;
xlabel(ax2,'Sampling Time [ms]','FontName','Calibri','FontSize',12);
% loglog(freqx,PSDx)
% xlim([0.1,25])
% ylim([1e-19,1e-6]/0.125e-3^2)
% xlabel("Freq (Hz)")
% ylabel('PSD [''^2/Hz]')
% title("Power Specturm Density for Position X")

% figure(2)
% loglog(freqy,PSDy)
% xlim([0.1,25])
% ylim([1e-19,1e-6]/0.125e-3^2)
% xlabel("Freq (Hz)")
% ylabel('PSD [''^2/Hz]')
% title("Power Specturm Density for Position Y")

% figure()
% hold on;
% plot(freqz,PSDz)
% xlabel("Freq (Hz)")
% ylabel('|P1(f)|')
% title("Power Specturm Density for Position Z")
toc
%DA Update 12/21/20
%Functions used in Script
function [AccX,AccY,AccZ,time] = cleandata(data)
    time=(data(:,3)-data(1,3))*86400; %Seconds, start time at 0
    x=isnan(data(:,4)); %NAN entries for X direction (they are the same)
    time(x)=[]; %Reomving time locations with NAN entries
    temp=data(:,4);AccX=temp(~isnan(temp)); %Removing NAN entries
    temp=data(:,5);AccY=temp(~isnan(temp));
    temp=data(:,6);AccZ=temp(~isnan(temp));
    AccX=(AccX-mean(AccX))*9.81; %m/s^2, convert from G's to m/s^2
    AccY=(AccY-mean(AccY))*9.81; %m/s^2, convert from G's to m/s^2
    AccZ=(AccZ-mean(AccZ))*9.81; %m/s^2, convert from G's to m/s^2
end
function [disp,calib_disp] = accel2disp(accel,time,calib,calib_t)
    order = 3; width = 251; %For sgolayfilter
    Calib_X0 = cumtrapz(calib_t, cumtrapz(calib_t, calib)); %Accel to Disp
    Calib_X = sgolayfilt(Calib_X0,order,width); %Filter out general trend
    calib_disp = Calib_X0-Calib_X; %Find actual Accel Data
    Disp_X0 = cumtrapz(time, cumtrapz(time, accel));
    Disp_X = sgolayfilt(Disp_X0,order,width);
    disp=(Disp_X0-Disp_X); %Same process for actual data
end
function [PSD, freq] = PSD2(Fs, disp, calib_disp)
    %FFT transform for Calib Data
    Nc = length(calib_disp);
    Yc = fft(calib_disp);
    P2c = abs(Yc/Nc).^2;
    P1c = P2c(1:Nc/2+1);
    P1c = 2 * P1c;  % Only half of the PSD is taken, 
                    % hence x2 to compensate for the other side.
    P1c = 4 * P1c;  % This is due to the fact that displacement of spot 
                    % doubled on reflection, hence x4 in PSD space.
%     P1c(2:end-1)=2*P1c(2:end-1);
%     freqc = 0:Fs/length(calib_disp):Fs/2;
    freqc = linspace(0,Fs/2,Nc/2+1); %0:Fs/length(calib_disp):Fs/2;
    %FFT transform for Actual Data
    N = length(disp);
    Y = fft(disp);
    P2 = abs(Y/N).^2;
    P1 = P2(1:N/2+1);
    P1 = 2 * P1;  % Only half of the PSD is taken, 
                  % hence x2 to compensate for the other side.
    P1 = 4 * P1;  % This is due to the fact that displacement of spot 
                  % doubled on reflection, hence x4 in PSD space.
%     P1(2:end-1)=2*P1(2:end-1);
%     freq = 0:Fs/length(disp):Fs/2;
    freq = linspace(0,Fs/2,N/2+1);
    %Subtract Calib data from Actual Data
    PSDc=interp1(freqc, P1c, freq)'; %Interpolate to make vector lengths =
    PSD=P1-PSDc;
    PSD(PSD<0)=0; %Set Power <0 equal to 0
end
%Fit composite function to calibrated data (later)
