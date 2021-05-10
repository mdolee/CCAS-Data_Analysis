clc; clear;
format long g 
tic
% scriptFolder = pwd;
% dataFolder = 'CCAS_vibration_data_sample/';
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
calib=xlsread('CX1_1027 - CX1_1027 - A - 20190726_210802.csv');
%data=xlsread('Inner_Tower - CX1_2612 - AC - 20201030_055108.csv');%Noisy Data
% data=xlsread('Inner_Tower - CX1_2612 - AC - 20210416_102111.csv');%Recent Data
data=xlsread('Mars_Platform - CX1_2613 - AC - 20210416_102100.csv');%Platform Data
% data=xlsread('Inner_Tower - CX1_2612 - AC - 20201101_015108.csv');%Calm Data

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
figure()
plot(freqx,PSDx)
xlabel("Freq (Hz)")
ylabel('|P1(f)|^2')
title("Power Specturm Density for Position X")
figure()
plot(freqy,PSDy)
xlabel("Freq (Hz)")
ylabel('|P1(f)|^2')
title("Power Specturm Density for Position Y")
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
    order = 3; width = 51; %For sgolayfilter
    Calib_X0 = cumtrapz(calib_t, cumtrapz(calib_t, calib)); %Accel to Disp
    Calib_X = sgolayfilt(Calib_X0,order,width); %Filter out general trend
    calib_disp = Calib_X0-Calib_X; %Find actual Accel Data
    Disp_X0 = cumtrapz(time, cumtrapz(time, accel));
    Disp_X = sgolayfilt(Disp_X0,order,width);
    disp=2*(Disp_X0-Disp_X); %Same process for actual data
end
function [PSD, freq] = PSD2(Fs, disp, calib_disp)
    %FFT transform for Calib Data
    Nc = length(calib_disp);
    Yc = fft(calib_disp);
    P2c = abs(Yc/Nc).^2;
    P1c = P2c(1:Nc/2+1);
    P1c(2:end-1)=2*P1c(2:end-1);
    freqc = 0:Fs/length(calib_disp):Fs/2;
    %FFT transform for Actual Data
    N = length(disp);
    Y = fft(disp);
    P2 = abs(Y/N).^2;
    P1 = P2(1:N/2+1);
    P1(2:end-1)=2*P1(2:end-1);
    freq = 0:Fs/length(disp):Fs/2;
    %Subtract Calib data from Actual Data
    PSDc=interp1(freqc, P1c, freq)'; %Interpolate to make vector lengths =
    PSD=P1-PSDc;
    PSD(PSD<0)=0; %Set Power <0 equal to 0
end
%Fit composite function to calibrated data (later)
