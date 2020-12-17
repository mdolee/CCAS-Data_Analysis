clc; clear; close all;
format long g 
tic
%subtract baseline from data PSD
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
data=xlsread('Inner_Tower - CX1_2612 - AC - 20201101_015108.csv');
%data=xlsread('Outer_Tower - CX1_2613 - AC - 20201101_022108.csv');
Calib_Ax=calib(:,4)*9.81; %m/s^2
Calib_Ay=calib(:,5)*9.81; %m/s^2
Calib_Az=calib(:,6)*9.81; %m/s^2
Calib_t=(calib(:,3)-calib(1,3))*86400; %seconds

Accel_X=data(:,4)*9.81; %m/s^2
Accel_Y=data(:,5)*9.81; %m/s^2
Accel_Z=data(:,6)*9.81; %m/s^2

time=(data(:,3)-data(1,3))*86400; %Seconds
x=find(isnan(Accel_X)); %NAN entries for X,Y,Z direction (they are the same)
time(x)=[]; %Reomving time locations with NAN entries

%Removing NAN entries
Accel_X=Accel_X(~isnan(Accel_X));
Accel_Y=Accel_Y(~isnan(Accel_Y));
Accel_Z=Accel_Z(~isnan(Accel_Z));

%Normalizing data around 0
Calib_Ax=Calib_Ax-mean(Calib_Ax);
Calib_Ay=Calib_Ay-mean(Calib_Ay);
Calib_Az=Calib_Az-mean(Calib_Az);
Accel_X=Accel_X-mean(Accel_X);
Accel_Y=Accel_Y-mean(Accel_Y);
Accel_Z=Accel_Z-mean(Accel_Z);
%DA Update 11/24/20
%Subtracting from average results in data starting/ending at zero

% HL update 11/11/20
% Integral of acceleration to velocity
Calib_Vx = cumtrapz(Calib_t, Calib_Ax);
Calib_Vy = cumtrapz(Calib_t, Calib_Ay);
Calib_Vz = cumtrapz(Calib_t, Calib_Az);
Vel_x = cumtrapz(time,Accel_X);
Vel_y = cumtrapz(time,Accel_Y);
Vel_z = cumtrapz(time,Accel_Z);

% HL update 11/11/20
% Integral of velocity to displacement
Calib_Xx0 = cumtrapz(Calib_t, Calib_Vx);
Calib_Xy0 = cumtrapz(Calib_t, Calib_Vy);
Calib_Xz0 = cumtrapz(Calib_t, Calib_Vz);
Disp_x0 = cumtrapz(time,Vel_x);
Disp_y0 = cumtrapz(time,Vel_y);
Disp_z0 = cumtrapz(time,Vel_z);

% DA update 11/24/20
% Disp_yft=fft(Disp_y);
% Disp_zft=fft(Disp_z);
% HL update 11/11/20
% Savitzky-Golay filtering to remove large scale trend which is likely
% intrinsic to the accelerometer rather than anything physical.
% order 3 and filter window size 51 samples. (Signal processing toolbox
% needed).
order = 3; width = 51;
Calib_Xx = sgolayfilt(Calib_Xx0,order,width);
Calib_Xy = sgolayfilt(Calib_Xy0,order,width);
Calib_Xz = sgolayfilt(Calib_Xz0,order,width);
Disp_x = sgolayfilt(Disp_x0,order,width);
Disp_y = sgolayfilt(Disp_y0,order,width);
Disp_z = sgolayfilt(Disp_z0,order,width);

% HL update 11/11/20
% Calibrate out the large scale trend
Calib_Xx = Calib_Xx0 - Calib_Xx;
Calib_Xy = Calib_Xy0 - Calib_Xy;
Calib_Xz = Calib_Xz0 - Calib_Xz;
Disp_x = Disp_x0 - Disp_x;
Disp_y = Disp_y0 - Disp_y;
Disp_z = Disp_z0 - Disp_z;

% DA update 12/8/20
% Conversion from Time Domain to Freq Domain
Fs=1/mean(diff(time)); %Sampling Frequency
[PSDCx, freqcx]=PSD(Fs, Calib_Xx);
[PSDCy, freqcy]=PSD(Fs, Calib_Xy);
[PSDCz, freqcz]=PSD(Fs, Calib_Xz);

[psdx, freqx]=PSD(Fs, Disp_x);
[psdy, freqy]=PSD(Fs, Disp_y);
[psdz, freqz]=PSD(Fs, Disp_z);

%Subtracting out PSD from Calibrated Accelerometer
psdcx=interp1(freqcx, PSDCx, freqx)';
psdcy=interp1(freqcy, PSDCy, freqy)';
psdcz=interp1(freqcz, PSDCz, freqz)';
PSDx=psdx - psdcx;
PSDy=psdy - psdcy;
PSDz=psdz - psdcz;

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

% % HL update 11/11/20
% figure()
% hold on; grid on;
% plot(time, Disp_x)
% xlabel("Time (sec)")
% ylabel("Position (m)")
% title("X Displacement vs Time")

% HL update 11/11/20
% figure()
% hold on; grid on;
% plot(time, Disp_y)
% xlabel("Time (sec)")
% ylabel("Position (m)")
% title("Y Displacement vs Time")

% HL update 11/11/20
% figure()
% hold on; grid on;
% plot(time, Disp_z)
% xlabel("Time (sec)")
% ylabel("Position (m)")
% title("Z Displacement vs Time")

% DA update 12/8/20
figure()
plot(freqx,PSDx)
xlabel("Freq (Hz)")
ylabel('|P1(f)|')
title("Power Specturm Density for Position X")

% DA update 12/8/20
figure()
plot(freqy,PSDy)
xlabel("Freq (Hz)")
ylabel('|P1(f)|')
title("Power Specturm Density for Position Y")

% DA update 12/8/20
figure()
plot(freqz,PSDz)
xlabel("Freq (Hz)")
ylabel('|P1(f)|')
title("Power Specturm Density for Position Z")

%Finding Integrated Power along Frequency Domain Range
PSDx=PSDx(PSDx>=0);
PSDy=PSDy(PSDy>=0);
PSDz=PSDz(PSDz>=0);

Powerx=sqrt(trapz(PSDx))*1e6;
tempx=sqrt(std(Disp_x).^2-std(Calib_Xx).^2)*1e6;
Powery=sqrt(trapz(PSDy))*1e6;
tempy=sqrt(std(Disp_y).^2-std(Calib_Xy).^2)*1e6;
Powerz=sqrt(trapz(PSDz))*1e6;
tempz=sqrt(std(Disp_z).^2-std(Calib_Xz).^2)*1e6;

fprintf("Integrated Power of Disp x: %1.8f\n",Powerx)
fprintf("Standard Dev of Disp x: %1.8f\n",tempx)
fprintf("Integrated Power of Disp y: %1.8f\n",Powery)
fprintf("Standard Dev of Disp y: %1.8f\n",tempy)
fprintf("Integrated Power of Disp z: %1.8f\n",Powerz)
fprintf("Standard Dev of Disp x: %1.8f\n\n",tempz)
toc
function [P1, freq] = PSD(Fs, data)
    N = length(data);
    Y = fft(data);
    P2 = abs(Y/N).^2;
    P1 = P2(1:N/2+1);
    P1(2:end-1)=2*P1(2:end-1);
    freq = 0:Fs/length(data):Fs/2;
end
%Random number between 0 and 2pi as phase, rand(0,2*pi,N)
%Sgolayfit through power spectrum to get A
%Ax = sqrt(PSDx)
%Fourier Signal = A*exp(i*phase)
%Inverse Fourier Tranform back to Time Domain
%Mock up displacement signal for 5 min
%Variance check agasint power you inputted
%Min, Max, Nominal vibration case of data -HL
