clc; clear; close all;
tic
%%% Import section
m = matfile('Mirrorframe.mat');FsM=m.FPS;
PSDxV=m.PSDx;PSDyV=m.PSDy;freqV=m.freq;timeV=m.time;
centroidxV=m.centroidx;centroidyV=m.centroidy;

m = matfile('Calibration.mat');FsC=m.FPS;
PSDxC=m.PSDx;PSDyC=m.PSDy;freqC=m.freq;timeC=m.time;
centroidxC=m.centroidx;centroidyC=m.centroidy;

m = matfile('PSDInnerAccel.mat');Fs=m.Fs;
PSDx=m.PSDx;PSDy=m.PSDy;freq=m.freqx;

m = matfile('PSDPlatform.mat');FsP=m.Fs;
PSDxP=m.PSDx;PSDyP=m.PSDy;freqP=m.freqx;

%%% Calucaltion Section
len=size(centroidxV);
for i=1:len(1)
    [tempPSDx, ~] = time2PSD(FsM, centroidxV(i,:)); %motion in Arcsec
    [tempPSDy, freqnew] = time2PSD(FsM, centroidyV(i,:));
    [tempPSDxC, ~] = time2PSD(FsM, centroidxC(i,:)); %motion in Arcsec
    [tempPSDyC, freqnewC] = time2PSD(FsM, centroidyC(i,:));
    PSDxNewV(i,:)=tempPSDx;
    PSDyNewV(i,:)=tempPSDy;
    PSDxNewC(i,:)=tempPSDxC;
    PSDyNewC(i,:)=tempPSDyC;
end
PSDxNewV=mean(PSDxNewV);PSDyNewV=mean(PSDyNewV);
PSDxNewC=mean(PSDxNewC);PSDyNewC=mean(PSDyNewC);
%%% Output Section
[maxPSDx,index]=max(PSDx);[maxPSDy,indexY]=max(PSDy);
[maxPSDxV,indexXV]=max(PSDxV);[maxPSDyV,indexYV]=max(PSDyV);
[maxPSDxNV,indexXNV]=max(PSDxNewV);[maxPSDyNV,indexYNV]=max(PSDyNewV);
[maxPSDxP,indexXP]=max(PSDxP);[maxPSDyP,indexYP]=max(PSDyP);
[maxPSDxC,indexXC]=max(PSDxNewC);[maxPSDyC,indexYC]=max(PSDyNewC);
fprintf("Calibration Values\n")
fprintf("The Max Power in X occurs at %1.4f Hz w/ a value of %1.4e\n",freqC(indexXC),maxPSDxC)
fprintf("The Max Power in Y occurs at %1.4f Hz w/ a value of %1.4e\n\n",freqC(indexYC),maxPSDyC)
fprintf("Video Values\n")
fprintf("The Max Power in X occurs at %1.4f Hz w/ a value of %1.4e\n",freqV(indexXV),maxPSDxV)
fprintf("The Max Power in Y occurs at %1.4f Hz w/ a value of %1.4e\n\n",freqV(indexYV),maxPSDyV)
fprintf("New Video Values\n")
fprintf("The Max Power in X occurs at %1.4f Hz w/ a value of %1.4e\n",freqnew(indexXNV),maxPSDxNV)
fprintf("The Max Power in Y occurs at %1.4f Hz w/ a value of %1.4e\n\n",freqnew(indexYNV),maxPSDyNV)
fprintf("Inner Tower Accelerometer Values\n")
fprintf("The Max Power in X occurs at %1.4f Hz w/ a value of %1.4e\n",freq(index),maxPSDx)
fprintf("The Max Power in Y occurs at %1.4f Hz w/ a value of %1.4e\n\n",freq(indexY),maxPSDy)
fprintf("Mars Platform Accelerometer Values\n")
fprintf("The Max Power in X occurs at %1.4f Hz w/ a value of %1.4e\n",freqP(indexXP),maxPSDxP)
fprintf("The Max Power in Y occurs at %1.4f Hz w/ a value of %1.4e\n\n",freqP(indexYP),maxPSDyP)

%%% Plotting Section
figure(1)
loglog(freqP,PSDxP)
hold on;
loglog(freq,PSDx)
xlabel("Freq (Hz)")
ylabel('|P1(f)|^2')
title("Power Specturm Density for X direction")
legend("Platform Acecl","Inner Tower Accel")
title("Accelerometer Comparison")
figure(2)
hold on;grid on;
plot(timeV,centroidxV,"-k")
plot(timeC,centroidxC,"-r")
xlabel("Time (sec)")
ylabel("X Position (arcsec)")
title("X Position of Centroid vs Time")
%%%%% Black is Video, Red is Calibration %%%%%
figure(3)
hold on;grid on;
plot(timeV,centroidyV,"-k")
plot(timeC,centroidyC,"-r")
%%%%% Black is Video, Red is Calibration %%%%%
xlabel("Time (sec)")
ylabel("Y Position (arcsec)")
title("Y Position of Centroid vs Time")
figure(4)
loglog(freqnew,PSDxNewV)
hold on;
loglog(freqC,PSDxNewC)
xlabel("Freq (Hz)")
ylabel('|P1(f)|^2')
title("Power Specturm Density Comparison in X Direction")
legend("Video","Calibration")
PSDxNewV=PSDxNewV-PSDxNewC;PSDxNewV(PSDxNewV<0)=0;
PSDyNewV=PSDyNewV-PSDyNewC;PSDyNewV(PSDyNewV<0)=0;
figure(5)
loglog(freqnew,PSDxNewV)
hold on;
loglog(freqP,PSDxP)
xlabel("Freq (Hz)")
ylabel('|P1(f)|^2')
title("Power Specturm Density for X direction")
legend("Video","Platform Accel")
figure(6)
loglog(freqnew,PSDyNewV)
hold on;
loglog(freqP,PSDyP)
xlabel("Freq (Hz)")
ylabel('|P1(f)|^2')
title("Power Specturm Density for Y direction")
legend("Video","Platform Accel")
figure(7)
hold on;
plot(freqnew,PSDxNewV)
plot(freqP,PSDxP)
plot(freqnew,PSDyNewV)
plot(freqP,PSDyP)
xlabel("Freq (Hz)")
ylabel('|P1(f)|^2')
title("Power Specturm Density")
legend("Video X","Platform Accel X","Video Y","Platform Accel Y")

index=find(freqnew>1,1,"first");
disp(sqrt(sum(PSDxNewV(index:end))))
disp(std(mean(centroidxV)))
disp(sqrt(sum(PSDyNewV(index:end))))
disp(std(mean(centroidyV)))
toc
function [PSD, freq] = time2PSD(Fs, disp)
    %FFT transform for Data
    N = length(disp);
    Y = fft(disp);
    P2 = abs(Y/N).^2;
    PSD = P2(1:floor(N/2)+1);%used floor ince N is odd
    PSD(2:end-1)=2*PSD(2:end-1);
    freq = 0:Fs/length(disp):Fs/2;
end