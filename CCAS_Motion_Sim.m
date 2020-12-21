clc;clear; close all;
format long g 
tic
%Load in PSD data files
m = matfile('PSD.mat');Fs=m.Fs;
PSDx = m.PSDx;freqx = m.freqx';
PSDy = m.PSDy;freqy = m.freqy';
PSDz = m.PSDz;freqz = m.freqz';
%Inverse FFT from Freq domain to Time Domain
min=30; %Length of simulated signal
[DispX_sig,t_sig]=invFFT(min,PSDx,Fs,freqx);
[DispY_sig,~]=invFFT(min,PSDy,Fs,freqy);
[DispZ_sig,~]=invFFT(min,PSDz,Fs,freqz);

%Spot Motion Sim
cordx=0;cordy=0;
for i=1:length(DispX_sig)
    newX=cordx(i)+DispX_sig(i);
    newY=cordy(i)+DispY_sig(i);
    cordx=[cordx; newX];cordy=[cordy; newY];
end
figure()
hold on;grid on;
plot(cordx,cordy)
plot(cordx(1),cordy(1),"ro")
xlabel("X Cordinates (m)")
ylabel("Y Cordinates (m)")
title("X-Y Motion of Dot")

%%% Plotting Section %%%
% figure()
% plot(t_sig,DispX_sig)
% xlabel("Time(sec)")
% ylabel("Disp (m)")
% title("Mock Up Signal for X from PSD Integration")
% figure()
% plot(t_sig,DispY_sig)
% xlabel("Time(sec)")
% ylabel("Disp (m)")
% title("Mock Up Signal for Y from PSD Integration")
% figure()
% plot(t_sig,DispZ_sig)
% xlabel("Time(sec)")
% ylabel("Disp (m)")
% title("Mock Up Signal for Z from PSD Integration")

%Check for how accurate Disp Signal is of PSD
fprintf("Integrated Power of PSDx: %1.8f microns\n",sqrt(trapz(PSDx))*1e6)
fprintf("Standard Dev of Created Disp x: %1.8f microns\n",sqrt(std(DispX_sig).^2)*1e6)
fprintf("Integrated Power of PSDy: %1.8f microns\n",sqrt(trapz(PSDy))*1e6)
fprintf("Standard Dev of Created Disp y: %1.8f microns\n",sqrt(std(DispY_sig).^2)*1e6)
fprintf("Integrated Power of PSDz: %1.8f microns\n",sqrt(trapz(PSDz))*1e6)
fprintf("Standard Dev of Created Disp z: %1.8f microns\n",sqrt(std(DispZ_sig).^2)*1e6)
toc
function [Disp_sig,t_sig] = invFFT(min,PSD,Fs,freq) %freq
    N=round(min*60*Fs,0);
    order = 3; width = 51; %play with?
    newX=linspace(0,Fs/2,N)';
    A=sqrt(interp1(freq,sgolayfilt(PSD,order,width),newX));
    Phase_x=2*pi*rand(N,1);
    fftx_signal=A.*exp(1i.*Phase_x);
    Disp_sig=real(ifft(fftx_signal))*N; %has imaginary parts but very small, ignore?
    t_sig=0:1/Fs:1/Fs*(N-1);
end