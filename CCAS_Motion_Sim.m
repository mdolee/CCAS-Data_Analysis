clc;clear; close all;
format long g 
tic
%Load in PSD data files
m = matfile('PSD.mat');Fs=m.Fs;
PSDx = m.PSDx;freqx = m.freqx';
PSDy = m.PSDy;freqy = m.freqy';
PSDz = m.PSDz;freqz = m.freqz';
%Inverse FFT from Freq domain to Time Domain
minute=30*60/60; %Length of simulated signal
[DispX_sig,t_sig]=invFFT(minute,PSDx,Fs,freqx);
[DispY_sig,~]=invFFT(minute,PSDy,Fs,freqy);
[DispZ_sig,~]=invFFT(minute,PSDz,Fs,freqz);

%Spot Motion Sim
sr=0.02; %sec, sampling rate of accelerometer, 20 ms
exposure=1; %seconds
cordx=DispX_sig/125e-6;cordy=DispY_sig/125e-6; %in arcsec, starting spot value
cordx=reshape(cordx,1/sr*exposure,minute*60/exposure);cordy=reshape(cordy,1/sr*exposure,minute*60/exposure);
figure()
hold on;grid on;
plot(cordx,cordy,"o")
xlabel("X Cordinates (arcsec)")
ylabel("Y Cordinates (arsec)")
title("X-Y Motion of Dot")

%Radial Distance Code
%pixel=0.044; %arcsec per pixel size
%r=sqrt(cordx.^2+cordy.^2); %center is 0,0
%f=histogram(r);
%y=0:pixel:1.1*max(r); %bin sizes
%f.BinEdges=y;
%hold on;
%mu=mean(r); %average of radii
%sigma=std(r); %standard deviation of radii
%func=exp(-(y-mu).^2./(2*sigma^2))*max(f.Values);%./(sigma*sqrt(2*pi));
%plot(y,func,'LineWidth',1.5)
%xlabel("Radii Distance (arcsec)")
%ylabel("Bin Count")
%title("Histogram of Spot Motion Radii")

%Image Distortion Animation
pixel=0.044; %arcsec, 1 pixel size
sigma=pixel+pixel*rand(1);%random width of bell from 1 to 2 pixels
grid_s=80; %50 square pixel grid
[X,Y]=meshgrid(-grid_s/2*pixel:pixel:grid_s/2*pixel, -grid_s/2*pixel:pixel:grid_s/2*pixel);
A=1; %Amplitude
B=size(cordx); %Size of reshaped cordx matrix
filename = '60secTestNoisy.gif';
f=figure(2);
axis tight manual
Z=zeros(length(X),length(Y),B(2)); %meshgrid zero
for i=1:B(2)
    for j=1:B(1)
        xo=cordx(j,i);yo=cordy(j,i);%amplitude
        Gaus_2d=A*exp(-(((X-xo).^2)/(2*sigma^2)+((Y-yo).^2)/(2*sigma^2)));
        Z(:,:,i)=Z(:,:,i)+Gaus_2d;
    end
    %Intesnity-weighted centroid of image
    centroidx(i)=sum(Z(:,:,i).*X)/sum(Z(:,:,i));
    centroidy(i)=sum(Z(:,:,i).*Y)/sum(Z(:,:,i));
    if i==1
        g=pcolor(X,Y,Z(:,:,i));
        title("Image Animation")
        xlabel("X cord (arcsec)")
        ylabel("Y cord (arcsec)")
        set(g, "EdgeColor", "none")
        drawnow();
        frame = getframe(f); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
    else 
        set(g, "CData", Z(:,:,i));
        set(g, "EdgeColor", "none")
        drawnow();
        frame = getframe(f); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
    end
end
figure()
hold on;grid on;
plot(centroidx,centroidy,"o")
xlabel("Xcord (arcsec)")
ylabel("Ycord (arcesec)")
title("Motion of Image Centroid")
%Spot Motion Animation Code
% h = figure;
% axis tight manual % this ensures that getframe() returns a consistent size
% filename = '1secTest.gif';
% for i=1:length(cordx) %length of data is same for all
%     plot(cordx(i),cordy(i),"o")
%     title("X-Y Animation of Dot")
%     xlabel("X Cordinates (arcsec)")
%     ylabel("Y Cordinates (arsec)")
%     xlim([1.1*min(cordx) 1.1*max(cordx)])
%     ylim([1.1*min(cordy) 1.1*max(cordy)])
%     drawnow 
%     % Capture the plot as an image 
%     frame = getframe(h); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if i == 1 
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
%     else 
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0); 
%     end
% end

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