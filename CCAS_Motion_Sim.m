clc;clear; close all;
format long g 
tic
% Report of Findings
% Fix 30 min duration to carve time to smaller/larger time
% Pick rndom points until signal adds up to lenght of signal
% Record animation length for only first minute

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

%Coordinate System Transformation
%(raw coord system data is in x,y,z, roll 55 deg CW)
theta=55; %deg
DispXnew=DispX_sig;                                  %[1    0      0];
DispYnew=DispY_sig*cosd(theta)-DispZ_sig*sin(theta); %[0 cos(t) -sin(t)]
DispZnew=DispY_sig*sind(theta)+DispZ_sig*cosd(theta);%[0 sin(t) cos(t)]
%Spot Motion now in (xnew,znew) plane, ynew affects size of dot

%Spot Motion Sim
sr=0.02; %sec, sampling rate of accelerometer, 20 ms
exposure=1/10; %seconds
cordx=DispXnew;cordz=DispZnew; %arcsec, starting spot value
cordx=reshape(cordx,1/sr*exposure,minute*60/exposure);cordz=reshape(cordz,1/sr*exposure,minute*60/exposure);
dy=(DispYnew-mean(DispYnew))/125e-6; %arcsec, converting to arcsec 

figure(1)
hold on;grid on;
plot(cordx,cordz,"o")
xlabel("X Cordinates (arcsec)")
ylabel("Z Cordinates (arsec)")
title("X-Z Motion of Dot")

%Image Distortion Animation
pixel=0.044; %arcsec, 1 pixel size
sigma=pixel+0.42*dy/1e6;%random width of bell from 1 to 2 pixels
grid_s=80; %80 square pixel grid
[X,Z]=meshgrid(-grid_s/2*pixel:pixel:grid_s/2*pixel, -grid_s/2*pixel:pixel:grid_s/2*pixel);
A=1; %Amplitude
B=size(cordx); %Size of reshaped cordx matrix
s=num2str(exposure);str1="Noisy";
filename=append(str1," ",s,"second exposure.gif");
f=figure(2);
axis tight manual
Y=zeros(length(X),length(Z),B(2)); %meshgrid zero
for i=1:B(2)
    for j=1:B(1)
        xo=cordx(j,i);zo=cordz(j,i);%amplitude
        Gaus_2d=A*exp(-(((X-xo).^2)/(2*sigma(i)^2)+((Z-zo).^2)/(2*sigma(i)^2)));
        Y(:,:,i)=Y(:,:,i)+Gaus_2d;
    end
    %Intesnity-weighted centroid of image
    centroidx(i)=sum(Y(:,:,i).*X)/sum(Y(:,:,i));
    centroidz(i)=sum(Y(:,:,i).*Z)/sum(Y(:,:,i));
    if i==1
        g=pcolor(X,Z,Y(:,:,i));
        title("Image Animation")
        xlabel("X cord (arcsec)")
        ylabel("Z cord (arcsec)")
        set(g, "EdgeColor", "none")
        drawnow();
        frame = getframe(f); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
    elseif i<180
        set(g, "CData", Y(:,:,i));
        set(g, "EdgeColor", "none")
        drawnow();
        frame = getframe(f); 
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
    else
        continue
    end
end
figure(3)
hold on;grid on;
plot(centroidx,centroidz,"o")
xlabel("Xcord (arcsec)")
ylabel("Zcord (arcsec)")
title("Motion of Image Centroid, Noisy")

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
% plot(t_sig,DispXnew)
% xlabel("Time(sec)")
% ylabel("Disp (m)")
% title("Mock Up Signal for X from PSD Integration")
% figure()
% plot(t_sig,DispYnew)
% xlabel("Time(sec)")
% ylabel("Disp (m)")
% title("Mock Up Signal for Y from PSD Integration")
% figure()
% plot(t_sig,DispZnew)
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
