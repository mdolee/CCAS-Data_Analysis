clc; clear; close all;
tic
%%% Reading in Video File
start=30; %Starting time
v=VideoReader("test1.mov","CurrentTime",start);
FPS=round(v.FrameRate); %Frame Rate of Video
time=v.Duration-start; %sec, Duration of Video
maxFrame=round(FPS*time);

%%% Recreating Movie for Analysis
vidWidth = v.Width; %Pixels, Height of Video
vidHeight = v.Height; %Pixels, Width of Video
mov = struct('cdata',zeros(vidHeight,vidWidth,3,'uint8'),...
    'colormap',[]);
k = 1;
while hasFrame(v)
    if k<301
        mov(k).cdata = readFrame(v);
        k = k+1;
    else
        break
    end
end

%%% If you wanted to show the movie %%%
% hf = figure;
% set(hf,'position',[0 0 vidWidth vidHeight]);
% movie(hf,mov,1,FPS);

%%%% Creating Array of Coordinates for Hexagonal Grid %%%%
%Inupt Parameters
Ring=5; %Number of Rings around Center Ring

angle=0; %Angle of hex addition/Hex Rotation
len=2; %Keeping track of number of hex's
L=21.875; %non-dim, length of side
long=L*sqrt(3); %ratio btw side length and vertical height
center_x=-12; %Center of Grid in X dir (-7/-10)
center_y=-14; %Center of Grid in Y dir (-14/-16)

% Center Ring
vert=linspace(0+deg2rad(angle),2*pi+deg2rad(angle),7);hexSum=1;
for i=1:Ring
    hexSum=6*i+hexSum;
end
xvert=zeros(hexSum,7);yvert=zeros(hexSum,7);
xvert(1,:)=(L)*cos(vert)+center_x;yvert(1,:)=(L)*sin(vert)+center_y;
%Surrounding Rings
for j=1:Ring*2+1
    if j==1
        for i=1:2*Ring-floor(j/2)
            if i<Ring+j
                xvert(len,:)=(L)*cos(vert)+i*long*sind(angle)+center_x;
                yvert(len,:)=(L)*sin(vert)+i*long*cosd(angle)+center_y;
                len=len+1;
            else 
                xvert(len,:)=(L)*cos(vert)-(i-Ring)*long*sind(angle)+center_x;
                yvert(len,:)=(L)*sin(vert)-(i-Ring)*long*cosd(angle)+center_y;
                len=len+1;
            end
        end
    else
        for i=1:2*Ring-floor(j/2)+1
            if i==1
                const_x=-floor(j/2)*long*sind(angle+60)*(-1)^(j);
                const_y=-floor(j/2)*long*cosd(angle+60)*(-1)^(j);
                xvert(len,:)=(L)*cos(vert)+center_x+const_x;
                yvert(len,:)=(L)*sin(vert)+center_y+const_y;
                len=len+1;
            elseif i<=Ring-floor(j/2)+1 %Fill into SW
                xvert(len,:)=(L)*cos(vert)+(i-1)*long*sind(angle)*(-1)^(j-1)+center_x+const_x;
                yvert(len,:)=(L)*sin(vert)+(i-1)*long*cosd(angle)*(-1)^(j-1)+center_y+const_y;
                len=len+1;
            else  %Fill into NE
                xvert(len,:)=(L)*cos(vert)+(i-(Ring-floor(j/2)+1))*long*sind(angle)*(-1)^(j)+center_x+const_x;
                yvert(len,:)=(L)*sin(vert)+(i-(Ring-floor(j/2)+1))*long*cosd(angle)*(-1)^(j)+center_y+const_y;
                len=len+1;
            end
        end
    end
end
%%% Apply Rotation to Grid(If Necessary) %%%
theta=-3; %3, deg CCW (OG=3)
R=[cosd(theta) sind(theta);
    -sind(theta) cosd(theta);];
for j=1:len-1
    for i=1:7
        Pos=[xvert(j,i)-center_x; yvert(j,i)-center_y;];
        Pos_R=R*Pos;
        Pos_P=Pos_R+[center_x; center_y];
        xvert(j,i)=Pos_P(1);yvert(j,i)=Pos_P(2);
    end
end
%%% Apply Scaling to Grid(If Necessary) %%%
% scaleX=1;scaleY=scaleX;
% xvert=(xvert-center_x)*scaleX+center_x;
% yvert=(yvert-center_y)*scaleY+center_y;

%%% Center of Each Hexagon %%%
xcenter=xvert(:,1)-L*cosd(theta); %First point is at 0 deg, subtract L for center
ycenter=yvert(:,1)+L*sind(theta); %First point is at 0 deg, inline with center

%r^2 relationship to scaling factor of grid to dots
[xvert, yvert]=Inv2(center_x,center_y,xvert,yvert,L);
[xcenter, ycenter]=Inv2(center_x,center_y,xcenter,ycenter,L);

% %%% Match Hex Grid to Hex Grid %%%
% framenum=1;
% C=struct2cell(mov(framenum));
% A=cell2mat(C(1));
% frame=A(:,:,2);  %green is page 2
% figure()
% image(frame,'Xdata',[1 size(frame,2)]-vidWidth/2,'Ydata',[1 size(frame,1)]-vidHeight/2)
% hold on;colorbar;
% for k=1:len-1
%     plot(xvert(k,:),yvert(k,:),"-r")
% %     plot(xcenter(k),ycenter(k),"r*")
% end
% xlabel("Video Width (Pixels)");ylabel("Video Height (Pixels)");title("Hex Grid Matching")

%%% Compare Dots to Hex Grid %%%
framenum=1;
C=struct2cell(mov(framenum));
A=cell2mat(C(1));
frame=A(:,:,2);  %green is page 2
% frame(frame<40)=0; %Filter out visual noise
figure()
image(frame,'Xdata',[1 size(frame,2)]-vidWidth/2,'Ydata',[1 size(frame,1)]-vidHeight/2)
centroidx=zeros(1,hexSum);centroidy=zeros(1,hexSum);
hold on;colorbar;
x=1:size(frame,2);y=1:size(frame,1);
[X,Y]=meshgrid(x-vidWidth/2,y-vidHeight/2);
grid=im2double(frame); %Convert Image to Double for Calculation
for k=1:len-1
    plot(xvert(k,:),yvert(k,:),"-r")
    in=inpolygon(X,Y,xvert(k,:),yvert(k,:));
    maxVal=max(grid(in));
    index=find(in==1);
    filtered=find(grid(index)>maxVal*0.5);
    centroidx(k)=sum(grid(index(filtered)).*X(index(filtered)))/sum(grid(index(filtered)));
    centroidy(k)=sum(grid(index(filtered)).*Y(index(filtered)))/sum(grid(index(filtered)));
%     dist(k)=sqrt((centroidx(k)-xcenter(k))^2+(centroidy(k)-ycenter(k))^2);
    plot(centroidx(k), centroidy(k),"*r")
    plot(xcenter(k),ycenter(k),"k*")
end
clear centroidx centroidy
xlabel("Video Width (Pixels)");ylabel("Video Height (Pixels)");title("Centroid of Light in Hex Grid")

%%% Moving Image Frame %%%
filename = 'CCAS Video.gif';
f=figure(3);
scaleW=1;scaleH=1; %Scale image by factors
vidWidth=vidWidth*scaleW;vidHeight=vidHeight*scaleH;
framenum1=0;framenum2=maxFrame;subt=framenum2-framenum1;
time=linspace(1/FPS,10000/FPS,10000);
centroidx=zeros(hexSum,length(mov));centroidy=zeros(hexSum,length(mov));
for i=1:length(mov)
    %%% Converting Frames into Image %%%
    C=struct2cell(mov(i));
    A=cell2mat(C(1));
    Newframe=A(:,:,2);  %green is page 2
%    Newframe = imresize(A(:,:,2), [vidHeight vidWidth], 'nearest');
    %%% Image Analysis %%%
    if i==1
        image(Newframe,'Xdata',[1 size(frame,2)]-vidWidth/2,'Ydata',[1 size(frame,1)]-vidHeight/2)
        axis image
        hold on;
        grid=im2double(Newframe); %Convert Image to Double for Calculation
        for k=1:len-1
            plot(xvert(k,:),yvert(k,:),"-r")
            in=inpolygon(X,Y,xvert(k,:),yvert(k,:));
            maxVal=max(grid(in));
            index=find(in==1);
            filtered=find(grid(index)>maxVal*0.5);
            centroidx(k,i)=sum(grid(index(filtered)).*X(index(filtered)))/sum(grid(index(filtered)));
            centroidy(k,i)=sum(grid(index(filtered)).*Y(index(filtered)))/sum(grid(index(filtered)));
%             dist(k,i)=sqrt((centroidx(k,i)-xcenter(k))^2+(centroidy(k,i)-ycenter(k))^2);
            plot(centroidx(k,i), centroidy(k,i),"*r")
        end
        title("Image Animation")
        xlabel("Image Width (pixels)")
        ylabel("Image Height (pixels)")
        drawnow();
        frameP = getframe(f);
        im = frame2im(frameP); 
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0); 
    elseif i<300
        image(Newframe,'Xdata',[1 size(frame,2)]-vidWidth/2,'Ydata',[1 size(frame,1)]-vidHeight/2)
        axis image
        hold on;
        grid=im2double(Newframe); %Convert Image to Double for Calculation
        for k=1:len-1
            plot(xvert(k,:),yvert(k,:),"-r")
            in=inpolygon(X,Y,xvert(k,:),yvert(k,:));
            maxVal=max(grid(in));
            index=find(in==1);
            filtered=find(grid(index)>maxVal*0.5);
            centroidx(k,i)=sum(grid(index(filtered)).*X(index(filtered)))/sum(grid(index(filtered)));
            centroidy(k,i)=sum(grid(index(filtered)).*Y(index(filtered)))/sum(grid(index(filtered)));
%             dist(k,i)=sqrt((centroidx(k,i)-xcenter(k))^2+(centroidy(k,i)-ycenter(k))^2);
            plot(centroidx(k,i), centroidy(k,i),"*r")
        end
        drawnow();
        frameP = getframe(f); 
        im = frame2im(frameP); 
        [imind,cm] = rgb2ind(im,256);
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);
    else
        grid=im2double(Newframe); %Convert Image to Double for Calculation
        for k=1:len-1
            in=inpolygon(X,Y,xvert(k,:),yvert(k,:));
            maxVal=max(grid(in));
            index=find(in==1);
            filtered=find(grid(index)>maxVal*0.5);
            centroidx(k,i)=sum(grid(index(filtered)).*X(index(filtered)))/sum(grid(index(filtered)));
            centroidy(k,i)=sum(grid(index(filtered)).*Y(index(filtered)))/sum(grid(index(filtered)));
%             dist(k,i)=sqrt((centroidx(k,i)-xcenter(k))^2+(centroidy(k,i)-ycenter(k))^2);
        end
    end
end
Pix2arc=12.5/long; %arcsec, side to side dist is 12.5 arcsec, long is pixel dist from s2s
% %%% Plotting of X/YMotion and Dist to Center vs Time %%%
% figure()
% hold on;
for i=1:len-1
    centroidx(i,:)=(centroidx(i,:)-mean(centroidx(i,:)))*Pix2arc;
%     plot(time,centroidx(i,:)) %X Motion
end
% xlabel("Time (Sec)");ylabel("X Position (Arcsec)");title("X Motion vs Time")
% figure()
% hold on;
for i=1:len-1
    centroidy(i,:)=(centroidy(i,:)-mean(centroidy(i,:)))*Pix2arc;
%     plot(time,centroidy) %X Motion
end
% xlabel("Time (Sec)");ylabel("Y Position (Arcsec)");title("Y Motion vc   s Time")
% figure()
% hold on;
% for i=1:len-1
%     plot(time,(dist(i,:)-mean(dist(i,:)))*Pix2arc) %X Motion
% end
% xlabel("Time (Sec)");ylabel("Distance to Center (Arcsec)");title("Distance to Center vs Time")
%%% Fourier Transform Conversion %%%
for i=1:len-1
    [tempPSDx, ~] = time2PSD(FPS, centroidx(i,:)); %motion in Arcsec
    [tempPSDy, freq] = time2PSD(FPS, centroidy(i,:));
    tempPSDx(tempPSDx>0.1)=0; %Filtering out jump from camera
    tempPSDy(tempPSDy>0.1)=0; %Filtering out jump from camera
    if i==1
        PSDx=tempPSDx;
        PSDy=tempPSDy;
    else
        PSDx=((i-1).*PSDx+tempPSDx)./i;
        PSDy=((i-1).*PSDy+tempPSDy)./i;
    end
end
index1=find(freq>1,1,"first");xsum=sum(PSDx(index1:end));ysum=sum(PSDy(index1:end));
%%% Plotting
figure()
loglog(freq, PSDx) %X Motion
xlabel("Freq (Hz)");ylabel("P1(f)^2");title("Power Specturm Density of X Motion")
figure()
loglog(freq, PSDy) %Y 
xlabel("Freq (Hz)");ylabel("P1(f)^2");title("Power Specturm Density of Y Motion")
figure()
loglog(freq, PSDx);hold on;loglog(freq, PSDy);legend("X Motion","Y Motion")
xlabel("Freq (Hz)");ylabel("P1(f)^2");title("Power Specturm Density of X&Y Motion")
%Save data to a matfile
% save('Calibration2.mat','FPS','PSDx','PSDy','freq','centroidx','centroidy',"time");
%%%%%% Dimensions of Previous Videos %%%%%%
%ScaleW: 0.875 (testcam2) 1 (testCC)
%L: 14.5 (testcam2) 11 (testCC)
%centerX: -8 (testcam2) -3 (testCC)
%centerY: -14 (testcam2) -6 (testCC)
%theta: -3 (testcam2) -3 (testCC)
%curvature function: [1/1300 1/1800 1] (testcam2) [1/900 1/1800 1] (testCC)
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
function [xcord, ycord]= Inv2(center_x,center_y,xcord,ycord,L)
    dist=sqrt(((xcord-center_x).^2)+((ycord-center_y).^2))/L;
    scale=-(dist.^2)/1000+dist/1800+1;
    xcord=xcord.*scale;ycord=ycord.*scale;
end