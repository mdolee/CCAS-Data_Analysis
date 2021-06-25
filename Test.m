clc; clear; close all;

L = .577; %meters **check
pixdensity = 20; 

[xvert, yvert, xcenter, ycenter, len, h, w] = hexGrid(5, 0, L, 0, 0, 0, false);

figure(1)
hold on;
drawnow();
axis('equal')
for i = 1:len
    
    plot(xvert(i,:),yvert(i,:), 'b')

end
plot(xcenter,ycenter, '*')


%where, or findall? 
 
sig = .1; %adjust, smaller
pix_y = ceil(pixdensity * h);
pix_x = ceil(pixdensity * w);


%meshgrid

[Xo,Yo] = meshgrid(linspace(-(w/2), w/2,pix_x), linspace(-(h/2),(h/2),pix_y));
% Zee = zeros(length(Xo),width(Xo));
for i = 1:len

[in, on] = inpolygon(Xo,Yo,xvert(i,:),yvert(i,:));

[xind, yind] = find(in);
 Xnew = zeros(length(xind),length(xind));
 Ynew = zeros(length(xind),length(xind));
 Zee = zeros(length(xind),length(xind));
tic
 for j = 1:length(xind)
    
    for k = 1:length(yind)
    Zee(j,k) = Gauss_2D(1, Xo(xind(j),yind(k)), xcenter(i), sig, Yo(xind(j),yind(k)), ycenter(i), sig);
    Xnew(j,k) = Xo(xind(j),yind(k));
    Ynew(j,k) = Yo(xind(j), yind(k));
    end
   
 end
toc
drawnow();
tic
 h = surf(Xnew,Ynew,Zee);
 toc
 
 set(h,'LineStyle','none')
 disp(i)
end


%pcolor

