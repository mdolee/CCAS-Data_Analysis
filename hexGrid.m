function [xvert, yvert, xcenter, ycenter, len, height, width] = hexGrid(Ring, angle, L, center_x, center_y, theta, warp)
%hexGrid
%   [xhex], [yhex], [xcentroid], [ycentroid], [hex#] = (ring#), ( 




len=2; %Keeping track of number of hex's
long=L*sqrt(3); %ratio btw side length and vertical height
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
if warp == true
[xvert, yvert]=Inv2(center_x,center_y,xvert,yvert,L);
[xcenter, ycenter]=Inv2(center_x,center_y,xcenter,ycenter,L);
end

height = (2*Ring + 1)*2*L*sind(60);
width = 2*((2*cosd(60) + L) + (Ring-1)*(cosd(60) + L));


len = len - 1;
end
