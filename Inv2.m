function [xcord, ycord, scale, dist]= Inv2(center_x,center_y,xcord,ycord,L)
%INV2 warps grid to match camera focus

    dist = sqrt(((xcord-center_x).^2) + ((ycord-center_y).^2))/ L; %pythag 
    scale= -(dist.^2) / 1000 + dist / 1800 + 1; %scale factor = -(r^2)/1000 + r/1800 + 1  
    xcord=xcord.*scale;
    ycord=ycord.*scale;
end

