function [G2d] = Gauss_2D(A, X, xo, sig_x, Y, yo, sig_y)
%GAUSS_2D Creates a 2 dimensional gaussian distribution
%   [2D Gaussian Distribution(z)] = (Amplitude), (Xinitial), (SigmaX), (Yinitial), (SigmaY)  




G2d = A.*exp(-(((X-xo).^2)./(2.*sig_x.^2)+((Y-yo).^2)./(2.*sig_y.^2)));

end

