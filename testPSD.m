% Testing PSD against sine wave
N = 10000;
Tdur = 10.0;
Arms = 0.126e-3;
t = linspace(0,Tdur,N);
y = Arms * sqrt(2) * sin(2*pi/5.*t);
f = linspace(0,1/(t(2)-t(1))/2,N/2);
df = f(2)-f(1);
Y = abs(fft(y)/N).^2;
psd = Y(1:N/2) * 2;
loglog(f,psd)

%%
t = linspace(0, 10,10001);
y = sin(2*pi*t/5) + randn(size(t))/10;
sgf = sgolayfilt(y, 3, 101);
figure
plot(t, y)
hold on
plot(t, sgf, '-r')
hold off
grid
L = numel(t);
Ts = mean(diff(t));
Fs = 1/Ts;
Fn = Fs/2;
FTy = fft(y)/L;
FTsgf = fft(sgf)/L;
Fv = linspace(0, 1, fix(L/2)+1)*Fn;
Iv = 1:numel(Fv);
figure
subplot(3,1,1)
plot(Fv, abs(FTy(Iv))*2)
grid
title('Fourier Transform: y')
subplot(3,1,2)
plot(Fv, abs(FTsgf(Iv))*2)
grid
title('Fourier Transform: Savitzky-Golay Filter of y')
subplot(3,1,3)
loglog(Fv, abs(FTsgf(Iv)./FTy(Iv))*2)
grid
title('Savitzky-Golay Transfer Function')