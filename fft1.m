clear all

fs = 150;
ts = 1/fs;
L = fs * 2;

t = (0:L-1)*ts;

mp = 10 * rand();
fp = 5;

yp(:,:,1) = mp * sin(2 * pi * fp * t);

zp(:,:,1) = fft(yp(:,:,1));

f = fs * (0:(L/2))/L;
P2 = abs(zp(:,:,1)/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2 * P1(2:end-1);


figure
plot(f, P1)