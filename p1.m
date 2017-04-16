
clear all

fs = 150;
ts = 1/fs;
L = fs * 2;

t = (0:L-1)*ts;

mp = 10 * rand();
fp = 5;

yp(:,:,1) = mp * sin(2 * pi * fp * t);

xp(:,:,1) = t.*(t>=0 & t<=1) + 2*t.*(t>1 & t<= 2);

gp(:,:,1) = yp(:,:,1) + xp(:,:,1);

%fft
zp(:,:,1) = fft(yp(:,:,1));

f = fs * (0:(L/2))/L;
P2 = abs(zp(:,:,1)/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2 * P1(2:end-1);


figure
% plot(f, P1)
plot(t, xp(:,:,1))
figure
plot(t, yp(:,:,1))
figure
plot(t, gp(:,:,1))



x=0:0.001:2;
y=x.*(x>=0 & x<1)+2*x.*(x>=1 & x<=2);
% figure
% plot(x,y)

