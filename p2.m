clear all
close all

testN = 10;
testT = 20;

fs = 256;
ts = 1/fs;
L = fs * testT;
t = (0:L-1)*ts;

for m=1:testN
    F = 100*rand();
    y(:,:,m) = rand()*sin(2*pi*F*t);
    z(:,:,m) = fft(y(:,:,m));
    B(:,:,m) = y(:,:,m) + rand()*sqrt(0.5)*randn(size(t));%signal with gaussian noise
end


for i=1:3
    %(a)
    N = randi([2 5]);
    Tk = 2^randi([0 2]);
    for n=1:N
        mp = 0.5 + 0.5*rand();
        fp = (0.25 + 1.75*rand())/Tk;
        php = 2*pi*fp*rand();
        pl(:,:,n) = mp*sin(2*pi*fp*t + php);
    end
    a(:,:,i) = pl(:,:,1);
    for n=2:N
        a(:,:,i) = a(:,:,i) + pl(:,:,n);
    end
    
    p(:,:,i) = a(:,:,i);
    
    %(b)
    b(:,:,i) = sqrt(0.5)*randn(size(t));
    
    p(:,:,i) = p(:,:,i) + b(:,:,i); %add gaussian noise
    
    %(c)
    mp = 0.2 * 1.8*rand();
    fp = (0.25 + 1.75*rand())/Tk;
    php = 2*pi*fp*rand();
    c(:,:,i) = + mp*sawtooth(2*pi*fp*t + php);
    
    p(:,:,i) = p(:,:,i) + c(:,:,i);
    
    %(d)
    for n4=1:10
       mp =  0.1 + 0.9*rand();
       fp = (20 + 80*rand())/Tk;
       dk(:,:,i,n4) = mp * sin(2*pi*fp*t);
       if n4 == 1
           d(:,:,i) = dk(:,:,i,n4);
       else
           d(:,:,i) = d(:,:,i) + dk(:,:,i,n4);
       end
    end
     p(:,:,i) = p(:,:,i) + d(:,:,i);
end

for n=1:testN
    tt = 0;
    while tt <= testT 
        mod = randi([1 4]);
        if mod == 4
            Tk = 3 * rand();
        else
            Tk = 2^randi([0 2]);
        end
        if Tk + tt >= testT
            break
        end
        tt = Tk + tt;
        if mod ~= 4
            %calcing a
            N = randi([2 5]);
            for n1=1:N
                mp = 0.5 + 0.5*rand();
                fp = (0.25 + 1.75*rand())/Tk;
                php = 2*pi*fp*rand();
                pl(:,:,n1) = mp*sin(2*pi*fp*t + php);
            end
            a(:,:,mod) = pl(:,:,1);
            for n1=2:N
                a(:,:,mod) = a(:,:,mod) + pl(:,:,n1);
            end
            %calcing b
            b(:,:,mod) = sqrt(0.5)*randn(size(t));
            B(:,:,n) = B(:,:,n) + (a(:,:,mod) + b(:,:,mod) + c(:,:,mod) + d(:,:,mod)).*(t>=tt & t< Tk+tt);
        end
        
    end
end

for n=1:testN
    
%     pRMS = rms(B(:,:,n))^2;
%     disp(pRMS);
    figure
    plot(t, B(:,:,n))

end

for n=1:3
    %fft
    Y = fft(d(:,:,n));
    f = fs * (0:(L/2))/L;
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);   
    P1(2:end-1) = 2 * P1(2:end-1);
 
    figure
    plot(f, P1);
end

