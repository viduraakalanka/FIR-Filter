
clear all
close all
clc

%--------------------------------------------------------------------------
%Initializing FIR bandpass filter variables 
A = 0;
B = 1;
C = 8;

%--------------------------------------------------------------------------
%Defining filter variables
ApBar   = 0.1 +0.01*A ;
AaBar   = 50 + B;
OmegaS  = 2*(100*C +1000);
OmegaP1 = 300 +100*C;
OmegaP2 = 600 +100*C;
OmegaA1 = 200 +100*C;
OmegaA2 = 750 +100*C;

%--------------------------------------------------------------------------
%calculating actual value of Ap
deltaP = (10^(0.05*ApBar)-1)/(10^(0.05*ApBar)+1);
deltaA = 10^(-0.05*AaBar);
delta  = min(deltaP,deltaA);
Aa     = -20*log10(delta);

%--------------------------------------------------------------------------
%Calculating value of alpha
if Aa <= 21 
    alpha = 0;
elseif  21 < Aa && Aa <= 50
    alpha = 0.5842*(Aa -21)^0.4-0.07886*(Aa -21);
else
    alpha = 0.1102*(Aa-8.7);
end

%--------------------------------------------------------------------------
%Calculating parameter D
if Aa <= 21 
    D = 0.9222;
else
    D = (Aa - 7.95)/14.36;
end

%--------------------------------------------------------------------------
%Calculating N

Bt      = min((OmegaP1-OmegaA1),(OmegaA2-OmegaP2));
Napprox = OmegaS*D/Bt +1;
val     = rem(ceil(Napprox), 2);

if val == 1
    N = ceil(Napprox);
else
    N = ceil(Napprox)+1;
end
 
%--------------------------------------------------------------------------
%Determining window function 

T    = 2*pi/OmegaS;
n    = floor(-(N-1)/2):1:floor((N-1)/2);
Beta = alpha *sqrt(1-(2*n/(N-1)).^2);

k       = 1;
limit   = 10^-20;
fact    = 1;
Io_alpha = 0;
 
while 1
    fact = fact*k ; 
    val  = (((alpha/2)^k)/fact)^2;
    if  val>= limit
        Io_alpha = Io_alpha +val;
        k=k+1;
    else
        break
    end    
end

Io_alpha = Io_alpha +1;


Io_beta = ones(1,length(n));
for j  = 1: length(Io_beta)
    k    = 1 ;
    fact = 1 ;
    while true
        fact = fact*k;
        val  = ((((Beta(j))/2)^k)/fact)^2;
        if  val>= limit
            Io_beta(j) = Io_beta(j) +val;
            k=k+1;
        else
            break
        end     
    end
end
figure 
wk = Io_beta/Io_alpha;
stem(n*T,wk,'blue');
title('Kaiser Window Function')
xlabel('Time(seconds)')
ylabel('w(nT)')
grid on


Omega = OmegaS*n/(N-1);

%--------------------------------------------------------------------------
%Determining  frequency response of the Ideal Bandpass filter

OmegaC1 = OmegaP1 - Bt/2;
OmegaC2 = OmegaP2 + Bt/2;


H_freq  = zeros(1,length(n));

for j = 1:length(n)
    if Omega(j)<= (-OmegaC1) && Omega(j) >= (-OmegaC2)
        H_freq(j) =1;
    elseif  (Omega(j))>= OmegaC1 && Omega(j) <= (OmegaC2)
         H_freq(j) =1;
    end
end

%--------------------------------------------------------------------------
%Determining time-domain  response of the Ideal bandpass filter
len = length(n);
h_time = zeros(1,len);
for i = n
    if i ==0
        h_time(i+(len+1)/2) =2*(OmegaC2-OmegaC1)/OmegaS;
    else
        j = i+(len+1)/2;
        h_time(j) =(sin(OmegaC2*n(j)*T)-sin(OmegaC1*n(j)*T))/(n(j)*pi);
    end
end

figure
stem(n,h_time)
title('Impulse Response'),xlabel('n'),ylabel('h(n)')
grid on

%--------------------------------------------------------------------------
%Plotting impulse response of the digital filter
h_FIR = h_time.*wk;
[hOmega,w1] = freqz(h_FIR);

%--------------------------------------------------------------------------
%Determining FIR filter magnitude and amplitude response in the frequency domain

figure
plot(w1/T,20*log(abs(hOmega)))
axis([OmegaP1 OmegaP2 -ApBar/2 ApBar/2])
title('Amplitude Response of the passband'),xlabel('Frequency(rad/s)'),ylabel('H(e(jwT)) FIR')
grid on

figure
plot(w1/T,20*log10(abs(hOmega)))
title('Magnitude Response of the digital filter'),xlabel('Frequency(rad/s)'),ylabel('|H(e(jwT))|')
grid on

%--------------------------------------------------------------------------
%Determining x(n)

Omega1 = OmegaA1/2;
Omega2 = (OmegaP1+OmegaP2)/2;
Omega3 = (OmegaA2+OmegaS/2)/2;
n_new= 1:500;
x_n = sin(Omega1*n_new*T)+sin(Omega2*n_new*T)+sin(Omega3*n_new*T);
X_Omega=fft(x_n,OmegaS/2);

%--------------------------------------------------------------------------
% Frequency Response and Time response of the Input Signal

figure
stem(n_new,x_n)
title('Signal x(n)'),xlabel('n'),ylabel('x(n)')

figure
plot(0:2:1800,abs(X_Omega(1:901)))
title('Frequency Response of the Input Signal')
ylabel('Magnitude(dB)')
xlabel('Frequency(rad/s)')
grid on

%--------------------------------------------------------------------------
% Compute the Output of the Ideal and FIR Bandpass Filter
H_Fil=fft(h_FIR,OmegaS/2);
Y_Omega=X_Omega.*H_Fil;
y_n=ifft(Y_Omega);

H_Fil=fft(h_time,OmegaS/2);
Y_Out=X_Omega.*H_Fil;
y_out=ifft(Y_Out);

%--------------------------------------------------------------------------
% Plot of Ideal Bandpass Filter Output and its frequency response
figure
stem(0:500,y_out(1:501))
title('Output of a Ideal Bandpass Filter'),xlabel('n'),ylabel('x_(n)*h(n)')
grid on

figure
plot(0:2:1800,abs(Y_Out(1:901)))
title('Frequency Response of the Output after the Ideal Bandpass Filter')
ylabel('Y(w)/dB')
xlabel(' Angular Frequency(rad/s)')
grid on

%--------------------------------------------------------------------------
%plot of FIR Bandpass Filter Output and its frequency response
figure
stem(0:1:500,y_n(1:501))
title('Output of the Desinged FIR Bandpass Filter')
xlabel('n'),ylabel('y(n)')
grid on

figure
plot(0:2:1800,abs(Y_Omega(1:901)))
title('Frequency Response of the Output after the  Desinged Bandpass Filter')
ylabel('Y(w)/dB'),xlabel('Angular Frequency(rad/s)')
grid on
