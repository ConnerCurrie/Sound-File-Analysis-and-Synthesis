%SYDE 252 Project 
%Part 1 - Decompose Singals 
clc
clear
clf
%Define file name
fname = 'conner_high.wav';
% return wav file as an array, x, and the sampling rate, Fs
[x, Fs] = audioread(fname); 
%Extract single column from stereo
x = x(:,1);
%Set the number of samples to analyse
sNum = length(x);
%Ensure even sample number
if (mod(sNum, 2) ~= 0)
    sNum = sNum -1;
end
%Creat a time vector:
t = (0:(length(x)-1))/Fs;
%Define the starting point of the sampling window:
t1 = 0;
%Sample to start at:
sNum1 = round(Fs*t1)+1;
%Take the fast fourier transform of x:
xdft = fft(x(sNum1:sNum1+sNum-1));
%Take one half of the frequency domain due to symmetry
plotxdft = xdft(1:sNum/2+1);  
%Calculate the phase in a range from -pi to +pi
phxdft = angle(plotxdft(1:sNum/2+1)); 
%Calculate the magnitude in dB
dbmagxdft = 20*log10(real(abs(plotxdft)));
%Array of frequencies
freq = 0:Fs/sNum:Fs/2;
%Create a figure window:
figure(1);
%Create the first subplot (2 rows, 1 column, 1st plot)
subplot(3,1,1);
% plot the time domain "waveform" of our signal
plot(t,x);
% to plot something else on top we need to tell the window to hold on
hold on;
%define axis limits for plotting to make it look nicer
title('Time Domain Waveform');
xlabel('Time (s)');
ylabel('Amplitude');

%Create the second subplot (3 rows, 1 column, 2nd plot)
subplot(3,1,2);
%plot(freq,10*log10(psdx));
plot(freq,dbmagxdft);
grid on;
%axis([20 Fs/2 -140 0]);
%Define the xlow xhigh ylow yhigh limits that the plot shows
axis([20 22050 -100 100]); 
title('Magnitude(X(jw))');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');

%Create the third subplot (3 rows, 1 column, 3rd plot)
subplot(3,1,3);
plot(freq,phxdft);
grid on;
%axis([20 Fs/2 -140 0]);
%Define the xlow xhigh ylow yhigh limits that the plot show
axis([20 5000 -pi pi]); 
title('Phase(X(jw))');
xlabel('Frequency (Hz)');
ylabel('Phase (rads)');

% Play back our file each time we run the script
player = audioplayer(x, Fs, 16);
play(player);
pause(6);

%Part 2 - Synthesize the sound ==========================================
[magArray indices] = sort(dbmagxdft);
%Create an array to store frequencies
topFreq = zeros(9,2);
%Create variable for spacing
currentFreq = 0;
i2 = 1;
%Search for highest magnitudes in the array
while (i2 <= 9)
    for i = length(magArray)/2:-1:1;
        if (i2 <= 9)
            if ((currentFreq - 50) < indices(i) && (currentFreq + 50) > indices(i))
            else
                topFreq(i2,1) = freq(indices(i));
                topFreq(i2, 2) = 10^(magArray(i)/20);
                i2 = i2 + 1;
                currentFreq = indices(i);
            end
        end
    end
end
%3 Cosines

%6 Cosines

%9 Cosines

% Define a sampling frequency, for audio use 44100 Hz 
Fs = 44100;
% set how long we want our sample sound to be in seconds
tTotal = 5;
% create an associated time vector
t = (0:1/Fs:tTotal);

% first frequency component in our signal
A1 = 1;
f1 = 1115;
ph1 = 0;
x1 = A1*cos(2*pi*f1*t+ph1);
% second frequency component in our signal
A2 = 0.251;
f2 = 2228;
ph2 = 0;
x2 = A2*cos(2*pi*f2*t+ph2);
% third frequency component in our signal
A3 = 0.057;
f3 = 3345;
ph3 = 0;
x3 = A3*cos(2*pi*f3*t+ph3);

xsum = x1+x2+x3;

% maximum time index when plotting
tlim = 5;
tattack = 0;
attack = 2;
decay = 2;

% define an envelope to apply to your signal
%env = ones(size(xsum));
%env = abs(sin(2*pi*0.2*t));
env = 0.9*[1/(tattack^attack)*t(t<tattack).^attack, 1/exp(1)^(-decay*tattack)*exp(1).^(-decay*t(t>=tattack))];

% multiply the signal by the envelope
x = xsum.*env; 

% Save the soundfile output
audiowrite('new_high_f.wav',x,Fs);

% Play back our file each time we run the script
player = audioplayer(x, Fs, 16);
play(player);










