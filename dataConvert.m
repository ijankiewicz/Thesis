clc
clear

# NOTE: nice results for dataToRead  = csvread('ltc3588_500_1.csv');

# load singal package required for data processing
pkg load signal

# accelerometer parameters
accelOffset = 1710;       # offset of the Z-axis output of accelerometer
accelSenitivity = 0.33;   # V/g
accelGain = 0.1;          # gain of the channel set on the logger board

# ADC parameters
adcResolution = 2^14;     # resolution of the logger board
voltageReference = 3.00;  # reference voltage provided for the ADC

# converter output parameters
voltGain = 1;	# gain of the channel set on the logger board
load = 5000;	# value of the load connectoed to the converter's out

# serial data acquisition channel paraemters
serialBaudrate        = 460800;   # serial baud rate
serialChnCount        = 2;        # current channel count
serialBytesPerSample  = 2;        # number of bytes per one sample
serialSamplePerSecond = serialBaudrate / ...
                        (8 * serialBytesPerSample * serialChnCount);
                                  # calculate number of samples per second
sampleTimeBase        = 1 /serialSamplePerSecond;
                  	# calculate time between consecutive samples

# data import
dataToRead  = csvread("ltc3588_5k.csv");     
                                              # select .csv file
dataNegOffset = 10;                           # rm samples at the beginning
dataPosOffset = 10;                           # rm samples at the end
voltRaw     = dataToRead(:,1);                # import raw - 1st channel
accelRaw    = dataToRead(:,2) - accelOffset;  # import raw - 2nd channel
sampleCount = size(voltRaw)(:,1);             # count raw samples

voltRawOffset  = dataToRead(dataNegOffset:(sampleCount - dataPosOffset),1);
             	# adjust the number of samples from 1st channel
accelRawOffset = dataToRead(dataNegOffset:(sampleCount ...
                    - dataPosOffset),2) - accelOffset;
             	# adjust the number of samples from 2nd channel
sampleCountOffset = size(voltRawOffset)(:,1); 
		# count samples after adjustment

# create time base for samples
sampleTimeDomain = zeros(sampleCountOffset,1);  # create emty array
for i = 1:sampleCountOffset
  sampleTimeDomain(i) = i * sampleTimeBase;     # create time points
endfor  


# convert raw data into proper values
voltConv    = voltRawOffset * voltageReference * voltGain / adcResolution;
                                  # calculate output voltage
currentConv = voltConv / load;    # calculate output current
powerConv   = voltConv .* currentConv;
                                  # calculate output power
accelConv   = accelRawOffset * voltageReference * ... 
              (1/accelGain) * accelSenitivity / adcResolution;
                                  # calculate acceletration (in g)

# calculate charge              
charge = 0;              
for j = 1:(sampleCountOffset)
  charge = (currentConv(i) * sampleTimeBase) + charge; 
                                  # use simple current * time relation
endfor

# band pass filter for accelerometer readings
accelLowPassOrder   = 5;
accelHighPassOrder  = 5;
accelLowPassFc      = 35;
accelHighPassFc     = 100;
           
[b,a]=butter( accelHighPassOrder, ...
              accelLowPassFc/(serialSamplePerSecond/2), ...
              'high');
accelPreaccelFiltered = filter(b,a,accelConv);

[d,c]=butter( accelLowPassOrder, ...
              accelHighPassFc/(serialSamplePerSecond/2), ...
              'low');
accelFiltered = filter(d,c,accelPreaccelFiltered);


# calculate acceleration associated with a hit
accelAvg = (abs(min(accelFiltered)) + max(accelFiltered)) / 2;


# generate FFT for acceletometer readings
Fs = serialSamplePerSecond;
Ts = 1/Fs;
L = sampleCountOffset;
t = (0:L-1)*Ts;
Y = fft(accelConv);
P2 = abs(Y/L);
P1 = P2(1:L/2);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(1:(L/2))/L;

Y2 = fft(accelFiltered);
P22 = abs(Y2/L);
P12 = P22(1:L/2);
P12(2:end-1) = 2*P12(2:end-1);
f = Fs*(1:(L/2))/L;


# find frequency peaks
[peak peakLocation] = findpeaks(P1,  "DoubleSided",           ...
                                     "MinPeakHeight",   0.01,  ...
                                     "MinPeakDistance", 10   );

for k = 1:size(peakLocation)(:,1)
    if (f(peakLocation(k)) > 30 && f(peakLocation(k)) < 100)
      vibrationFrequency = f(peakLocation(k));
    endif
endfor  

                               
figure(1)
h1 = figure(1);
set (h1,'papertype', '<custom>')
set (h1,'paperunits','inches');
set (h1,'papersize',[9 5.5])
set (h1,'paperposition', [0,0,[8 5.5]])
set (h1,'defaultaxesposition', [0.15, 0.15, 0.75, 0.75])
set (0,'defaultaxesfontsize', 10)

subplot(3,1,1)

[ax11 hline111 hline112] = plotyy(sampleTimeDomain, accelConv ...
                                , sampleTimeDomain, voltConv);
set(hline112, 'linewidth', 1);
title('Acceleration and output voltage in time domain');
xlabel('time [s]');
ylabel(ax11(1), 'Acceleration [g]');
ylabel(ax11(2), 'Output voltage [V]');
l11 = legend ('Acceleration','Output voltage');
set (l11, 'interpreter', 'tex', 'fontsize', 6, "location", "northeast" ...
     , "interpreter", "latex");
grid on;

subplot(3,1,2)
[ax12 hline121 hline122] = plotyy(sampleTimeDomain, accelConv ...
                                , sampleTimeDomain, currentConv * 1000);
set(hline122, 'linewidth', 1);
title('Acceleration and output current in time domain');
xlabel('time [s]');
ylabel(ax12(1), 'Acceleration [g]');
ylabel(ax12(2), 'Output current [mA]');
l12 = legend ('Acceleration','Output current');
set (l12, 'interpreter', 'tex', 'fontsize', 6, "location", "northeast" ...
     , "interpreter", "latex");
grid on;

subplot(3,1,3)
[ax13 hline131 hline132] = plotyy(sampleTimeDomain, accelConv ...
                                , sampleTimeDomain, powerConv * 1000);
set(hline132, 'linewidth', 1);
title('Acceleration and output power in time domain');
xlabel('time [s]');
ylabel(ax13(1), 'Acceleration [g]');
ylabel(ax13(2), 'Output power [mW]');
l13 = legend ('Acceleration','Output power');
set (l13, 'interpreter', 'tex', 'fontsize', 6, "location", "northeast" ...
     , "interpreter", "latex");
grid on;

print('timedomain.pdf','-dpdf')


figure(2)
h2 = figure(2);
set (h2,'papertype', '<custom>')
set (h2,'paperunits','inches');
set (h2,'papersize',[8 4.5])
set (h2,'paperposition', [0,0,[8 4.5]])
set (h2,'defaultaxesposition', [0.15, 0.15, 0.75, 0.75])
set (0,'defaultaxesfontsize', 10)


subplot(2,2,1)
plot(sampleTimeDomain, accelConv)
title('The unfiltered accelerometer sginal')
grid minor;
xlabel('time [s]')
ylabel('acceletation [g]')

subplot(2,2,2)
plot(sampleTimeDomain, accelFiltered)
title('The filtered accelerometer sginal')
grid minor;
xlabel('time [s]')
ylabel('acceletation [g]')

subplot(2,2,3)
semilogx(f,P1, "k") 
grid minor;
title('FFT of the unfiltered accelerometer sginal')
xlabel('f[Hz]')
ylabel('|P(f)|')


subplot(2,2,4)
semilogx(f,P12, "k")
grid minor;
title('FFT of the filtered accelerometer sginal')
xlabel('f[Hz]')
ylabel('|P(f)|')
ylim([0 0.12]);

print('spectrum.pdf','-dpdf')


# RESULTS
printf("Generated charge:    %.3fuC\r\n", charge*1e6);
printf("Vibration frequency: %.2fHz\r\n", vibrationFrequency);
printf("Accelaration:        %.3fg\r\n", accelAvg);