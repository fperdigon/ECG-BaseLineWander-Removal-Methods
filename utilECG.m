function utilECG = utilECG
    utilECG.detectRPeaks_FDeriv = @detectRPeaks_FDeriv;
    utilECG.detectRPeaks = @detectRPeaks;
    utilECG.detectRPeaksTh = @detectRPeaksTh;
    utilECG.detectPQ = @detectPQ;
    utilECG.detectTP = @detectTP;
    utilECG.detectST = @detectST;
    utilECG.instVariance = @instVariance;
    utilECG.addBaseLine = @addBaseLine;
    utilECG.addLineNoise = @addLineNoise;
    utilECG.lowPassFilter = @lowPassFilter;
    utilECG.highPassFilter = @highPassFilter;
    utilECG.readECG_Discardia = @readECG_Discardia;
    utilECG.FFTplot = @FFTplot;
    utilECG.ECGx = @ECGx;
end

function [serieRRx,serieRRy] = detectRPeaks_FDeriv(ecgy)
%  Detection of R wave using the first derivate
%
%  ecgy:        ECG signal
%
%  serieRRx:    X axis (time) position of R waves
%  serieRRy:    Y axis (amplitude) of R waves
%
%  Reference:
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    Decg = diff(ecgy);
    [M,N]=size(ecgy);

    c= 1; x= 1; ymax= max(Decg); flag=0;
    y= ymax * 0.5;
    for i=2 : M-1
        if Decg(i) > y
            y = Decg(i);
            x = i;
            flag = 1;
        end
        if Decg(i) <= 0.2*ymax && flag == 1;
            serieRRx(c)= i + 1 ;
            serieRRy(c)= ecgy(i + 1);
%             serieRRx(c)= i;
%             serieRRy(c)= ecgy(i);
            y= ymax * 0.5;
            flag= 0;
            c = c + 1;
        end
    end
    serieRRx = serieRRx';
    serieRRy = serieRRy';
end

function [serieRRx,serieRRy] = detectRPeaks(ecgy, Vecg)
  %  Detection of R wave using the instant varince
  %
  %  ecgy:        ECG signal
  %
  %  serieRRx:    X axis (time) position of R waves
  %  serieRRy:    Y axis (amplitude) of R waves
  %
  %  Reference:
  %
  %  implemented by: Francisco Perdigon Romero
  %  email: fperdigon88@gmail.com

    VecgMax = max(max(Vecg));
    [peaks_y,peaks_x] = findpeaks(Vecg, 'MINPEAKHEIGHT', VecgMax * 0.6, 'MINPEAKDISTANCE', 10);
    l = length(peaks_x);
    l1 = peaks_x(2) - peaks_x(1);
    l2 = peaks_x(3) - peaks_x(2);
    c = 1;
    if (l2 < l1)
        for i= 2:2:l-2
            tmp = Vecg(peaks_x(i):peaks_x(i+1));
            [p_y, p_x] = findpeaks(-tmp);
            serieRRx(c) = peaks_x(i) + p_x(1) - 1; %Hay que restarle una muestra
            serieRRy(c) = ecgy(serieRRx(c));
            c = c + 1;
        end
    else
        for i= 1:2:l-2
            tmp = Vecg(peaks_x(i):peaks_x(i+1));
            [p_y, p_x] = findpeaks(-tmp);
            serieRRx(c) = peaks_x(i) + p_x;
            serieRRy(c) = ecgy(serieRRx(c));
            c = c + 1;
        end
    end
end

function [serieRRx,serieRRy] = detectRPeaksTh(ecgy)
  %  Detection of R wave using Threshold (0.8 of the max)
  %
  %  ecgy:        ECG signal
  %
  %  serieRRx:    X axis (time) position of R waves
  %  serieRRy:    Y axis (amplitude) of R waves
  %
  %  Reference:
  %
  %  implemented by: Francisco Perdigon Romero
  %  email: fperdigon88@gmail.com
    ecgMax = max(max(ecgy));
    [serieRRy,serieRRx] = findpeaks(ecgy, 'MINPEAKHEIGHT',ecgMax * 0.8, 'MINPEAKDISTANCE', 10);
end

function [seriePQx,seriePQy] = detectPQ(ecgy, serieRRx, serieRRy, Vecg)
%  Detection of PQ interval using the instant variance and
%  the R points as reference
%
%  ecgy:        ECG signal
%  serieRRx:    X axis (time) position of R waves
%  serieRRy:    Y axis (amplitude) of R waves
%  Vecg:        Instant variance of the ECG signal TODO: calculate inside the method
%
%  seriePQx:    X axis (time) position of the center of PQ segment
%  seriePQy:    Y axis (amplitude) of the center PQ segment
%
%  Reference:
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    l = length(serieRRx);
    seriePQx = zeros(l,1);
    seriePQy = zeros(l,1);
    meanDrr = mean(diff(serieRRx));

    c= 1;
    % NMu  -> samples number before the R peak
    NMu = 10;
    for i=1:length(serieRRx)%
        Vtmp = Vecg(serieRRx(i)-floor(meanDrr/4) : (serieRRx(i) - NMu));
        [peaks_y,peaks_x] = findpeaks(-Vtmp);
        t = length(peaks_x);
        seriePQx(c) = (serieRRx(i) - NMu) - (length(Vtmp) - peaks_x(t));
        seriePQy(c) = ecgy(seriePQx(c));
        c = c + 1;
    end

end

function Vecg = instVariance(ecgy, m, scale)
%  Calculate the instant variance of the signal
%  in a windows size = m * 2 + 1

    l = length(ecgy);
    Vecg = zeros(l,1);
    [M,N]=size(ecgy);

    for i = m+1 : M - m
        Vecg(i) = var(ecgy(i-m : i+m));
    end

    Vecg = Vecg * scale;

end

function ecgy_BL = addSineBaseLine(ecgy, F, Fs, Amp)
%  Add artificial (sine) BLW to the ecg signal
%
%  ecgy:        ECG signal
%  F:           Frequency of the BLW
%  Fs:          Sampling Frequency
%  Amp:         The BLW  amplitude in proportion to the ECG (values between 0-1)
%
%  ecgy_BL:     ECG signal plus the artificial BLW
%
%  Reference:
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com


    [M,N] = size(ecgy);
    maxx = max([M,N]);
    ecgx = 0:maxx;
    ecgx = ecgx/Fs;
    AA = (max(ecgy)-min(ecgy)) * Amp;
    bl = AA * sin(2 * pi * F * ecgx);
    ecgy_BL = ecgy + bl;
end

function ecgy_LN = addLineNoise(ecgy, F, Fs)
  %  Add artificial electric line to the ecg signal
  %
  %  ecgy:        ECG signal
  %  F:           Frequency of the line noise (50 or 60)
  %  Fs:          Sampling Frequency
  %
  %  ecgy_LN:     ECG signal contaminated with artificial
  %               electric line noise
  %
  %  Reference:
  %
  %  implemented by: Francisco Perdigon Romero
  %  email: fperdigon88@gmail.com

    [M,N] = size(ecgy);
    maxx = max([M,N]);
    ecgx = zeros(maxx);
    AA = max(ecgy) * 0.005;
    ln = AA * sin(2 * pi * F/Fs * ecgx);
    ecgy_LN = ecgy + ln;
end

function ecgy_filt = lowPassFilter(ecgy, F, Fs)
%  Low pass IIR filter
%
%  ecgy:        ECG signal
%  F:           Cut-off Frequency for the filter
%  Fs:          Sampling Frequency
%
%  ecgy_filt:     ECG signal contaminated with artificial
%               electric line noise
%
%  Reference:
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    [b,a] = butter(4,F/Fs,'low');
    ecgy_filt = filtfilt (b,a,ecgy); % bidirectional filteroing
end

function ecgy_filt = highPassFilter(ecgy, F, Fs)
  %  High pass IIR filter
  %
  %  ecgy:        ECG signal
  %  F:           Cut-off Frequency for the filter
  %  Fs:          Sampling Frequency
  %
  %  ecgy_filt:     ECG signal contaminated with artificial
  %               electric line noise
  %
  %  Reference:
  %
  %  implemented by: Francisco Perdigon Romero
  %  email: fperdigon88@gmail.com

    [b,a] = butter(4,F/Fs,'high');
    ecgy_filt = filtfilt (b,a,ecgy); % bidirectional filteroing
end

function [signal t]= readECG_Discardia(filename)
%  Credits to Dicardia team
%  [signal t] = ReadECG(filename)
%
%  This function reads an ECG record from the DICARDIA database. It opens
%  the desired file and stores the 8 ECG derivations in the independent
%  variable 'signal'. The function also returns a time vector 't' that is
%  used to plot the signals.
%
%  Input:
%
%  'filename' is a sring containing the full path to the ECG record that
%   will be read
%
%  Output:
%
%  'signal' is an array of doubles containing 8 rows that correspond to
%  the ECG derivations contained in the file. The rows contain the DI,
%  DII, V1, V2, V3, V4, V5 and V6 derivations respectively. These signals
%  have a 500Hz sampling frequency
%
%  't' is a time axis that is coherent with the signals. It may be used to
%  plot the signals having a time reference.
%
%  This routine was developed by the Simon Bolivar University's Applied
%  Biophysics and Bioengineering Group (GBBA-USB) for free and public use.
%  It is protected under the terms of the Creative Commons Attribution-Non
%  Commercial 4.0 International License, no profit may come from any
%  application that uses this function and proper credit must be given to
%  the GBBA-USB. For further information consult
%  http://creativecommons.org/licenses/by-nc/4.0/

    freq = 500;
    ConversionResolution = 3.06e-3;
    NumberBits = 12;

    ders= 1:8;
    fid = fopen(filename,'r');
    A = fread(fid,'short');

    for i=1:8
        x = [A(i:8:end)];
        signal(i,:) = x(5000:end);
    end

    t = 0:1/freq:(length(signal)-1)/freq;
end

function [Y,f,h] = FFTplot(signal, Fs)
    T = 1/Fs;                     % Sample time
    L =  max(size(signal));                 % Length of signal


    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
    Y = fft(signal,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);

    % Plot single-sided amplitude spectrum.
    h = figure();
    set(gca, 'fontsize', 20, 'FontName','arial')

    plot(f,2*abs(Y(1:NFFT/2+1))) ;
    title('Single-Sided Amplitude Spectrum')
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    %xlim([0 5])
end

function [ecgx] = ECGx(ecgy, Fs)
    len = max(size(ecgy));
    ecgx = [0:1/Fs:len]
end
