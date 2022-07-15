%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This source file contains several Baseline Wonder removal
%  methods:
%
%  - Based on Cubic SPlines
%  - Based on FIR filters
%  - Based on IRR filters
%  - Based on LMS adaptative filters
%  - Based on moving-average filter
%  - Based on Independent Components Analysis (ICA)
%  - Based on Interpolation and Successive Subtraction of Median Values (ISSM)
%  - Based on Empirical Mode Decomposition (EMD)
%  - Based on Wavelet Transform
%
%  All these methods were programmed according the literature information.
%  The reference to these papers appear in the header information of each method
%
%  Author: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com
%  year: 2018
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BLremover = BLremover
    BLremover.SPLRemoveBL = @SPLRemoveBL;
    BLremover.FIRRemoveBL = @FIRRemoveBL;
    BLremover.IIRRemoveBL = @IIRRemoveBL;
    BLremover.IIR_b_RemoveBL = @IIR_b_RemoveBL;
    BLremover.FARemoveBL = @FARemoveBL;
    BLremover.MAFRemoveBL = @MAFRemoveBL;
    BLremover.ICARemoveBL = @ICARemoveBL;
    BLremover.ICARemoveBL_v2 = @ICARemoveBL_v2;
    BLremover.ISSMRemoveBL = @ISSMRemoveBL;
    BLremover.EMDRemoveBL = @EMDRemoveBL;
    BLremover.TWRemoveBL = @TWRemoveBL;
end

function [ECG_Clean] = TWRemoveBL(ecgy, Fs, Fc)
%  BLW removal method based on Wavelets Transform filter
%
%  ecgy:        the contamined signal
%  Fs:          sample frequiency
%  Fc:          cut-off frequency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Mozaffary B, Tinati MA. ECG Baseline Wander Elimination using Wavelet Packets.
%  World Acad Sci Eng Technol [Internet]. 2005
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    % Wavelet type
    w='sym10';
    % Threshold type
    thr_met='s';

    % Perform decomposition at level 5
    %[c,l] = wavedec(s,5,w);

    % The decomposition leve is estimated based in the cut-off frequency
    lev = ceil(log2(Fs/Fc));

    BL = wden(ecgy,'heursure',thr_met,'one',lev, w);

    ECG_Clean = ecgy - BL;

end

function [ECG_Clean] = ICARemoveBL(ecgy)
%  BLW removal method based on Independent Componet Analysis
%
%  ecgy:        the contamined signal
%  ECG_Clean :  processed signal without BLW
%  
%  FastICA source code and install:
%  http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml
%
%  Reference:
%  Barati Z, Ayatollahi A. Baseline Wandering Removal by Using Independent Component
%  Analysis to Single-Channel ECG data. 2006 Int Conf Biomed Pharm Eng. 2006;(1):152–6.
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    % The signal is shifted (11 to 20 sampled in random) to create a
    % multichanel signal of 60 channels.

    addpath('FastICA_25/')
    t_ecg = max(size(ecgy));
    ECGMix = zeros(60, t_ecg + 20);
    ECGMix(1,1:t_ecg) = ecgy';
    for i = 2:60
        delay = randi([11,20]);
        ECGMix(i,delay:delay + t_ecg-1) = ecgy';
    end

    %preform ica to unmix signal
    [ICA, A, W] = fastica(ECGMix);

    LICA = min(size(ICA));
    % BL have negative kurtosis, eliminate them
    for i=1:LICA
        k(i) = kurtosis(ICA(i, :)) - 3;

        if k(i) < 0
            ICA(i,:) = zeros(1, t_ecg + 20);
        end
    end

    ECGCleanMix =  (ICA' * W)';
    ECG_Clean = (ECGCleanMix(1,1:t_ecg))';

end

function [ECG_Clean] = ICARemoveBL_v2(ecgy)
%  BLW removal method based on Independent Componet Analysis
%
%  ecgy:        the contamined signal
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Barati Z, Ayatollahi A. Baseline Wandering Removal by Using Independent Component
%  Analysis to Single-Channel ECG data. 2006 Int Conf Biomed Pharm Eng. 2006;(1):152–6.
%  
%  FastICA source code and install:
%  http://research.ics.aalto.fi/ica/fastica/code/dlcode.shtml
%
%  This version of the implementation implementa a bias automatic removal (do not affect the signal shape)
%  This version is the used for the comparative estudy
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com


    % The signal is shifted (11 to 20 sampled in random) to create a
    % multichanel signal of 60 channels.

    addpath('FastICA_25/')
    t_ecg = max(size(ecgy));
    ECGMix = zeros(60, t_ecg + 20);
    ECGMix(1,1:t_ecg) = ecgy';
    for i = 2:60
        delay = randi([11,20]);
        ECGMix(i,delay:delay + t_ecg-1) = ecgy';
    end

    %preform ica to unmix signal
    [ICA, A, W] = fastica(ECGMix);

    LICA = min(size(ICA));
    % BL have negative kurtosis, eliminate them
    for i=1:LICA
        k(i) = kurtosis(ICA(i, :)) - 3;

        if k(i) < 0
            ICA(i,:) = zeros(1, t_ecg + 20);
        end
    end

    ECGCleanMix =  (ICA' * W)';
    ECG_Clean = (ECGCleanMix(1,1:t_ecg))';
    %%% Out adjust
    mean_in = mean(ecgy);
    mean_out =  mean(ECG_Clean);
    FacAmp = mean_out/mean_in; % To remove a big biases
    bias =  0.0;
    ECG_Clean = (ECG_Clean/FacAmp) - bias;
end

function [ECG_Clean] = IIRRemoveBL(ecgy,Fs)
%  BLW removal method based on IIR Filter
%
%  ecgy:        the contamined signal
%  Fs:          sample frequiency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Pottala EW, Bailey JJ, Horton MR, Gradwohl JR. Suppression of baseline wander in the ECG
%  Using a bilinearly transformed, null-phase filter. J Electrocardiol [Internet].
%  1990 Jan [cited 2016 May 6];22:243–7.
%  This the same implementation an coficients used by the autors, the results were not good
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    % Damping coefficient
    delta = 0.6;
    % Cut frequency
    fc =0.65;
    g = tan(pi*fc*(1/Fs)*0.001);
    g1 = 1 + (2*delta*g) + (g^2);
    g2 = (2*g^2) - 2;
    g3 = 1 - (2*delta*g) + (g^2);
    b1 = g^2/g1;
    b2 = 2*g^2/g1;
    b3 = b1;
    a2 = g2/g1;
    a3 = g3/g1;

    % Forward Pass
    % Initials conditions
    ECG_Clean(1) = (b1 + b2 + b3 - a2 -a3) * ecgy(1);
    ECG_Clean(2) = b1 * ecgy(2) +(b2 + b3 -a3) * ecgy(1) - a2 * ECG_Clean(1);
    for i = 3: max(size(ecgy))
       ECG_Clean(i) = b1 * ecgy(i) + b2 * ecgy(i-1) + b3 * ecgy(i-2) - a2 * ECG_Clean(i-1) - a3 * ECG_Clean(i-2);
    end

    % Backward Pass
    ECG_Clean = fliplr(ECG_Clean);
    % Initials conditions
    ECG_Clean(1) = (b1 + b2 + b3 - a2 -a3) * ecgy(1);
    ECG_Clean(2) = b1 * ecgy(2) +(b2 + b3 -a3) * ecgy(1) - a2 * ECG_Clean(1);
    for i = 3: max(size(ecgy))
       ECG_Clean(i) = b1 * ecgy(i) + b2 * ecgy(i-1) + b3 * ecgy(i-2) - a2 * ECG_Clean(i-1) - a3 * ECG_Clean(i-2);
    end

    ECG_Clean = ecgy - ECG_Clean';

end

function [ECG_Clean] = IIR_b_RemoveBL(ecgy,Fs,Fc)
%  BLW removal method based on IIR Filter
%
%  ecgy:        the contamined signal
%  Fc:          cut-off frequency
%  Fs:          sample frequiency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  https://www.mathworks.com/help/signal/ref/butter.html
%
%  This variant use the Matlab butterword calculation coficients for a IIR filter
%  the results were better that the original implementation of the paper
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    [b,a] = butter(4,Fc/(Fs/2)),'high');

    ECG_Clean = filtfilt (b,a,ecgy); % filtrado bidirecional
end

function [ECG_Clean] = FIRRemoveBL(ecgy,Fs,Fc)
%  BLW removal method based on FIR Filter
%
%  ecgy:        the contamined signal
%  Fc:          cut-off frequency
%  Fs:          sample frequiency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Van Alsté JA, Schilder TS. Removal of base-line wander and power-line interference from the ECG
%  by an efficient FIR filter with a reduced number of taps. IEEE Trans Biomed Eng [Internet].
%  1985 Dec [cited 2016 May 6];32(12):1052–60.
%
%  implemented by: Francisco Perdigon Romero
%  modify by: David Castro Pi�ol 17/03/2020
%       modification regarding the expansion of the signal to ensure the
%       signal's length be greater than three times the filter order before
%       using filtfilt function
%
%  email: fperdigon88@gmail.com

    fcuts = [(Fc-0.07) (Fc)];
    mags = [0 1];
    devs = [0.005 0.001];
    len = length(ecgy);

    [n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,Fs);
    
    b = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
    
    if(n*3>length(ecgy))
        diff = n*3 - length(ecgy);
        ecgy = [ecgy;flip(ecgy);ecgy(1)*ones(diff,1)];
    end
       
    a = 1;
    ECG_Clean = filtfilt(b,a,ecgy);
    ECG_Clean = ECG_Clean(1:len,1);
end

function [ECG_Clean] = EMDRemoveBL(ecgy,Fs, Fc)
%  BLW removal method based on Empirical Mode Decomposition
%  EMD toolbox developed by Rilling and Flandrin.
%
%  EMD library source code and install:
%  http://perso.ens-lyon.fr/patrick.flandrin/emd.html
%
%  ecgy:        the contamined signal
%  Fc:          cut-off frequency
%  Fs:          sample frequiency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Blanco-Velasco M, Weng B, Barner KE. ECG signal denoising and baseline wander
%  correction based on the empirical mode decomposition. Comput Biol Med. 2008;38(1):1–13.
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    % Import utilECG resources
    uE = utilECG;

    addpath('EMD_Baseline/package_emd/')

    % Always install the package just in case ;)
    install_emd
    EMD_ecg = emd(ecgy);

    % Visualize IMF
    % emd_visu(ecgy,EMD_ecg)

    [N_imf, L] = size(EMD_ecg);

    for i = N_imf: -1: N_imf - 8 % N_imf-9
        EMD_ecg(i,:) = uE.highPassFilter(EMD_ecg(i,:), Fc, Fs);
    end

    Final_ecg = zeros(L,1);
    for i = 1 : N_imf
        Final_ecg = Final_ecg + EMD_ecg(i,:)';
    end

    ECG_Clean = Final_ecg;
end

function [ECG_Clean] = SPLRemoveBL(ecgy,Fs)
%  BLW removal method based on cubic SPlines
%
%  ecgy:        the contamined signal
%  Fs:          sample frequiency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Meyer CR, Keiser HN. Electrocardiogram baseline noise estimation and removal using cubic splines
%  and state-space computation techniques. Comput Biomed Res [Internet]. 1977
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com


    % Import utilECG resources
    uE = utilECG;

    [serieRRx,serieRRy] = uE.detectRPeaks_FDeriv(ecgy);

    len = max(size(ecgy));
    ecgx = 0:len-1;
    %Obtencion de las series PQ, segun el paper 66ms antes del R
    % 66 ms en muestras
    dif_x = ceil((66/1000)*Fs);
    seriePQx = serieRRx - dif_x;
    lenPQx = max(size(seriePQx));
    for i = 2 : lenPQx
        seriePQy(i) = ecgy(seriePQx(i));
    end
    % Add the begining and end
    seriePQx = [1; seriePQx; len];
    seriePQy = [ecgy(1); seriePQy; ecgy(len)];

    ebl_spline = spline(seriePQx, seriePQy, ecgx)';

    ECG_Clean = ecgy - ebl_spline;
end

function [ECG_Clean] = MAFRemoveBL(ecgy,Fs, Fc)
%  BLW removal method based on moving-average filter
%
%  ecgy:        the contamined signal
%  Fc:          cut-off frequency
%  Fs:          sample frequiency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Canan S, Ozbay Y, Karlik B. A method for removing low varying frequency trend
%  from ECG signal. In: Biomedical Engineering Days, 1998 Proceedings of the 1998
%  2nd International Conference. IEEE; 1997. p. 144–5.
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    N = max(size(ecgy));
    m = ceil(Fs/(2*Fc));
    if (mod(m,2) == 0)
        m = m + 1;
    end
    r = (m - 1) / 2;
    y = zeros(1,N);

    for n=(r+1):(N-r)
        sum = 0;
        for i=n-r:n+r
            sum = sum + ecgy(i);
        end
        % stimated BLW
        y(n) =  sum/m;
    end

    ECG_Clean = ecgy - y';
end

function [ECG_Clean] = ISSMRemoveBL(ecgy,Fs,Fc)
%  BLW removal method based on interpolation and successive
%  subtraction of median values in RR intervals
%
%  ecgy:        the contamined signal
%  Fc:          cut-off frequency
%  Fs:          sample frequiency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Chouhan VS, Mehta SSS. Total removal of baseline drift from ECG signal.
%  In: Computing: Theory and Applications, 2007 ICCTA’07 International Conference
%  on [Internet]. IEEE; 2007
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com

    % Import utilECG resources
    uE = utilECG;
    [serieRRx,serieRRy] = uE.detectRPeaks_FDeriv(ecgy);
    
    N = max(size(ecgy));
    
    y = zeros(1,N);
    
    % On the original paper authors dont stablish a number if iterations or a stop criteria,
    % we set the number of iterations = 100, then stop the algorithm
    iterations = 100;
    for iter=1:iterations  
        % Fill y values from elemt 1 up to first R peak
        tempRR_size = serieRRx(n) - 1 
        y(1 : serieRRx(n + 1)) = ones(1, tempRR_size) * median(ecgy(1 : serieRRx(n));
        
        % Filly values for the rest of the signal lengt using median values between RR intervals
        for n=1:(size(serieRRx) -1)
            tempRR_size = serieRRx(n + 1) - serieRRx(n)
            y(serieRRx(n) : serieRRx(n + 1)) = ones(1,tempRR_size) * median(ecgy(serieRRx(n) : serieRRx(n + 1));
        end
        % Substract the estimated BLW to the signal
        ECG_Clean = ecgy - y'
    end
end

function [ECG_Clean] = FARemoveBL(ecgy,Fs)
%  BLW removal method based on LMS adaptavive filters
%
%  ecgy:        the contamined signal
%  Fc:          cut-off frequency
%  Fs:          sample frequiency
%  ECG_Clean :  processed signal without BLW
%
%  Reference:
%  Laguna P, Jane R, Caminal P. Adaptive filtering of ECG baseline wander. In: Engineering in Medicine and
%  Biology Society, 1992 14th Annual International Conference of the IEEE. 1992. p. 508–9.
%
%  implemented by: Francisco Perdigon Romero
%  email: fperdigon88@gmail.com


    l = max(size(ecgy));
    x1 = ones(l,1);
    
    mu = 0.001;            % LMS step size.
    ha = adaptfilt.lms(32,mu);
    [y1,e1] = filter(ha,x1,ecgy);

    uE = utilECG;
    [serieRRx,serieRRy] = uE.detectRPeaks_FDeriv(e1);

    l2 = max(size(serieRRx));

    x2 = zeros(l,1);

    for i=1:l2
        x2(serieRRx(i)) = 1;
    end

    mu = 0.05;            % LMS step size.
    DserieRRx = diff(serieRRx);
    meanRR = mean(DserieRRx);
    ha = adaptfilt.lms(ceil(meanRR),mu);
    [y2,e2] = filter(ha,x2,e1);


    %figure()
    %plot(time,y2)
    ECG_Clean = y2;

end
