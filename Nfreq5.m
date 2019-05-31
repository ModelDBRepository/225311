function [freq_signal, f] = Nfreq5(v,tv)
% Finds frequecy spectrum of time series v (at time pts tv)

v = v-mean(v);

% Interpolation
difftv = diff(tv); % In some cases the same data point occurs twice. Remove duplicates:
keepind = find(difftv);
tv = tv(keepind);
v = v(keepind);

NFFT = 2^(nextpow2(length(v)));    % Next power of 2 from length of y. This optimizes the fft.

tin = linspace(min(tv), max(tv), NFFT); % interpolate the signal so it has this length
v = interp1(tv, v, tin);
dt = tin(5)-tin(4); % All time points were equally spaced

% Sample
samp_freq = 1/dt;
freq_signal = fft(v);
freq_signal = fft(v,NFFT)/length(v);


% TO ONLY LOOK AT FIRST HALF OF FREQ. SPEC.
f = samp_freq/2*linspace(0,1,NFFT/2+1);
freq_signal = abs(freq_signal(1:NFFT/2+1));
