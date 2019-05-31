function PLoSfreqplot(comp, S1, S2)
% Compares two datafiles in temporal development
% S1: Diffusion on
% S2: Diffusion off

Ncomp = 1;
flatten = 1; % Reduce number of data points (1 per 0.1 log unit)
dopower = 1; % 1 for power spectum (0 for amplitude spetrum)
howmany = 4; % split signal into this number of parts

ton = S1.Simdata.t(2:end-1);
Von = S1.Simdata.V(2:end-1,3);
toff = S2.Simdata.t(2:end-1);
Voff = S2.Simdata.V(2:end-1,3);

Von = Von*1000; % convert to mV
Voff = Voff*1000;


% FIG 1: Compare case No diff between different intervals
Figcompf = figure;
set(Figcompf, 'Color', [1 1 1]);
myaxis = [-1.5 3.5 -10 -4];

for i = 1:howmany
    AAoff = 1 + floor((i-1)*(length(Voff)/howmany));
    BBoff = floor(i*(length(Voff)/howmany));
    AAon = 1 + floor((i-1)*(length(Von)/howmany));
    BBon = floor(i*(length(Von)/howmany));
    if BBoff>length(toff)
        BBoff = length(toff);
    end
    if BBon>length(ton)
        BBon = length(ton);
    end    
    offint = AAoff:BBoff;
    onint = AAon:BBon;
    myton = ton(onint);
    myVon = Von(onint,:);
    mytoff = toff(offint);
    myVoff = Voff(offint,:);

   [poff, foff] = Nfreq5(myVoff,mytoff);
   [pon, fon] = Nfreq5(myVon,myton);
   
   if flatten
       [poff, foff] = Fflatten(poff, foff); 
       [pon, fon] = Fflatten(pon, fon);
   end
   
   if dopower
       poff = poff.^2;
       pon = pon.^2;
   end       
   figure(Figcompf); subplot(1,howmany,i); hold on;
   plot(log10((foff(2:end))),log10((poff(2:end))),'color', 'r');
   plot(log10((fon(2:end))),log10((pon(2:end))),'color', 'b');
   axis(myaxis);
 %  title(['Interval ', num2str(i)])
   offtestsum = sum(poff)
   ontestsum = sum(pon)
end

figure(Figcompf);
tpts = [0 21 42 63 84];
letters = 'ABCD';
for i = 1:howmany
subplot(1,howmany,i)
xlabel('log_{10}(f)')
ylabel('log_{10}(PSD)')
titlestring = [letters(i), ') ', 't = ', num2str(tpts(i)), ' - ', num2str(tpts(i+1)), ' s'];
title(titlestring);
end
subplot(1,howmany,howmany)
legend('diff off', 'diff on', 'neurons off');



function [pp, ff] = Fflatten(pp, ff)
% Filter signal
% Interpolate to get one data point per 0.1 log unit
lff = log10(ff);
lff(1) = [];
lminf = min(lff);
lmaxf = max(lff);
lnewff = lminf:0.1:lmaxf;
newff = 10.^lnewff;
newpp = interp1(ff, pp, newff);
pp = newpp;
ff = newff;
