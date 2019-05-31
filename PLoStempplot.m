function PLoStempplot(Sdiff, Snodiff)
% Compares two datafiles in temporal development

c1a = 'r';
c1b = 'r--';
c2a = 'b';
c2b = 'b--';

h = figure; set(h, 'Tag', 'temporalplot', 'color', [1 1 1]);

cplotit(Snodiff, c1a, c1b)
cplotit(Sdiff, c2a, c2b)

vstep = 1/4+0.028;
hstep = 0.47;

gca;
legend('soma nodiff', 'apical nodiff', 'soma diff', 'apical diff');


% hh1 = uicontrol('Style','text', 'String', 'A');
% set(hh1, 'Units', 'Normalized', 'BackgroundColor', [1 1 1]);
% set(hh1, 'Position', [0.035 0.89 0.03 0.05], 'Fontsize', 14);
% 
% hh2 = uicontrol('Style','text', 'String', 'B');
% set(hh2, 'Units', 'Normalized', 'BackgroundColor', [1 1 1]);
% set(hh2, 'Position', [0.035 0.89-vstep*1 0.03 0.05], 'Fontsize', 14);
% 
% hh3 = uicontrol('Style','text', 'String', 'C');
% set(hh3, 'Units', 'Normalized', 'BackgroundColor', [1 1 1]);
% set(hh3, 'Position', [0.035 0.89-vstep*2 0.03 0.05], 'Fontsize', 14);


function cplotit(S, ca, cb)
t = S.Simdata.t'; %(2:end-1); % just to eliminate some endpoint bugs
tvr = S.Simdata.tvr';
cV = S.Simdata.V';% (2:end-1,:);
cVr = S.Simdata.Vr';
cK = S.Simdata.cK'; 
cNa = S.Simdata.cNa'; 
cCa = S.Simdata.cCa'; 
cX = S.Simdata.cX'; 
Avox = S.geometry.Avox;
diffon = S.diffon;
lambda_0 = S.diffconsts(1);
D_K = S.diffconsts(2);
D_Na = S.diffconsts(3);
D_Ca = S.diffconsts(4);
D_X = S.diffconsts(5);
deltax = S.geometry.deltax;

% Some constants:
F= 96485.3365; % C/mol
T = 300; % K
R = 8.3; % J/mol/K
Npsi = R*T/F; % V


% Extracellular fluxes: (mol/s)
jKEd = - diffon*Avox*( D_K*(cK(2:end,:) - cK(1:end-1,:))/deltax);
jKEf = - Avox*D_K/Npsi*(cK(2:end,:) + cK(1:end-1,:))/2.*(cV(2:end,:) - cV(1:end-1,:))/deltax; % N+1 elements
jNaEd =  - diffon*Avox*( D_Na*(cNa(2:end,:) - cNa(1:end-1,:))/deltax);
jNaEf = - Avox*D_Na/Npsi*(cNa(2:end,:) + cNa(1:end-1,:))/2.*(cV(2:end,:) - cV(1:end-1,:))/deltax; % N+1 elements
jCaEd = - diffon*Avox*( D_Ca*(cCa(2:end,:) - cCa(1:end-1,:))/deltax);
jCaEf = - Avox*2*D_Ca/Npsi*(cCa(2:end,:) + cCa(1:end-1,:))/2.*(cV(2:end,:) - cV(1:end-1,:))/deltax; % N+1 elements
jXEd =  - diffon*Avox*( D_X*(cX(2:end,:) - cX(1:end-1,:))/deltax);
jXEf = + Avox*D_X/Npsi*(cX(2:end,:) + cX(1:end-1,:))/2.*(cV(2:end,:) - cV(1:end-1,:))/deltax; % N+1 elements

% Extracellular currents: (A)
iEd = F*(2*jCaEd + jKEd + jNaEd - jXEd);
iEf = F*(2*jCaEf + jKEf + jNaEf - jXEf);
%iEd = iEd/F; iEf = iEf/F;

iEdin = iEd(1,:);
iEdout = iEd(end,:);
iEfin = iEf(1,:);
iEfout = iEf(end,:);
netdin = iEdin-iEdout;
netfin = iEfin-iEfout;
netin = netdin + netfin;

cK = cK'; cNa = cNa'; cCa = cCa'; cX = cX'; 

comps = 1:size(cNa,1); % number of compartments
Nsoma = 3; % selected compartment
Ntrunk = 7;
Napical = 13; % selected compartment

tms = t*1000; % ms
briefint = [2000 2040];
briefind1 = find(tms>briefint(1) & tms<briefint(2)); % to plot only 20 ms
briefint = [30000 30040];
briefind2 = find(tms>briefint(1) & tms<briefint(2)); % to plot only 20 ms
briefint = [83000 83040];
briefind3 = find(tms>briefint(1) & tms<briefint(2)); % to plot only 20 ms
tms = t; 

ICfac = 1/Avox; % Convert A to A/m^2.
JCfac = 1/Avox*1e6; % Convert mol/s to micromol/m^2/s




%%%%%%%% TIME DEVELOPMENT OF V and I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = findobj('Tag', 'plotteporal');
set(h, 'Color', [1 1 1], 'Position', [94 311 818 520]);

subplot(3,4,1); hold on; 
plot(tms, 1000*cV(Nsoma, :), ca);
plot(tms, 1000*cV(Napical, :), cb); 
title('A1', 'FontSize', 14)
ylabel('V (mV)', 'FontSize', 12)
axis([0 85 -3 2]);


subplot(3,4,2); hold on; 
plot(tms(briefind1), 1000*cV(Nsoma, briefind1), ca); 
plot(tms(briefind1), 1000*cV(Napical, briefind1), cb); 
title('A2', 'FontSize', 14)
%axis([2 2.04 -2 0.5]);

subplot(3,4,3); hold on; 
plot(tms(briefind2), 1000*cV(Nsoma, briefind2), ca); 
plot(tms(briefind2), 1000*cV(Napical, briefind2), cb); 
title('A3', 'FontSize', 14)
%axis([35 35.04 -2 0.5]);

subplot(3,4,4); hold on; 
plot(tms(briefind3), 1000*cV(Nsoma, briefind3), ca); 
plot(tms(briefind3), 1000*cV(Napical, briefind3), cb); 
title('A4', 'FontSize', 14)
%axis([82 82.04 -2 0.5]);


%ylabel('V (mV)')

subplot(3,4,5); hold on;
plot(tms, ICfac*iEf(Nsoma,:), ca);
plot(tms, ICfac*iEf(Napical,:), cb);
ylabel('i^f (A/m^2)', 'FontSize', 12)
title('B1', 'FontSize', 14)
axis([0 85 -20 10]);

%axis([0 8000 -10 5]);

subplot(3,4,6); hold on;
plot(tms(briefind1), ICfac*iEf(Nsoma,briefind1), ca);
plot(tms(briefind1), ICfac*iEf(Napical,briefind1), cb);
title('B2', 'FontSize', 14)
%axis([2 2.04 -10 5]);

subplot(3,4,7); hold on;
plot(tms(briefind2), ICfac*iEf(Nsoma,briefind2), ca);
plot(tms(briefind2), ICfac*iEf(Napical,briefind2), cb);
title('B3', 'FontSize', 14)
%axis([35 35.04 -10 5]);

subplot(3,4,8); hold on;
plot(tms(briefind3), ICfac*iEf(Nsoma,briefind3), ca);
plot(tms(briefind3), ICfac*iEf(Napical,briefind3), cb);
title('B4', 'FontSize', 14)
%axis([82 82.04 -10 5]);


%axis([briefint -10 5]);

subplot(3,4,9); hold on;
plot(tms, ICfac*iEd(Nsoma,:), ca);
plot(tms, ICfac*iEd(Napical,:), cb);
ylabel('i^d (A/m^2)', 'FontSize', 12)
title('C1', 'FontSize', 14)
axis([0 85 -0.2 0.6]);
xlabel('t(s)', 'FontSize', 12)

subplot(3,4,10); hold on;
plot(tms(briefind1), ICfac*iEd(Nsoma,briefind1), ca);
plot(tms(briefind1), ICfac*iEd(Napical,briefind1), cb);
title('C2', 'FontSize', 14)
%axis([2 2.04 -0.05 0.6]);
xlabel('t(s)', 'FontSize', 12)


subplot(3,4,11); hold on;
plot(tms(briefind2), ICfac*iEd(Nsoma,briefind2), ca);
plot(tms(briefind2), ICfac*iEd(Napical,briefind2), cb);
title('C3', 'FontSize', 14)
%axis([35 35.04 -0.05 0.6]);
xlabel('t(s)', 'FontSize', 12)

subplot(3,4,12); hold on;
plot(tms(briefind3), ICfac*iEd(Nsoma,briefind3), ca);
plot(tms(briefind3), ICfac*iEd(Napical,briefind3), cb);
title('C4', 'FontSize', 14)
%axis([82 82.04 -0.05 0.6]);
xlabel('t(s)', 'FontSize', 12)
