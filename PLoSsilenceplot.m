function PLoSsilenceplot(S1, S2, S3)
% Plot results for case with neurons were turned off midways in simulation
% Compare with default cases
%
% S1: Diffusion on
% S2: Diffusion off
% S3: Neurons turned off midways in simulation

h = figure('Tag','simple');
subplot(2,9,1);
set(h,'Color', [1 1 1]);
%set(h, 'Position', [40 370 1200  430]);
plotcprof(S3);
plotfreqs(S1,S2,S3);

figure(h);
subplot(2,9,1); title('A1', 'Fontsize', 14);
subplot(2,9,[2 3]); title('A2', 'Fontsize', 14);
xlabel('[Na^+] (mM)', 'Fontsize', 12);

subplot(2,9,[4 5]); title('A3', 'Fontsize', 14);
xlabel('[K^+] (mM)', 'Fontsize', 12);

subplot(2,9,[6 7]); title('A4', 'Fontsize', 14);
xlabel('[Ca^{2+}] (mM)', 'Fontsize', 12);

subplot(2,9,[8 9]); title('A5', 'Fontsize', 14);
xlabel('[X^-] (mM)', 'Fontsize', 12);

subplot(2,9,10); title('B1', 'Fontsize', 14);
subplot(2,9,[12 13]); title('B2', 'Fontsize', 14);
xlabel('V (mV)', 'Fontsize', 12);

subplot(2,9,[15 16]); title('B3', 'Fontsize', 14);
subplot(2,9,[17 18]); title('B4', 'Fontsize', 14);


function plotcprof(S)
h = findobj('Tag', 'simple');
Dcol = 'k';

% Useful parameters
C_m = 1.00e-2; % Membrane capacitance (Farad/m^2);
F = 96485.3365; % C/mol
T = 300; % K
R = 8.3; % J/mol/K
psi = R*T/F; % V


N = S.Nvox;
deltax = S.geometry.deltax;
Avox = S.geometry.Avox;

t = S.Simdata.t;
keepinds = find(t>ceil(max(t)/2));
t = t(keepinds);
tv = t;

cNa = S.Simdata.cNa(keepinds,:);
cK = S.Simdata.cK(keepinds,:);
cCa = S.Simdata.cCa(keepinds,:);
cX = S.Simdata.cX(keepinds,:);
V = S.Simdata.V(keepinds,:);

diffconsts = S.diffconsts; 
lambda_o = diffconsts(1);
D_K = diffconsts(2);
D_Na = diffconsts(3);
D_Ca = diffconsts(4);
D_X = diffconsts(5);
diffon = S.diffon; 

tv = tv(2:end-1); % just to eliminate some endpoint bugs
V = V(2:end-1,:);
V = [V(1,:); V(:,:); V(end,:)]';
Npsi = 0.0258;

comps = 1:size(cNa,2); % number of compartments
timepts = round(linspace(2,length(t),6)); % Selected time pts.
zeroline = zeros(1,15);

%%%%%%%%% PLOT SPATIAL PROFILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Spatial profiles %%%%%%
h = findobj('Tag','simple');

subplot(2,9,1); 
hold on;
% make a drawing of the neuron
plot(3,3,'^', 'markersize', 8, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
plot([3 1.7], [3 1.7], '-', 'Linewidth', 1.5, 'Color', Dcol);
plot([3 4.3], [3 1.7], '-', 'Linewidth', 1.5, 'Color', Dcol)
plot([3 3], [3 13], '-', 'Linewidth', 2.5, 'Color', Dcol)
plot([3 1.5], [13 14], '-', 'Linewidth', 1.5, 'Color', Dcol)
plot([3 4.5], [13 14], '-', 'Linewidth', 1.5, 'Color', Dcol)
axis([0 6 1 15]);
%set(gca, 'Yticklabel', '')
set(gca, 'Xticklabel', '')

subplot(2,9,10); 
hold on;
% make a drawing of the neuron
plot(3,3,'^', 'markersize', 8, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
plot([3 1.7], [3 1.7], '-', 'Linewidth', 1.5, 'Color', Dcol);
plot([3 4.3], [3 1.7], '-', 'Linewidth', 1.5, 'Color', Dcol)
plot([3 3], [3 13], '-', 'Linewidth', 2.5, 'Color', Dcol)
plot([3 1.5], [13 14], '-', 'Linewidth', 1.5, 'Color', Dcol)
plot([3 4.5], [13 14], '-', 'Linewidth', 1.5, 'Color', Dcol)
axis([0 6 1 15]);
%set(gca, 'Yticklabel', '')
set(gca, 'Xticklabel', '')


subplot(2,9,[2 3]);
plot(cNa(timepts,:)', comps);
set(gca, 'Yticklabel', '')
haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);

subplot(2,9,[4 5]); 
plot(cK(timepts,:)', comps);
haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);
set(gca, 'Yticklabel', '')

subplot(2,9,[6 7]); 
plot(cCa(timepts,:)', comps);
haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);
set(gca, 'Yticklabel', '')

subplot(2,9,[8 9]); 
plot(cX(timepts,:)', comps);
haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);
set(gca, 'Yticklabel', '')
legt = round(10*t(timepts))/10; mylegend = [num2str(legt), ['ssssss']'];
legend(mylegend);
set(gca, 'Yticklabel', '')


% Plot V-profiles
%avgV = zeros(15,length(timepts)-1);
%for ii = 2:length(timepts)
%avgV(:,ii-1) = mean(V(:,timepts(ii-1):timepts(ii)),2);
%end
avgV=V(:,timepts); %just plot at time points (no filtering)

subplot(2,9,[12 13]);
arnie = get(gca, 'ColorOrder');
%set(gca,'ColorOrder', arnie(2:end,:));
hold on;
plot(1000*avgV', comps);
%set(gca, 'Yticklabel', '')
haxis = axis; haxis(1) = -0.2; haxis(2) = 0.1; haxis(3) = 1; haxis(4) = 15; axis(haxis);

%legt2 = round(10*t(timepts(2:end)))/10; 
%legt1 = [42; legt2(1:end-1)];
%mylegend = [num2str(legt1), ['     ']', ['-----']', ['     ']', num2str(legt2), ['sssss']'];
legt = round(10*t(timepts))/10; mylegend = [num2str(legt), ['ssssss']'];
legend(mylegend);

legend(mylegend);
set(gca, 'Yticklabel', '')






function plotfreqs(S1,S2,S3)
Ncomp = 1;
comp = 3; %soma
flatten = 1;
dopower = 1;

ton = S1.Simdata.t;
ton(end) = []; %endpoint bug
Von = S1.Simdata.V(:,comp)*1000; %mV
Von(end) = [];
%clear S1;

toff = S2.Simdata.t;
toff(end) = [];
Voff = S2.Simdata.V(:,comp)*1000;
Voff(end) = [];
%clear S2;

tdead = S3.Simdata.t;
tdead(end) = [];
Vdead = S3.Simdata.V(:,comp)*1000;
Vdead(end) = [];
%clear S3;

AA1 = 1 + floor(1/2*(length(Voff)));
BB1 = floor(3/4*(length(Voff)));
AA2 = 1 + floor(3/4*(length(Voff)));
BB2 = floor(length(Voff))-1;

myint1 = AA1:BB1;
myint2 = AA2:BB2;

myVon1 = Von(myint1,:);
myVoff1 = Voff(myint1,:);
mytoff1 = toff(myint1);
myton1 = ton(myint1);
myVon2 = Von(myint2,:);
myVoff2 = Voff(myint2,:);
mytoff2 = toff(myint2);
myton2 = ton(myint2);
myVdead1 = Vdead(myint1,:);
mytdead1 = tdead(myint1);
myVdead2 = Vdead(myint2,:);
mytdead2 = tdead(myint2);


[poff1, foff1] = Nfreq5(myVoff1,mytoff1);
[pon1, fon1] = Nfreq5(myVon1,myton1);
[poff1, foff1] = Fflatten(poff1, foff1); 
[pon1, fon1] = Fflatten(pon1, fon1);

[poff2, foff2] = Nfreq5(myVoff2,mytoff2);
[pon2, fon2] = Nfreq5(myVon2,myton2);
[poff2, foff2] = Fflatten(poff2, foff2); 
[pon2, fon2] = Fflatten(pon2, fon2);

[pdead1, fdead1] = Nfreq5(myVdead1,mytdead1);
[pdead2, fdead2] = Nfreq5(myVdead2,mytdead2);
[pdead1, fdead1] = Fflatten(pdead1, fdead1); 
[pdead2, fdead2] = Fflatten(pdead2, fdead2); 

poff1 = poff1.^2; %convert from amplitude to power
pon1 = pon1.^2;
poff2 = poff2.^2;
pon2 = pon2.^2;
pdead1 = pdead1.^2;
pdead2 = pdead2.^2;

Figcompf = findobj('Tag', 'simple'); 
figure(Figcompf);
myaxis = [-1.5 3.5 -10 -4];

subplot(2,9,[15:16]);
hold on;
plot(log10((foff1(2:end))),log10((poff1(2:end))),'color', 'r');
plot(log10((fon1(2:end))),log10((pon1(2:end))),'color', 'b');
plot(log10((fdead1(2:end))),log10((pdead1(2:end))),'color', 'k');
xlabel('log(f)')
ylabel('log(power)')
legend('diff off', 'diff on', 'neurons off');
axis(myaxis);

subplot(2,9,[17:18]);
hold on;
plot(log10((foff2(2:end))),log10((poff2(2:end))),'color', 'r');
plot(log10((fon2(2:end))),log10((pon2(2:end))),'color', 'b');
plot(log10((fdead2(2:end))),log10((pdead2(2:end))),'color', 'k');
xlabel('log(f)')
%ylabel('log(power)')
legend('diff off', 'diff on', 'neurons off');
axis(myaxis);

%titlestring = [letters(i), ') ', 't = ', num2str(tpts(i)), ' - ', num2str(tpts(i+1)), ' s'];
%title(titlestring);


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