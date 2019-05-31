function PLoSdynplot(Sdiff,Snodiff)
%close all; 
c1 = 'r';
c2 = 'b';

h = figure;
set(h, 'Tag', 'vprofiles');
set(h,'Color', [1 1 1]);
set(h, 'Position', [40 370 1200  430]);

h2 = figure;
set(h2, 'Tag', 'cprofiles2');
set(h2,'Color', [1 1 1]);
set(h2, 'Position', [40 370 1200  430]);

hv = figure;
set(hv, 'Tag', 'vvprofiles');
set(hv,'Color', [1 1 1]);
set(hv, 'Position', [40 370 1200  430]);

hh = figure; 
set(hh, 'Tag', 'jprofiles');
set(hh,'Color', [1 1 1]);
set(hh, 'Position', [57    49   914   767]);

hhhh = uicontrol('Style','text');
set(hhhh, 'Units', 'Normalized', 'BackgroundColor', [1 1 1]);
set(hhhh, 'String', 'A', 'Position', [0.06 0.9 0.05 0.05], 'Fontsize', 14);

hhhh = uicontrol('Style','text');
set(hhhh, 'Units', 'Normalized', 'BackgroundColor', [1 1 1]);
set(hhhh, 'String', 'B1', 'Position', [0.06 0.67 0.05 0.05], 'Fontsize', 14);

hhhh = uicontrol('Style','text');
set(hhhh, 'Units', 'Normalized', 'BackgroundColor', [1 1 1]);
set(hhhh, 'String', 'B2', 'Position', [0.06 0.45 0.05 0.05], 'Fontsize', 14);

hhhh = uicontrol('Style','text');
set(hhhh, 'Units', 'Normalized', 'BackgroundColor', [1 1 1]);
set(hhhh, 'String', 'B3', 'Position', [0.06 0.22 0.05 0.05], 'Fontsize', 14);


plotcprof(Snodiff, 1, c1);
plotcprof(Sdiff, 2, c2);
figure(hh);

% Make some text on top
subplot(4,6,2:6);
separator0 = '          ';
separator = '               ';
separator2 = '             ';
separator3 = '                        ';
fmytext = [separator0, 'j^{Na}', separator, 'j^{K}', separator, 'j^{Ca}', separator2, ' j^{X}', separator3, 'i/F'];
haxis = axis;
text(haxis(1),haxis(4)*1.2, fmytext, 'Fontsize', 16)


figure(h);
subplot(1,2,1);
title(['A'], 'Fontsize', 14);
subplot(1,2,2);
xlabel('i^M (nA)', 'Fontsize', 12);
title(['B'], 'Fontsize', 14);

figure(hv);
subplot(1,3,1);
title(['A'], 'Fontsize', 14);
subplot(1,3,2);
title(['B (diff. off)'], 'Fontsize', 14);
xlabel('V (mV)', 'Fontsize', 12);
subplot(1,3,3);
title(['C: (diff. on)'], 'Fontsize', 14);
xlabel('V (mV)', 'Fontsize', 12);


figure(h2);
for i = 1:5
subplot(2,5,i);
title(['A',num2str(i)], 'Fontsize', 14);
subplot(2,5,5+i);
title(['B',num2str(i)], 'Fontsize', 14);
end

subplot(2,5,7); xlabel('[Na^+] (mM)', 'Fontsize', 12);
subplot(2,5,8); xlabel('[K^+] (mM)', 'Fontsize', 12);
subplot(2,5,9); xlabel('[Ca^{2+}] (mM)', 'Fontsize', 12);
subplot(2,5,10); xlabel('[X^{-}] (mM)', 'Fontsize', 12);


function plotcprof(S, crow, Dcol);
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
tv = S.Simdata.t;
cNa = S.Simdata.cNa;
cK = S.Simdata.cK;
cCa = S.Simdata.cCa;
cX = S.Simdata.cX;
V = S.Simdata.V;
diffconsts = S.diffconsts; 
lambda_o = diffconsts(1);
D_K = diffconsts(2);
D_Na = diffconsts(3);
D_Ca = diffconsts(4);
D_X = diffconsts(5);
diffon = S.diffon;

jk = S.Neurondata.jk;
jna = S.Neurondata.jna;
jca = S.Neurondata.jca;
jx = S.Neurondata.jx;
icap = S.Neurondata.icap;
%isyn = S.Neurondata.isyn;
imemb = S.Neurondata.imemb;
times = S.Neurondata.times; 

tv = tv(2:end-1); % just to eliminate some endpoint bugs
V = V(2:end-1,:);
V = [V(1,:); V(:,:); V(end,:)]';
cK = cK'; cNa = cNa'; cCa = cCa'; cX = cX'; 
Npsi = 0.0258;

% Extracellular fluxes: (mol/s)
jKEd = - diffon*Avox*( D_K*(cK(2:end,:) - cK(1:end-1,:))/deltax);
jKEf = - Avox*D_K/Npsi*(cK(2:end,:) + cK(1:end-1,:))/2.*(V(2:end,:) - V(1:end-1,:))/deltax; % N+1 elements
jNaEd =  - diffon*Avox*( D_Na*(cNa(2:end,:) - cNa(1:end-1,:))/deltax);
jNaEf = - Avox*D_Na/Npsi*(cNa(2:end,:) + cNa(1:end-1,:))/2.*(V(2:end,:) - V(1:end-1,:))/deltax; % N+1 elements
jCaEd = - diffon*Avox*( D_Ca*(cCa(2:end,:) - cCa(1:end-1,:))/deltax);
jCaEf = - Avox*2*D_Ca/Npsi*(cCa(2:end,:) + cCa(1:end-1,:))/2.*(V(2:end,:) - V(1:end-1,:))/deltax; % N+1 elements
jXEd =  - diffon*Avox*( D_X*(cX(2:end,:) - cX(1:end-1,:))/deltax);
jXEf = + Avox*D_X/Npsi*(cX(2:end,:) + cX(1:end-1,:))/2.*(V(2:end,:) - V(1:end-1,:))/deltax; % N+1 elements

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

comps = 1:size(cNa,2); % number of compartments
timepts = round(linspace(2,length(t),6)); % Selected time pts.

%avgisyn = mean(isyn, 2);
%avgisyn = [0; avgisyn; 0];
ints = 5;
intt = 7/5;
imtpts = 0:intt:7;

inti = round(length(imemb)/ints);
avgimemb = zeros(ints, size(imemb,1)+2);

for i = 1:ints;
    myint = ((i-1)*inti+1):(i*inti-1);
    avgimemb(i,:) = [0;mean(imemb(:,myint),2);0];

end
zeroline = zeros(1,size(imemb,1)+2);


Mycolors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880];

%%%%%%%%% PLOT SPATIAL PROFILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% FIGURE %%%%%%%
%%%% Spatial profiles %%%%%%
h = findobj('Tag', 'vprofiles');
figure(h);

subplot(1,2,1);
hold on;
% make a drawing of the neuron
plot(3,3,'^', 'markersize', 20, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
plot([3 1.7], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol);
plot([3 4.3], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol)
plot([3 3], [3 13], '-', 'Linewidth', 4, 'Color', Dcol)
plot([3 1.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
plot([3 4.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
%plot(3,3,'^', 'markersize', 8, 'Linewidth', 8, 'Color', Dcol);

axis([0 6 1 15]);
%set(gca, 'Yticklabel', '')
set(gca, 'Xticklabel', '')

subplot(1,2,2); 
set(gca, 'ColorOrder', Mycolors);      
hold on;
%plot(avgisyn*1e9, comps, 'Color', [0.5 0.5 0.5]);
plot(avgimemb*1e9, comps);
plot(zeroline, comps, 'k:');
%haxis = axis; haxis(1) = -0.5; haxis(2) = 1; haxis(3) = 1; haxis(4) = 15; axis(haxis);
Lfroms = num2str(imtpts(1:end-1)');
Ltos = num2str(imtpts(2:end)');
mylegend = [Lfroms, ['     ']', ['-----']', ['     ']', Ltos, ['sssss']'];
legend(mylegend);
set(gca, 'Yticklabel', '')



%%%% FIGURE %%%%%%%
%%%% Spatial profiles ONLY CONCENTRATIONS%%%%%%
h2 = findobj('Tag', 'cprofiles2');
figure(h2);

subplot(2,5,5*(crow-1)+1); 
hold on;
% make a drawing of the neuron
plot(3,3,'^', 'markersize', 10, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
plot([3 1.7], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol);
plot([3 4.3], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol)
plot([3 3], [3 13], '-', 'Linewidth', 4, 'Color', Dcol)
plot([3 1.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
plot([3 4.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
axis([0 6 1 15]);
%set(gca, 'Yticklabel', '')
set(gca, 'Xticklabel', '')

subplot(2,5,5*(crow-1)+2); 
plot(cNa(timepts,:)', comps);
set(gca, 'Yticklabel', '')
haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);

subplot(2,5,5*(crow-1)+3); 
plot(cK(timepts,:)', comps);
haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);
set(gca, 'Yticklabel', '')

subplot(2,5,5*(crow-1)+4); 
plot(cCa(timepts,:)', comps);
haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);
set(gca, 'Yticklabel', '')

subplot(2,5,5*(crow-1)+5); 
plot(cX(timepts,:)', comps);
haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);
set(gca, 'Yticklabel', '')

if crow == 1
    legt2 = round(10*t(timepts(2:end)))/10; 
    legt1 = [0; legt2(1:end-1)];
    legt = [0;legt2];

    subplot(2,5,5); 
    mylegend = [num2str(legt), ['ssssss']'];
    legend(mylegend);
    set(gca, 'Yticklabel', '');
end



%%%% SPATIAL PROFILES ONLY V
hv = findobj('Tag', 'vvprofiles');
figure(hv);
subplot(1,3,1)
hold on;
Dcol = 'k';
% make a drawing of the neuron
plot(3,3,'^', 'markersize', 20, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
plot([3 1.7], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol);
plot([3 4.3], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol)
plot([3 3], [3 13], '-', 'Linewidth', 4, 'Color', Dcol)
plot([3 1.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
plot([3 4.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
axis([0 6 1 15]);
%set(gca, 'Yticklabel', '')
set(gca, 'Xticklabel', '')

if crow == 1
    subplot(1,3,2)
end
if crow == 2
    subplot(1,3,3)
end
avgV = zeros(15,length(timepts)-1);
for ii = 2:length(timepts)
avgV(:,ii-1) = mean(V(:,timepts(ii-1):timepts(ii)),2);
end
hold on; 
plot(1000*avgV', comps);
    plot(zeroline, comps, 'k:');

set(gca, 'Yticklabel', '')
haxis = axis; haxis(1) = -1.5; haxis(2) = 0.2; haxis(3) = 1; haxis(4) = 15; axis(haxis);
%haxis = axis; haxis(3) = 1; haxis(4) = 15; axis(haxis);
set(gca, 'Yticklabel', '')

if 1 
if crow == 1
    subplot(1,3,3)
    avgV = mean(V(:,:),2);
    hold on;
    plot(1000*avgV', comps, 'k', 'LineWidth', 1.5);
%    plot(zeroline, comps, 'k:');
    set(gca, 'Yticklabel', '')
end    
end

if crow == 2
    subplot(1,3,3)
    legt2 = round(10*t(timepts(2:end)))/10;
    legt1 = [0; legt2(1:end-1)];
    mylegend = [num2str(legt1), ['     ']', ['-----']', ['     ']', num2str(legt2), ['sssss']'];
    mylegend = ['  diff. off '; mylegend];
    legend(mylegend);
    plot(zeroline, comps, 'k:');
    set(gca, 'Yticklabel', '')
end




%%%% FIGURE %%%%%%%
%%%% SPATIAL PROFILE OF FLUXES
avgjNaEd = zeros(14,length(timepts)-1);
avgjNaEf = zeros(14,length(timepts)-1);
avgjKEd = zeros(14,length(timepts)-1);
avgjKEf = zeros(14,length(timepts)-1);
avgjCaEd = zeros(14,length(timepts)-1);
avgjCaEf = zeros(14,length(timepts)-1);
avgjXEd = zeros(14,length(timepts)-1);
avgjXEf = zeros(14,length(timepts)-1);


for ii = 2:length(timepts)
avgjNaEd(:,ii-1) = mean(jNaEd(:,timepts(ii-1):timepts(ii)),2);
avgjNaEf(:,ii-1) = mean(jNaEf(:,timepts(ii-1):timepts(ii)),2);
avgjKEd(:,ii-1) = mean(jKEd(:,timepts(ii-1):timepts(ii)),2);
avgjKEf(:,ii-1) = mean(jKEf(:,timepts(ii-1):timepts(ii)),2);
avgjCaEd(:,ii-1) = mean(jCaEd(:,timepts(ii-1):timepts(ii)),2);
avgjCaEf(:,ii-1) = mean(jCaEf(:,timepts(ii-1):timepts(ii)),2);
avgjXEd(:,ii-1) = mean(jXEd(:,timepts(ii-1):timepts(ii)),2);
avgjXEf(:,ii-1) = mean(jXEf(:,timepts(ii-1):timepts(ii)),2);
end

% Extracellular currents: (A)
avgiEd = F*(2*avgjCaEd + avgjKEd + avgjNaEd - avgjXEd);
avgiEf = F*(2*avgjCaEf + avgjKEf + avgjNaEf - avgjXEf);
avgiEd = avgiEd/F; avgiEf = avgiEf/F; 
avgiE = avgiEd + avgiEf;

fcomps = comps(1:end-1); % One less

Mycolors = [0    0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880];


%%%%% FIELD STUFF
hh = findobj('Tag', 'jprofiles');
figure(hh);

if crow == 2
    subplot(4,6,7);
    hold on;
    % make a drawing of the neuron
    Dcol = 'b';
    plot(3,3,'^', 'markersize', 10, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
    plot([3 1.7], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol);
    plot([3 4.3], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol)
    plot([3 3], [3 13], '-', 'Linewidth', 4, 'Color', Dcol)
    plot([3 1.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
    plot([3 4.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
    axis([0 6 1 15]);
    set(gca, 'Xticklabel', '')
    
    
    haxis = [-2.0000e-14 1.5400e-13 0.5 14.5];
    %xsep = 1.5e-14;
%    xsep = 2.8e-14;
    xsep = 3.2e-14;
    subplot(4,6,8:12);
    set(gca, 'ColorOrder', Mycolors);      
    hold on;
    plot(avgjNaEf, fcomps);
    plot([0 0], [1 15], 'k--');
    
    plot(avgjKEf+xsep, fcomps);
    plot([xsep xsep], [1 15], 'k--');
    
    plot(avgjCaEf+2*xsep, fcomps);
    plot([2*xsep 2*xsep], [1 15], 'k--');
    
    plot(avgjXEf+3*xsep, fcomps);
    plot([3*xsep 3*xsep], [1 15], 'k--');
    
    plot(avgiEf+4.5*xsep, fcomps);
    plot([4.5*xsep 4.5*xsep], [1 15], 'k--');
    
    axis(haxis);
    set(gca, 'Yticklabel', '')
    set(gca, 'Xticklabel', '')
    set(gca, 'Xtick',1000)
    set(gca, 'Visible', 'off')
    
    
    plot([-1 1],[0.5 0.5], 'k-');
    plot([-1 1],[14.5 14.5], 'k-');
    text(-2e-14, 11.5, 'field', 'Fontsize', 14)
    
    %set(gcf,'Position',[ 27 82 1001 733]);
    
    
    %%%%% DIFFUSIVE STUFF
    subplot(4,6,13);
    hold on;
    Dcol = 'b';
    % make a drawing of the neuron
    plot(3,3,'^', 'markersize', 10, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
    plot([3 1.7], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol);
    plot([3 4.3], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol)
    plot([3 3], [3 13], '-', 'Linewidth', 4, 'Color', Dcol)
    plot([3 1.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
    plot([3 4.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
    axis([0 6 1 15]);
    %set(gca, 'Yticklabel', '')
    set(gca, 'Xticklabel', '')
    
    subplot(4,6,14:18);
    set(gca, 'ColorOrder', Mycolors);      
    hold on;
    
    plot(avgjNaEd, fcomps); %title('jNaEf')
    plot([0 0], [1 15], 'k--');

    plot(avgjKEd+xsep, fcomps); %title('jKEf')
    plot([xsep xsep], [1 15], 'k--');
    
    plot(avgjCaEd+2*xsep, fcomps); %title('jCaEf')
    plot([2*xsep 2*xsep], [1 15], 'k--');
    
    plot(avgjXEd+3*xsep, fcomps); %title('jXEf')
    plot([3*xsep 3*xsep], [1 15], 'k--');
    
    plot(avgiEd+4.5*xsep, fcomps); %title('iEf')
    plot([4.5*xsep 4.5*xsep], [1 15], 'k--');
    
    plot([-1 1],[0.5 0.5], 'k-');
    plot([-1 1],[14.5 14.5], 'k-');
    text(-2e-14, 11.5, 'diff', 'Fontsize', 14)
    
    axis(haxis);
    set(gca, 'Yticklabel', '')
    set(gca, 'Xticklabel', '')
    set(gca, 'Xtick',1000)
    set(gca, 'Visible', 'off')
    
    
    %%%%% TOTAL STUFF
    subplot(4,6,19);
    hold on;
    % make a drawing of the neuron
    Dcol = 'b';
    plot(3,3,'^', 'markersize', 10, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
    plot([3 1.7], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol);
    plot([3 4.3], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol)
    plot([3 3], [3 13], '-', 'Linewidth', 4, 'Color', Dcol)
    plot([3 1.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
    plot([3 4.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
    axis([0 6 1 15]);
    %set(gca, 'Yticklabel', '')
    set(gca, 'Xticklabel', '')
    
    subplot(4,6,20:24);
    set(gca, 'ColorOrder', Mycolors);      
    hold on;
    plot(avgjNaEd + avgjNaEf, fcomps); %title('jNaE')
    plot([0 0], [1 15], 'k--');
    plot(avgjKEd  + avgjKEf + xsep, fcomps); %title('jKE')
    plot([xsep xsep], [1 15], 'k--');
    plot(avgjCaEd + avgjCaEf +2*xsep, fcomps); %title('jCaE')    
    plot([2*xsep 2*xsep], [1 15], 'k--');
    plot(avgjXEd + avgjXEf + 3*xsep, fcomps); %title('jXE')
    plot([3*xsep 3*xsep], [1 15], 'k--');
    plot(avgiE + 4.5*xsep, fcomps); %title('iE')
    plot([4.5*xsep 4.5*xsep], [1 15], 'k--');    
    plot([-1 1],[0.5 0.5], 'k-');
    plot([-1 1],[14.5 14.5], 'k-');
    text(-2e-14, 11.5, 'tot', 'Fontsize', 14)
    

    
    axis(haxis);
    set(gca, 'Yticklabel', '')
    set(gca, 'Xticklabel', '')
    set(gca, 'Xtick',1000)
    set(gca, 'Visible', 'off')
    
end
if crow == 1
    %%%%% TOTAL STUFF
    subplot(4,6,1);
    hold on;     
    haxis = [-2.0000e-14 1.5400e-13 0.5 14.5];
    xsep = 3.2e-14;

    % make a drawing of the neuron
    Dcol = 'r';
    plot(3,3,'^', 'markersize', 10, 'Linewidth', 0.5, 'MarkerFaceColor', Dcol ,'Color', Dcol);
    plot([3 1.7], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol);
    plot([3 4.3], [3 1.7], '-', 'Linewidth', 2, 'Color', Dcol)
    plot([3 3], [3 13], '-', 'Linewidth', 4, 'Color', Dcol)
    plot([3 1.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
    plot([3 4.5], [13 14], '-', 'Linewidth', 2, 'Color', Dcol)
    axis([0 6 1 15]);
    %set(gca, 'Yticklabel', '')
    set(gca, 'Xticklabel', '')
    
    subplot(4,6,2:6);
    set(gca, 'ColorOrder', Mycolors);      
    hold on;
    plot(avgjNaEd + avgjNaEf, fcomps); %title('jNaE')
    plot([0 0], [1 15], 'k--');
    
    plot(avgjKEd  + avgjKEf + xsep, fcomps); %title('jKE')
    plot([xsep xsep], [1 15], 'k--');
    
    plot(avgjCaEd + avgjCaEf +2*xsep, fcomps); %title('jCaE')
    plot([2*xsep 2*xsep], [1 15], 'k--');
    
    plot(avgjXEd + avgjXEf + 3*xsep, fcomps); %title('jXE')
    plot([3*xsep 3*xsep], [1 15], 'k--');
    
    plot(avgiE + 4.5*xsep, fcomps); %title('iE')
    plot([4.5*xsep 4.5*xsep], [1 15], 'k--');
    
    legt2 = round(10*t(timepts(2:end)))/10;
    legt1 = [0; legt2(1:end-1)];
    mylegend = [num2str(legt1), ['     ']', ['-----']', ['     ']', num2str(legt2), ['sssss']'];
    
    hmfr = legend(mylegend, 'orientation', 'horizontal');
    set(hmfr, 'Fontsize', 12);
    plot([-1 1],[0.5 0.5], 'k-');
    plot([-1 1],[14.5 14.5], 'k-');    
    text(-2e-14, 11.5, 'field', 'Fontsize', 14)

%scale bar
% For 10 mu mol/m^2/s
%    barlength = 10e-6*Avox
%    plot([4.4*xsep 4.4*xsep+barlength], [2 2] , '-k', 'LineWidth', 3); %scalebar 10 mu mol/m^2/s
%    text(4.4*xsep, 3.5, '10 {\mu}mol/(m^2s)' );

    barlength = 10e-6*Avox
    plot([2.02*xsep 2.02*xsep+barlength], [10 10] , '-k', 'LineWidth', 2.5); %scalebar 10 mu mol/m^2/s
    plot([1.98*xsep-barlength 1.98*xsep], [10 10] , '-k', 'LineWidth', 2.5); %scalebar 10 mu mol/m^2/s

    text(2.01*xsep, 8, '10 {\mu}mol/(m^2s)', 'Fontsize',12);
    text(2.01*xsep + 0.2*barlength, 12, '+', 'Fontsize', 13 );
    text(1.99*xsep-0.7*barlength, 12, '-', 'Fontsize', 13 );
    
    axis(haxis);
    set(gca, 'Yticklabel', '')
    set(gca, 'Xticklabel', '')
    set(gca, 'Xtick',1000)
    set(gca, 'Visible', 'off')
end


