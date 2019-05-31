function PLoSdataplot(S)
%%% PLOT THE DATA1SETS
F = 96485.3365; % C/mol
CVF = F*1e9; % Convert flux densities to nA

ik = S.Neurondata.jk*CVF;
ina = S.Neurondata.jna*CVF;
ica = S.Neurondata.jca*2*CVF;
ix = -S.Neurondata.jx*CVF;
icap = S.Neurondata.icap*1e9; %convert A to nA
imembt = S.Neurondata.imemb*1e9;
times = S.Neurondata.times*1000; %convert s to ms
iion = ik + ina + ica + ix; % nA

% total currents
inat = sum(ina,1); % summed over all depths, as fcn of time
ikt = sum(ik,1);
icat = sum(ica,1);
icapt = sum(icap,1);
ixt = sum(ix,1);
iiont = sum(iion,1);
imembtt = sum(imembt,1);

% time interval to plot:
tinterv = [0 7000];
tindexes = find(times>=tinterv(1) & times <=tinterv(2));

% compartments to plot
NN1 = 2; % soma 
NN2 = 7; % trunk
NN3 = 12; % apical
% Note: Neuronal output for 13 compartments (size(ina) = #timepts x 13)
% Later we add a voxel on top and bottom to get 15 compartments
% Then e.g. soma is voxel 3.


h = figure; hold on;

myc = 'k';

set(gcf, 'Color', [1 1 1]);
set(gcf, 'Position', [21 141 1050 697]);

%%% PLOT SOMA
subplot(6,4,1)
plot(times(tindexes), ina(NN1,tindexes), myc);

subplot(6,4,5)
plot(times(tindexes),ik(NN1,tindexes), myc);

subplot(6,4,9)
plot(times(tindexes),ica(NN1, tindexes), myc);

subplot(6,4,13)
plot(times(tindexes),ix(NN1,tindexes), myc);


subplot(6,4,17)
plot(times(tindexes),icap(NN1,tindexes), myc);

subplot(6,4,21)
plot(times(tindexes),imembt(NN1,(tindexes)), myc);
xlabel('t(ms)',  'FontSize', 12);


%%% PLOT TRUNK
subplot(6,4,2)
plot(times(tindexes), ina(NN2,tindexes), myc);

subplot(6,4,6)
plot(times(tindexes),ik(NN2,tindexes), myc);

subplot(6,4,10)
plot(times(tindexes),ica(NN2, tindexes), myc);

subplot(6,4,14)
plot(times(tindexes),ix(NN2,tindexes), myc);

subplot(6,4,18)
plot(times(tindexes),icap(NN2,tindexes), myc);

subplot(6,4,22)
plot(times(tindexes),imembt(NN2,(tindexes)), myc);
xlabel('t(ms)',  'FontSize', 12);


%%% PLOT APICAL
subplot(6,4,3)
plot(times(tindexes), ina(NN3,tindexes), myc);

subplot(6,4,7)
plot(times(tindexes),ik(NN3,tindexes), myc);

subplot(6,4,11)
plot(times(tindexes),ica(NN3, tindexes), myc);

subplot(6,4,15)
plot(times(tindexes),ix(NN3,tindexes), myc);

subplot(6,4,19)
plot(times(tindexes),icap(NN3,tindexes), myc);

subplot(6,4,23)
plot(times(tindexes),imembt(NN3,(tindexes)), myc);
xlabel('t(ms)',  'FontSize', 12);


%%% PLOT TOTAL
subplot(6,4,4)
plot(times(tindexes), inat(tindexes), myc);

subplot(6,4,8)
plot(times(tindexes),ikt(tindexes), myc);

subplot(6,4,12)
plot(times(tindexes),icat(tindexes), myc);

subplot(6,4,16)
plot(times(tindexes),ixt(tindexes), myc);

subplot(6,4,20)
plot(times(tindexes),icapt(tindexes), myc);

subplot(6,4,24)
plot(times(tindexes),imembtt(tindexes), myc);
xlabel('t(ms)',  'FontSize', 12);


%%% Remove xticks
for i = 1:20
    subplot(6,4,i);
    set(gca, 'XtickLabel', '');
end


%%% Get right time-axes
for i = 1:24
    subplot(6,4,i)
    haxis = axis;
    haxis(1) = tinterv(1); haxis(2) = tinterv(2);
    axis(haxis);
end



%%% Insert titles
subplot(6,4,1); title('A    Soma', 'FontSize',14)
subplot(6,4,2); title('B    Trunk', 'FontSize',14)
subplot(6,4,3); title('C    Apical', 'FontSize',14)
subplot(6,4,4); title('D    Total', 'FontSize',14)

% Insert y-labels
subplot(6,4,1); ylabel('i_{Na}(nA)', 'FontSize', 12);
subplot(6,4,5); ylabel('i_K(nA)', 'FontSize', 12);
subplot(6,4,9); ylabel('i_{Ca}(nA)', 'FontSize', 12);
subplot(6,4,13); ylabel('i_X(nA)', 'FontSize', 12);
subplot(6,4,17); ylabel('i_{cap}(nA)', 'FontSize', 12);
subplot(6,4,21); ylabel('i_{tot}(nA)', 'FontSize', 12);

%%% Number rows
vstep = 1/7;
for i = 1:6
hhhh = uicontrol('Style','text');
set(hhhh, 'Units', 'Normalized', 'BackgroundColor', [1 1 1]);
strnr = num2str(i);
set(hhhh, 'String', strnr, 'Position', [0.04 0.87-vstep*(i-1) 0.03 0.05], 'Fontsize', 14);
end



%figure;
%plot(mean(isyn,2), 1:13)
