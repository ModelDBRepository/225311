% combinemattomat_fixeddt.m
% A MATLAB script for combining the files saved by running "calcsumcurr_manyareagsynmediumtau_parts_fixeddt.py 20 0.025 0.000042 10000 10000 2 myseed 200"
% Expects that the variable "myseed" has been initialized
% Tuomo Maki-Marttunen, 2014-2016

synloctype = 2;
nsegs = 20;
dt = 0.025;
tstop = 10000;
Nsynlocs = 10000;
singleSimT = 200;

syngmaxes = [nan 0.000042 nan];
syngmax = syngmaxes(synloctype);

Nsims = floor(1.0*tstop/singleSimT+0.9999);
Nparts = 13;
ina = [];
ik = [];
ica = [];
ih = [];
il = [];
VtimesA = [];
imemb = [];
Vsoma = [];
icap = [];
times = [];
dt_int = 0.1;

tic
for isim=1:Nsims
  disp(['isim=' num2str(isim)]);
  disp(['Loading isim=' num2str(isim) ', myseed=' num2str(myseed) ', toc=' num2str(toc)]);
  A = load(['currsums_parts_' num2str(Nsynlocs) 'areagsynsmediumtau_fixeddt_type' num2str(synloctype) '_amp' num2str(syngmax) '_tstop' num2str(tstop) '.0_nseg' num2str(nsegs) '_dt' num2str(dt) '_seed' num2str(myseed) '_sim' num2str(isim-1) 'x200.0.mat']);

  times = [times;A.times];
  ina = [ina, A.ina];
  ik = [ik, A.ik];
  ica = [ica, A.ica];
  ih = [ih, A.ih];
  il = [il, A.il];
  VtimesA = [VtimesA, A.VtimesA];
  imemb = [imemb, A.imemb];
  icap = [icap, A.icap];
  Vsoma = [Vsoma; A.Vsoma];

  if isim==1
    ts_syn = A.ts_syn;
    part_syn = A.part_syn;
  end
end

clear A
clear isim
clear Nsims
clear dt_int
clear ts


save(['currsums_parts_' num2str(Nsynlocs) 'areagsynsmediumtau_fixeddt_type' num2str(synloctype) '_amp' num2str(syngmax) '_tstop' num2str(tstop) '.0_nseg' num2str(nsegs) '_dt' num2str(dt) '_seed' num2str(myseed) '_comb200.0.mat']);
