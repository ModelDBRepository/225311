
function runit(times, jna, jk, jca, jx, icap, imemb, deltax, Avox, VECS, N, Nnrns, constsigma, diffon, nrnon, Datafile, Filename);
% Run simulation
% Prepares variables for simulation
% Runs simulation by calling function NVOXELS
% Collects data in structure S and saves it to file.

global tcounter
global Vout
global Iout

% Allocate memory. 
Nentries = length(times)*0.6
% Note: This is more or less the number of time points that Vout is
% estimated at in the RK-based numerical scheme.
tcounter = 1;
Vout = zeros(Nentries+50000,N+3); % N for nrns, 2 edges, 1 for time
Iout = zeros(Nentries+50000,N+3); % Added 50.000 elements just to make sure we have enough

% Diffusion constants
% Nano and Molecular Electronics Handbook: Sergey Edward Lyshevski. CRC Press, Taylor & Francis group 2007
lambda_o = 1.6; % ECS tortuousity, Chen & Nicholson 2000;
D_K = 1.96e-9/lambda_o^2; % Diffusion coefficients (m^2/s);
D_Ca=0.71e-9/lambda_o^2;  % 0,8 used by Gardner2015, 0.6 by Lewin2012, 0.71 in Lyshevski
D_Na=1.33e-9/lambda_o^2; 
D_X = 2.0e-9/lambda_o^2; % Using diffusion constant for Cl-
diffconsts = [lambda_o, D_K, D_Na, D_Ca, D_X];

% Initial conditions
Ns = 5; % 4 ion species +  V
cK0 = 3;
cNa0 = 150;
cCa0 = 1.4;
cX0 = cK0 + cNa0 + 2*cCa0;
V0 = 0;
basal = [cK0; cNa0; cCa0; cX0; V0];


[t,Y] = nvoxels(times, jna, jk, jca, jx, icap, imemb, deltax, Avox, VECS, N, Ns, basal, diffconsts, constsigma, diffon, nrnon);
%[t,Y] = nvoxels_nodiff(times, jna, jk, jca, jx, icap, deltax, Avox, VECS, N, Ns, basal);

% Unwrap variables:
Y = reshape(Y, [length(Y) N Ns]);
cKedge = cK0*ones(length(t),1);
cNaedge = cNa0*ones(length(t),1);
cCaedge = cCa0*ones(length(t),1);
cXedge = cX0*ones(length(t),1);
cK = [cKedge, Y(:,:,1), cKedge]; % mol/m^3 (or mM)
cNa = [cNaedge, Y(:,:,2), cNaedge]; % mol/m^3 (or mM)
cCa = [cCaedge, Y(:,:,3), cCaedge]; % mol/m^3 (or mM)
cX = [cXedge, Y(:,:,4), cXedge]; % mol/m^3 (or mM)
clear Y

% Interpolate to get V at same time points as cK:
% Note: With the numerical scheme used here, tv & V comes out at 1.5 times 
% as many interpolation points as t. 

tvr = Vout(:,1);
Vr = Vout(:,2:end);
tir = Iout(:,1);
Ir = Iout(:,2:end);
clear Vout;
clear Iout;

keepind = find(diff(tvr)); % some time steps are zero. Remove these
tvr = tvr(keepind); % "raw" data for global variables
Vr = Vr(keepind,:); % "raw" data
keepind = find(diff(tir)); % some time steps are zero. Remove these
tir = tir(keepind); % "raw" data
Ir = Ir(keepind,:); % "raw" data

V = interp1(tvr,Vr,t); % same #indices as concentration data
Im = interp1(tir,Ir,t);


% Prepare structure of output data
S.Datafile = Datafile;
S.Nionspecies = Ns;
S.Nvox = N+2;
S.diffconstnames = 'lambda_0 D_K D_Na D_Ca D_X';
S.diffconsts = [lambda_o, D_K, D_Na, D_Ca, D_X];
S.basalnames = 'cK0 cNa0 cCa0 cX0 V0';
S.basal = [cK0 cNa0 cCa0 cX0 V0];
S.diffon = diffon;
S.nrnon = nrnon;
S.constsigma = constsigma; 

S.geometry.deltax = deltax;
S.geometry.ECSfrac = 0.2;
S.geometry.Nnrns = Nnrns; 
S.geometry.Avox = Avox;
S.geometry.Vvox = VECS;

S.Neurondata.Nnrns = Nnrns;
S.Neurondata.jk = jk; 
S.Neurondata.jna = jna; 
S.Neurondata.jca = jca; 
S.Neurondata.jx = jx; 
S.Neurondata.icap = icap; 
S.Neurondata.imemb = imemb; 
S.Neurondata.times = times;

S.Simdata.t = t;
S.Simdata.cK = cK; 
S.Simdata.cNa = cNa; 
S.Simdata.cCa = cCa; 
S.Simdata.cX = cX;

S.Simdata.V = V; 
S.Simdata.Im = Im;

S.Simdata.tvr = tvr; % Also save the "raw" global data
S.Simdata.tir = tir;
S.Simdata.Vr = Vr;
S.Simdata.Ir = Ir;


save(Filename, 'S');