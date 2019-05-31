function [T,Y] = nvoxels(times, jna, jk, jca, jx, icap, imemb, deltax, Avox, VECS, N, Ns, basal, diffconsts, constsigma, diffon, nrnon)
% Geir Halnes, oct. 2014

% Initial conditions (in mM = mol/m^3);
bcs = basal; % basal concentrations
cK = basal(1)*ones(1,N);
cNa = basal(2)*ones(1,N);
cCa = basal(3)*ones(1,N);
cX = basal(4)*ones(1,N);
V = basal(5)*ones(1,N);
ics = [cK cNa cCa cX V]; % Must be vector, length N*Ns

simstart = times(1);
simstop = times(end);
if nrnon
    shutup = simstop;
else
shutup = simstop/2;
end
mymaxstep = 1e-3;
options = odeset('RelTol',1e-3,'AbsTol',[1e-4], 'MaxStep', mymaxstep)
[T,Y] = ode45(@(t,x)threevox(t,x,times, jk, jna, jca, jx, icap, imemb, deltax, Avox, VECS, diffconsts, N, Ns, bcs, constsigma, diffon, shutup),[simstart simstop],ics,options);


function [dx] = threevox(t,x, times, jk, jna, jca, jx, icap, imemb, deltax, Avox, VECS, diffconsts, N, Ns, bcs, constsigma, diffon, shutup)
global Vout
global Iout
global tcounter

if rand<0.001
disp(['Biological time: ', num2str(t),  ' of 84 seconds' ]);
disp(['Processing time: ', num2str(toc),  ' seconds' ]);
end

% Some constants:
F= 96485.3365; % C/mol
T = 300; % K
R = 8.3; % J/mol/K
Npsi = R*T/F; % V
lambda_o = diffconsts(1);
D_K = diffconsts(2);
D_Na = diffconsts(3);
D_Ca = diffconsts(4);
D_X = diffconsts(5);


% Read out membrane currents from data file:
if t>shutup
    hoho = 0; % For removing neuronal output after time shutup
else
    hoho = 1;
end

% DO SOME INPERPOLATION HERE
index = min(find(times>=t)); % To get input data at right time.
if index < 1
    index = 1;
end
jk = hoho*jk(:,index)';
jna = hoho*jna(:,index)';
jca = hoho*jca(:,index)';
jx = hoho*jx(:,index)';
icap = hoho*icap(:,index)';
imemb = hoho*imemb(:,index)';


% Constant background
cK0 = bcs(1);
cNa0 = bcs(2);
cCa0 = bcs(3);
cX0 = bcs(4);
V0 = bcs(5);

% Unwrap variables
xx = reshape(x,[N, Ns])';
cK = [cK0 xx(1,:) cK0]; % N+2 elements (N with neural output, 2 edges)
cNa = [cNa0 xx(2,:) cNa0];
cCa = [cCa0 xx(3,:) cCa0];
cX = [cX0 xx(4,:) cX0];
%V = [V0 xx(5,:) V0];


% Diffusive currents (NP: idiff = F*z_k*D_k*dc_k/dx, muliply with Avox to get net current)
idiff = -Avox/deltax*F*( D_Na*(cNa(2:end) - cNa(1:end-1)) + D_K*(cK(2:end) - cK(1:end-1))...
    + 2*D_Ca*(cCa(2:end) - cCa(1:end-1)) - D_X*(cX(2:end) - cX(1:end-1))); % N+1 elements
idiff = idiff*diffon;


% Define conductivities:
sigma = F/2/Npsi*(D_Na.*(cNa(1:end-1) + cNa(2:end)) + D_K*(cK(1:end-1) + cK(2:end))...
    + 4*D_Ca*(cCa(1:end-1) + cCa(2:end)) + D_X*(cX(1:end-1) + cX(2:end))); % % N+1 elements (ohm*m)^-1
if constsigma
    sigma(1:end) = F/Npsi*(D_Na.*cNa0 + D_K*cK0 + 4*D_Ca*cCa0 + D_X*cX0);
end


%%%%%% FIND V
% Solve equation set on the form:
% M_mn * V_n = A_n ---> V_n = inv(M_mn)*A_n;

% Define Matrix M_mn:
MNdiag0 = 1;
MNdiagmid = -(sigma(1:end-1) + sigma(2:end)); %N elements
MNdiagend = -sigma(end);
MNdiag = [MNdiag0 MNdiagmid MNdiagend]; % N+2 elements

MNabove0 = 0;
MNabovemid = sigma(2:end); %N elements
MNabove = [MNabove0 MNabovemid]; % N+1 elements

MNbelowmid = sigma(1:end-1); % N elements
MNbelowend = sigma(end);
MNbelow = [MNbelowmid MNbelowend]; % N+1 elements

Mmn = diag(MNdiag) + diag(MNbelow,-1) + diag(MNabove,1); % (N+2)*(N+2)

%Ionic membrane currents (ik = F*z_k*j_k)
iion = F*(jna + jk + 2*jca - jx); % N elements

% Define matrix A_n
AAnmid = (-icap - iion - idiff(1:end-1) + idiff(2:end))*deltax/Avox; % N elements
AAn0 = 0; % 1 elements: V0 = 0
AAnend = - idiff(end)*deltax/Avox; % 1 elements
AAn = [AAn0 AAnmid AAnend];

%%%% Extracellular voltage
VVn = inv(Mmn)*AAn'; % N elements
V = VVn';

%%%% Update global variables
Vout(tcounter,:) = [t V];
Iout(tcounter,:) = [t 0 imemb 0];
tcounter = tcounter + 1;


% Extracellular fluxes: (mol/s)
jKE = - Avox*( diffon*D_K*(cK(2:end) - cK(1:end-1))/deltax ...
+ D_K/Npsi*(cK(2:end) + cK(1:end-1))/2.*(V(2:end) - V(1:end-1))/deltax); % N+1 elements

jNaE = - Avox*( diffon*D_Na*(cNa(2:end) - cNa(1:end-1))/deltax ...
+ D_Na/Npsi*(cNa(2:end) + cNa(1:end-1))/2.*(V(2:end) - V(1:end-1))/deltax); % N+1 elements

jCaE = - Avox*( diffon*D_Ca*(cCa(2:end) - cCa(1:end-1))/deltax ...
+ 2*D_Ca/Npsi*(cCa(2:end) + cCa(1:end-1))/2.*(V(2:end) - V(1:end-1))/deltax); % N+1 elements

jXE = - Avox*( diffon*D_X*(cX(2:end) - cX(1:end-1))/deltax ...
- D_X/Npsi*(cX(2:end) + cX(1:end-1))/2.*(V(2:end) - V(1:end-1))/deltax); % N+1 elements

% Divide by ECS volume to get changes in ion conc (mol/m^3/s = mM/s)
dKE = ( -(jKE(2:end)-jKE(1:end-1)) + jk)/VECS;
dNaE = ( -(jNaE(2:end)-jNaE(1:end-1)) + jna)/VECS; % N elements
dCaE = ( -(jCaE(2:end)-jCaE(1:end-1)) + jca)/VECS;
dXE = ( -(jXE(2:end)-jXE(1:end-1)) + jx)/VECS;
dVE = zeros(size(dXE)); % Not used
dx = [dKE dNaE dCaE dXE dVE]';
