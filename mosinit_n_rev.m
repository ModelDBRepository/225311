function mosinit_n_rev
% Simulates extracellular dynamics for given neuronal output, and plots results.
% 
% Data files: 
%   mosinit_n_rev requires that the file 'revdata.mat' is located in the 
%   current folder. 'revdata.mat' contains neuronal transmembrane data
%   and is used as input time series to the extracellular space.
%
%   mosinit_n_rev generates three data files, representing tree independent 
%   simulations:
%       Sim_diffon_rev.mat: Extracellular diffusion included.
%       Sim_diffoff_rev.mat: Extracellular diffusion set to zero.
%       Sim_dead_rev.mat: Extracellular diffusion included, but neurons
%       'turned off' midways in the simulation.
%   These three files are saved to the current folder. 
%   
%   The user is given two options, dosimulations and doplot (see below).
%   dosimulations = 1 must be chosen the first time this file is run. 
%   Later, dosimulation = 0 can be chosen to plot already saved data-files. 


% Options
dosimulations = 1; % To run extracellular simulations (results saved to file).
doplot = 1; % To plot data


% Some simulation choices
constsigma = 0; % keep conductance constant


% Geometrical parameters
deltax = 100e-6; % Length of cylinder section
ECSfrac = 0.2; % Fraction of tissue being extracellular space
Nnrns = 10; % We simulate 10 neurons
Anrn = 300e-12; % Area per nrn. Linden et al. 2011: 1 nrn per 300 mum^2 of column.
Avox = Anrn*ECSfrac*Nnrns; % Compartment surface area of extracellular space (m^2)
VECS = Avox*deltax; % Compartment volume (m^3)




% Simulate Extracellular Dynamics (3 scenarios) and save to files.
if dosimulations
    % NEURON output data
    Datafile = 'revdata.mat'; % This data files was produced elswhere (folder nrndata) and must be availabe.
    load('revdata.mat');

    tic
    disp('Running simulation: Diffusion on:')
    diffon = 1;
    nrnon = 1;
    Sims1 = 'Sim_diffon.mat';
    runit(ts, jna, jk, jca, jx, icap, imemb, deltax, Avox, VECS, N, Nnrns, constsigma, diffon, nrnon, Datafile, Sims1);
    disp('took seconds:')
    timetook = toc

    disp('Running simulation: Diffusion off')
    diffon = 0;
    nrnon = 1;
    Sims2 = 'Sim_diffoff.mat';
    runit(ts, jna, jk, jca, jx, icap, imemb, deltax, Avox, VECS, N, Nnrns, constsigma, diffon, nrnon, Datafile, Sims2);
    disp('took seconds:')
    timetook = toc

    disp('Running simulation: Neuron turned off midways.')
    diffon = 1;
    nrnon = 0;
    Sims3 = 'Sim_dead.mat';
    runit(ts, jna, jk, jca, jx, icap, imemb, deltax, Avox, VECS, N, Nnrns, constsigma, diffon, nrnon, Datafile, Sims3);
    disp('took seconds:')
    timetook = toc    
else
    disp('No simulations run. Using saved simulation data.')
    Sims1 = 'Sim_diffon_def.mat';
    Sims2 = 'Sim_diffoff_def.mat';
    Sims3 = 'Sim_dead_def.mat';
end


% GENERATE PLOTS
if doplot
    close all;
    S1 = load(Sims1); S1 = S1.S;
    S2 = load(Sims2); S2 = S2.S;
    S3 = load(Sims3); S3 = S3.S;
    PLoSdataplot(S1);
    PLoSdynplot(S1,S2);
    PLoStempplot(S1,S2);
    PLoSfreqplot(3, S1, S2);
    PLoSsilenceplot(S1,S2,S3);
end
