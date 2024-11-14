%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function: mono_volume
%%% Calculates the volume of each polymeric unit of PNIPAM, PHPA, and
%%% PDMAPS in Angstrom^3, as well as the degree of polymerization
%%% Input: (1) polymer: 'PNIPAM','PHPA','PDMAPS'
%%%        (2) MWp: polymer molecular weight in g/mol
%%% Output: (1) v, volume [Angstrom^3]
%%%         (2) N, degree of polymerization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v,N] = mono_volume(polymer,MWp)
% Assume ideal 
Na = 6.022e23; %Avogadro's number
if strcmp(polymer,'PNIPAM')
    M = 113.16; %NIPAM, g/mol
    rho = 1.05; %PNIPAM, g/mL
elseif strcmp(polymer,'PHPA')
    M = 130.14; %HPA, g/mol
    rho = 1.16; %PHPA, g/mL
elseif strcmp(polymer,'PDMAPS')
    M = 279.35; %DMAPS, g/mol
    rho = 1.37; %PDMAPS, g/mL
elseif strcmp(polymer,'POEGA')
    M = 150.86; %OEGA, g/mol, based on C6.9 H11.9 O3.5 (Chang 2014)
    rho = 1.05; %POEGA, g/mL
end

rho = rho/1e21; %convert to g/nm^3
v = M/(rho*Na); %nm^3
N = MWp/M; %degree of polymerization

% Convert to A^3
v = v*1000;