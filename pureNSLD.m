%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function: pureNSLD
%%% A database for relevant NSLDs of the hydration project
%%% Input: compound: the compound that you want to calculate the NSLD for
%%%        options: D2O, H2O, PNIPAM, PDMAPS, PHPA
%%% Output: the scattering length density in 1/A^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NSLD = pureNSLD(compound)
% All NSLD in 1/A^2
if strcmp(compound, 'D2O')
    NSLD = 6.39E-06; %with density 1.11 g/mL or cm^3
elseif strcmp(compound, 'H2O')
    NSLD = -5.61E-07; %with density 1 g/mL
elseif strcmp(compound, 'PNIPAM')  
    NSLD = 0.777e-6; %with density 1.05 g/mL % I changed this from 0.814 (from density 1.1 g/mL)
elseif strcmp(compound, 'PDMAPS')
    NSLD = 1.057e-6; % with density 1.37 g/mL % I changed this from 1.202 (from wrong molecular formula C12 instead of C11)
elseif strcmp(compound, 'PHPA')
    NSLD = 1.0698e-6; % with density 1.16 g/mL
end