function tissue = makeTissueList(nm)
%function tissueProps = makeTissueList(nm)
%   Returns the tissue optical properties at the wavelength nm:
%       tissueProps = [mua; mus; g]';
%       global tissuenames(i).s
%   Uses 
%       SpectralLIB.mat

%% Load spectral library
load spectralLIB.mat
%   muadeoxy      701x1              5608  double              
%   muamel        701x1              5608  double              
%   muaoxy        701x1              5608  double              
%   muawater      701x1              5608  double              
%   musp          701x1              5608  double              
%   nmLIB         701x1              5608  double              
MU(:,1) = interp1(nmLIB,muaoxy,nm);
MU(:,2) = interp1(nmLIB,muadeoxy,nm);
MU(:,3) = interp1(nmLIB,muawater,nm);
MU(:,4) = interp1(nmLIB,muafat,nm);
MU(:,5) = interp1(nmLIB,muamel,nm);
LOADED = 1;

%% Create tissueList

j=1;
tissue(j).name  = 'escape';
tissue(j).mua   = 0.0001;
tissue(j).mus   = 1.0;
tissue(j).g     = 1.0;

j=2;
tissue(j).name  = 'air';
tissue(j).mua   = 0.001;
tissue(j).mus   = 10;
tissue(j).g     = 1.0;

j=3;
tissue(j).name = 'water';
tissue(j).mua = MU(:,3);
tissue(j).mus = 0.003; %from https://s.campbellsci.com/documents/br/technical-papers/obs_light_absorption.pdf
tissue(j).g  = 1.0; %nog checken!! maar lijkt me wel logisch

j=4;
tissue(j).name  = 'blood';
B       = 1.00;         %blood volume fraction
S       = 0.75;         %oxygen saturation of hemoglobin
W       = 0.95;         %water volume fraction
F       = 0;            %fat volume fraction
M       = 0;            %volume fraction of melasonomes
musp500 = 10;           %reduced scattering coefficient at 500 nm
fray    = 0.0;          %fraction of Rayleigh scattering at 500 nm
bmie    = 1.0;          %scatter power for mie scattering
gg      = 0.90;         %anisotropy factor
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W F M]';
tissue(j).mua = MU*X;   %absorption coefficient, oxyblood, deoxyblood water fat melansomes]
tissue(j).mus = musp/(1-gg);%scattering coefficient
tissue(j).g   = gg;     %scattering anisotropy

j=5;
tissue(j).name  = 'epidermis';
B = 0;
S = 0.75;
W = 0.75;
F = 0;
M = 0.03;
musp500 = 40;
fray    = 0.0;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W F M]';
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

j=6; %waardes vab magnetic resonance guided NIRT of breast 2004, Brooksby, geldig voor 785 nm!!
tissue(j).name  = 'adipose';
gg=0.9;
tissue(j).mua = MU(:,4);
tissue(j).mus = 0.93/(1-gg);
tissue(j).g   = gg;

j=7; %waardes vab magnetic resonance guided NIRT of breast 2004, Brooksby, geldig voor 785 nm!!
tissue(j).name  = 'figroglandular';
gg=0.9;
tissue(j).mua = 0.006;
tissue(j).mus = 1.12/(1-gg);
tissue(j).g   = gg;

j = 8;
tissue(j).name = 'dermis';
B = 0.002; 
S = 0.67;
W = 0.65;
F = 0;
M = 0;
musp500 = 42.4;
fray    = 0.62;
bmie    = 1.0;
gg      = 0.90;
musp = musp500*(fray*(nm/500).^-4 + (1-fray)*(nm/500).^-bmie);
X = [B*S B*(1-S) W F M]';
tissue(j).mua = MU*X;
tissue(j).mus = musp/(1-gg);
tissue(j).g   = gg;

disp(sprintf('---- tissueList ------ \tmua   \tmus  \tg  \tmusp'))
for i=1:length(tissue)
    disp(sprintf('%d\t%15s\t%0.4f\t%0.1f\t%0.3f\t%0.1f',...
        i,tissue(i).name, tissue(i).mua,tissue(i).mus,tissue(i).g,...
        tissue(i).mus*(1-tissue(i).g)))
end
disp(' ')

