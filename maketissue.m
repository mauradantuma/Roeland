% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       myname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.c, [mua, mus, g], for each tissue type specified
%   in myname_T.bin. This listing is saved in
%       myname_H.mci = the input file for use by mcxyz.c.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use,
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%

clear all
close all

%format compact
clc
savefirsttime=1;

for fiber_nr=[2];

%loading phantom and setting the axis right
load(strcat(['space_fiber',num2str(fiber_nr),'.mat'])); %Txyz

T=permute(T,[3 1 2]); %Txyz --> Tzxy
T=imrotate(T,180); %spiegelen in zx vlak
T=permute(T,[2 3 1]); %Tzxy --> Txyz, maar dan nu gespiegeld over de z as

[r,k,d]=size(T);
max_dim=max([r,k,d]);


%% Set angle of incidence in degrees

% only used if mcflag == 0 or 1 (uniform or Gaussian beam, not isotropic pt.)
radius = 2.46;   % 1/e radius of beam 
waist = 0.02;    %  wordt niet gebruikt, ik ga ervanuit dat de focus 1 punt is. 

angle_x=[5]; %between 0 and 180 degrees;
angle_y=[90]; %between 0 and 180 degrees;
focal_length=(radius/tand(21.2)); %collimated beam --> focal length = inf;
% 
% computingmodel=2; %1 for loaded image, 2 for generated model

%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;            % 1 = save myname_T.bin, myname_H.mci 0 = don't save. Just check the program.
myname      = strcat(['fiber_', num2str(fiber_nr)]); % name for files: myname_T.bin, myname_H.mci
time_min    = 0.1;      	% time duration of the simulation [min], not used
reqphotons  = 1e6;%8.23e4;       % Requested Photons, overrules time_min
nm          = 755;   	    % desired wavelength of simulation
% Nbins       = 500;    	    % # of bins in each dimension of cube
binsize     = 0.02;       %[cm]	% size of each bin, eg. [cm] or [mm]

% Set Monte Carlo launch flags
mcflag      = 0;            % launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt.
launchflag  = 0;            % 0 = let mcxyz.c calculate launch trajectory % 1 = manually set launch vector.
boundaryflag = 1;       %    0 = no boundaries, 1 = escape at boundaries 2 = escape at surface only. No x, y, bottom z boundaries

%% Set position of Source and Focus
% Sets position of source center
xs=6;       %x source
ys=0;       %y source
zs=4.6874;      %z source

if isinf(focal_length); focal_length=1e12; end % als collimated beam, set focal length to 1e12

% Set position of focus, so mcxyz can calculate launch trajectory
xfocus=xs+focal_length*sind((90-angle_x));
yfocus=ys+focal_length*sind((90-angle_y));
zfocus=(zs+focal_length)-(focal_length*sind(90-angle_x)*tand((90-angle_x)/2))-(focal_length*sind(90-angle_y)*tand((90-angle_y)/2));

%unit vector for direction of beam, only used if launchflag=1;
ux0         = cosd(angle_x);      % trajectory projected onto x axis
uy0         = cosd(angle_y);      % trajectory projected onto y axis
uz0         = sqrt(1 - ux0^2 - uy0^2); % such that ux^2 + uy^2 + uz^2 = 1

%% Prepare Monte Carlo
% Create tissue properties
tissue = makeTissueList_fat(nm);    % make tissue properties
Nt = length(tissue);            % number of tissues

%making vectors containing mua mus and g;
for i=1:Nt;
    muav(i)  = tissue(i).mua;
    musv(i)  = tissue(i).mus;
    gv(i)    = tissue(i).g;
end

% Specify Monte Carlo parameters
Nx = max_dim;
Ny = max_dim;
Nz = max_dim;
dx = binsize;
dy = binsize;
dz = binsize;
x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);


%%
if SAVEON
    if savefirsttime==1;
        save(sprintf('tissuematrix_%d.mat',size(T,2)),'T','-v7.3')
        %save tissueproperties.mat tissue
%         save(sprintf('absorpmap_%d.mat',size(absorpmap,2)),'absorpmap')
        savefirsttime=0;
    end
    
    
    %tic;
    %     convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Nx*Ny*Nz,1));
    
    % WRITE FILES
    %     Write myname_H.mci file, this is used by mcxyz.c
    %     which contains the Monte Carlo simulation parameters
    %     and specifies the tissue optical properties for each tissue type.

    disp(sprintf('Creating %s ',myname))
    filename = sprintf('%s_H.mci',myname);
    fid = fopen(filename,'w');
    %run parameters
    fprintf(fid,'%0.2f\n',time_min);
    fprintf(fid,'%d\n'   ,reqphotons);
    fprintf(fid,'%d\n'   ,Nx);
    fprintf(fid,'%d\n'   ,Ny);
    fprintf(fid,'%d\n'   ,Nz);
    fprintf(fid,'%0.4f\n',dx);
    fprintf(fid,'%0.4f\n',dy);
    fprintf(fid,'%0.4f\n',dz);
    %launch parameters
    fprintf(fid,'%d\n'   ,mcflag);
    fprintf(fid,'%d\n'   ,launchflag);
    fprintf(fid,'%d\n'   ,boundaryflag);
    fprintf(fid,'%0.4f\n',xs); %12
    fprintf(fid,'%0.4f\n',ys);
    fprintf(fid,'%0.4f\n',zs);
    fprintf(fid,'%0.4f\n',xfocus); %15
    fprintf(fid,'%0.4f\n',yfocus);
    fprintf(fid,'%0.4f\n',zfocus);
    fprintf(fid,'%0.4f\n',ux0); % if manually setting ux,uy,uz
    fprintf(fid,'%0.4f\n',uy0);
    fprintf(fid,'%0.4f\n',uz0);
    fprintf(fid,'%0.4f\n',radius);
    fprintf(fid,'%0.4f\n',waist);
    fprintf(fid,'%0.4f\n',focal_length);
    fprintf(fid,'%0.4f\n',angle_x);
    fprintf(fid,'%0.4f\n',angle_y);
    %tissue optical properties
    fprintf(fid,'%d\n',Nt);
    for i=1:Nt
        fprintf(fid,'%0.4f\n',muav(i));
        fprintf(fid,'%0.4f\n',musv(i));
        fprintf(fid,'%0.4f\n',gv(i));
    end
    fclose(fid);
    
    % write myname_T.bin file
    filename = sprintf('%s_T.bin',myname);
    disp(['Writing ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);
    
    %toc
end % SAVEON

%% Look at structure of Tzx at iy=Ny/2
% Txzy = permute(T,[1 3 2]);   % Txyz --> Txzy
% Tzx  = Txzy(:,:,round(Ny/2))'; % Tzx slice op de helft van y
% clear Txzy;

% %% tissue matrix laten zien op y=Ny/2
% figure(2); clf
% fsz = 18;  % font size
% imagesc(x,z,Tzx,[1 Nt])
% hold on
% set(gca,'fontsize',fsz)
% xlabel('x [cm]')
% ylabel('z [cm]')
% colorbar
% cmap = makecmap(Nt);
% colormap(cmap)
% set(colorbar,'fontsize',1)
% % label colorbar
% zdiff = zmax-zmin;
% 
% %%%
% for i=1:Nt
%     yy = (Nt-i)/(Nt-1)*Nz*dz;
%     text(Nx*dx*1.2,yy, tissue(i).name,'fontsize',12)
% end
% 
% text(xmax*0.9,zmin - zdiff*0.06, 'Tissue types','fontsize',18)
% axis equal image
% axis([xmin xmax zmin zmax])
% 
% 
% %% draw launch, trajectory lijnen plotten
% 
% for j=-1:0.5:1;
%     plot([xs+j*radius xfocus],[zs zfocus],'k-'); hold on
% end

%% mcxyz.c script runnen en resultaten laten zien
disp('Launching Monte Carlo')
tic
system('compile.bat');
system(strcat(['carlo' num2str(fiber_nr) '.bat']));
toc

disp('Launching lookmcyz.mat')
lookmcxyz(myname,nm);
end


