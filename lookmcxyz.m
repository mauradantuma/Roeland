function lookmcxyz(myname, nm);

% lookmcxyz.m
%   Looks at myname_F.bin, created by mcxyz.c
%   where myname is the name of the run: myname_T.bin, myname_H.mci
%
% lookmcxyz_alc2.m, July 23
%       Absorption array is OFF.
%
%   Simulates angiolight catheter within 4-mm-dia. vessel
% with the vessel wall at a particular musp value.
%   For this run, musp = 200 cm^-1.
%
% Reads 8 mcxyz runs, for 8 catheter positions, zs = 0.2 to 0.6 cm.
% For each run:
%   alc#_H.mci --> header:
%       timin,Nx,Ny,Nz,dy,dx,dz,xs,ys,zx,Nt,muav(),musv(),gv()
%   alc#_F.bin --> F(y,x,z) = relative fluence rate [1/cm^2]
%   alc#_T.bin --> T(y,x,z) = tissue types
%
% Displays
%   Tzx = end-view of tissue structure
%   Fzx = end-view of vessel @ source, ys = -0.15 cm
%   Fzy = side-view along length of vessel
%
% Saves
%   Fzy_data4.mat = Fzy y z zzs Fdet
%       Fzy(400,400,8) = 8 z,y images
%       Fdet(8,1) = signal [1/cm^2] @ detector fiber
%
% home
% % clear
% format compact
% commandwindow

cc = 'rbgm'; % color

disp(sprintf('------ mcxyz %s -------',myname))

% Load header file
filename = sprintf('%s_H.mci',myname);
disp(['Loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
time_min = A(1);
reqphotons = A(2);
Nx = A(3);
Ny = A(4);
Nz = A(5);
dx = A(6);
dy = A(7);
dz = A(8);
mcflag = A(9);
launchflag = A(10);
xs = A(12);
ys = A(13);
zs = A(14);
xfocus = A(15);
yfocus = A(16);
zfocus = A(17);
ux0 = A(18);
uy0 = A(19);
uz0 = A(20);
radius = A(21);
waist = A(22);
focal_length = A(23);
angle_x=A(24);
angle_y=A(25);
Nt = A(26);
j = 26;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
end

reportHmci(myname);

%% Load Fluence rate F(y,x,z)
filename = sprintf('%s_F.bin',myname);
disp(['Loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)

%load deposited energy matrix (A)
filename = sprintf('%s_A.bin',myname);
disp(['Loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
A = reshape(Data,Ny,Nx,Nz); % F(y,x,z)

save(strcat([myname, '_A.mat']),'A','-v7.3');
save(strcat([myname, '_F.mat']),'F','-v7.3');

% % Load tissue structure in voxels, T(y,x,z) 
% filename = sprintf('%s_T.bin',myname);
% disp(['loading ' filename])
% tic
%     fid = fopen(filename, 'rb');
%     [Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
%     fclose(fid);
% toc
% T = reshape(Data,Ny,Nx,Nz); % T(y,x,z)

clear Data

%%
x = ([1:Nx]-Nx/2-1/2)*dx;
y = ([1:Ny]-Ny/2-1/2)*dx;
z = ([1:Nz]-1/2)*dz;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

%% Look at structure, Tzx
% Tzx = reshape(T(round(Ny/2),:,:),Nx,Nz)';
% tissue = makeTissueList(nm);
% Nt = length(tissue);
% 
% figure(1);clf
% imagesc(x(ux),z(uz),Tzx(uz,ux),[1 Nt])
% hold on
% cmap = makecmap(Nt);
% colormap(cmap)
% colorbar
% set(gca,'fontsize',18)
% set(colorbar,'fontsize',1)
% xlabel('x [cm]')
% ylabel('z [cm]')
% title('Tissue types','fontweight','normal')
% for i=1:Nt
%     yy = zmin + (Nt-i)/(Nt-1)*zdiff;
%     text(xmin + xdiff*1.1,yy, sprintf('%d %s',i,tissue(i).name),'fontsize',9)
% end
% 
% % draw launch, trajectory lijnen plotten
% 
% for j=-1:0.5:1;
%     plot([xs+j*radius xfocus],[zs zfocus],'r-'); hold on
% end
% 
% name = sprintf('%s_tissue.jpg',myname);
% print('-djpeg','-r300',name)


% Look at Fluence Fzx @ launch point
% fluence
Fzx = reshape(F(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source

figure(2);clf; subplot(1,2,1);
imagesc(x,z,log10(Fzx),[-3 3])
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [J/cm^2/J_{delivered}] ','fontweight','normal')
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])

name = sprintf('%s_Fzx.jpg',myname);
print('-djpeg','-r300',name)

%deposited energy
Azx = reshape(A(round(Ny/2),:,:),Nx,Nz)'; % in z,x plane through source
figure(2); subplot(1,2,2); 
imagesc(x,z,log10(Azx),[-3 3])
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'log_{10}(E)','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('x [cm]')
ylabel('z [cm]')
title('Deposited energy [J/cm^3/J_{delivered}] ','fontweight','normal')
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])

name = sprintf('%s_Azx.jpg',myname);
print('-djpeg','-r300',name)

% Look at Fluence Fzy @ launch point
% fluence
Fzy = reshape(F(:,round(Ny/2),:),Ny,Nz)'; % in z,x plane through source

figure(3);clf; subplot(1,2,1);
imagesc(x,z,log10(Fzy),[-3 3])
hold on
text(max(y)*0.9,min(z)-0.04*max(z),'log_{10}( \phi )','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [J/cm^2/J_{delivered}] ','fontweight','normal')
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])

name = sprintf('%s_Fzy.jpg',myname);
print('-djpeg','-r300',name)

%deposited energy
Azy = reshape(A(:,round(Ny/2),:),Ny,Nz)'; % in z,x plane through source
figure(3); subplot(1,2,2); 
imagesc(y,z,log10(Azy),[-3 3])
hold on
text(max(y)*0.9,min(z)-0.04*max(z),'log_{10}(E)','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('x [cm]')
ylabel('z [cm]')
title('Deposited energy [J/cm^3/J_{delivered}] ','fontweight','normal')
colormap(makec2f)
axis equal image
%axis([min(x) max(x) min(z) max(z)])

name = sprintf('%s_Azy.jpg',myname);
print('-djpeg','-r300',name)

% look Fxy

%fluence
Fzy = reshape(F(:,:,round(Nz/4)),Ny,Nx)';
figure(4);clf; subplot(1,2,1);
imagesc(y,x,log10(Fzy),[-1 1]*3)
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'log_{10}(\phi)','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('y [cm]')
ylabel('x [cm]')
title('Fluence \phi [J/cm^2/J_{delivered}] ','fontweight','normal')
colormap(makec2f)
axis equal image

name = sprintf('%s_Fzy.jpg',myname);
print('-djpeg','-r300',name)


%deposited energy
Azy = reshape(A(:,:,round(Nz/4)),Ny,Nx)';
figure(4); subplot(1,2,2);
imagesc(y,x,log10(Azy),[-1 1]*3)
hold on
text(max(x)*0.9,min(z)-0.04*max(z),'log_{10}(E)','fontsize',18)
colorbar
set(gca,'fontsize',18)
xlabel('Y [cm]')
ylabel('X [cm]')
title('Deposited energy [J/cm^3/J_{delivered}] ','fontweight','normal')
colormap(makec2f)
axis equal image

name = sprintf('%s_Azy.jpg',myname);
print('-djpeg','-r300',name)

drawnow

disp('done')



end
