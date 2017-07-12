%% dit script heb ik gebruikt om de fluence maps van de individuele fibers 
% te combineren tot 1 totale fluence map. Er zit hier nog een bug in, iets
% met de interpolatie tijdens het imrotaten.. als je kijkt naar de totale fluence 
% map die ik heb opgeslagen (M_fluence of M_absorb) zie je dat het soms aan
% de randen van bepaalde structuren niet helemaal goed gaat. Dus daar moet
% nog even naar gekeken worden. 

clear all; close all; clc;
% M_absorb=[];
M_fluence=[];

for fiber_nr=[2 3 5 6 8 9 11 12 14];
display('loading matrix')
tic
% load(strcat(['fiber_' num2str(fiber_nr) '_A.mat']));
load(strcat(['fiber_' num2str(fiber_nr) '_F.mat']));
toc
display('loading tissue matrix')
tic
load(strcat(['space_fiber' num2str(fiber_nr) '.mat']));
toc

if fiber_nr == 1 || fiber_nr==2
    fiber_angle=0;
else
   fiber_angle=1*(fiber_nr-2)*27.69;
end

%make binary mask
T(T==3)=0;
T(T~=0)=1;
T=permute(T,[3 1 2]); %Txyz --> Tzxy
T=imrotate(T,180); %spiegelen in zx vlak
T=permute(T,[2 3 1]); %Tzxy --> Txyz, maar dan nu gespiegeld over de z as

F=F.*T;
clear T;

display('making sparse matrix')
for i=1:819;
    F(:,:,i)=sparse(F(:,:,i));
end

display(strcat(['Rotating matrix of fiber' num2str(fiber_nr)]))
F=imrotate(F,fiber_angle);

[r,k,d]=size(F);
center_x=round(r/2);
xmin = center_x-400; %matrix weer 719x719x719 maken
xmax = center_x+400;
center_y = round(k/2);
ymin = center_y-400;
ymax = center_y+400;
F=F(xmin:xmax,ymin:ymax,:);

% multiplication factor - fraction of light coming out of that fiber
if fiber_nr==1;
    factor=2/3;
else 
    factor=1/27;
end

if fiber_nr==2;
    M_absorb=zeros(size(F));
end

% M_absorb=M_absorb+factor.*A;

M_fluence=M_fluence+factor.*F;

clear F;
end

%saving final matrices
% save('M_absorb.mat','M_absorb','-v7.3');
save('M_fluence.mat','M_fluence','-v7.3');


% i=1;
% for fiber_nr=[3 5 6 8 9 11 12 14];
%    
% load(strcat(['fiber_' num2str(fiber_nr) '_A.mat']));
% figure(1); subplot(3,3,i);
% imshow(log10(A(:,:,300)),[]);
% i=i+1;
% end
% 
% i=1;
% for fiber_nr=[3 5 6 8 9 11 12 14];
%    
% load(strcat(['space_fiber' num2str(fiber_nr) '.mat']));
% T=flipud(fliplr(T));
% T=permute(T,[3,2,1]);
% figure(2); subplot(3,3,i);
% imshow(log10(T(:,:,300)),[]);
% i=i+1;
% end

