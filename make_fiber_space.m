%% dit script heb ik gebruikt om de tissue matrices te maken voor de verschillende fibers
% voor elke top fiber heb ik de fiber positie steeds constant gehouden,
% maar het fantoom gedraaid onder een hoek van 27.69 graden (=hoek tussen
% fibers)
clear all; close all;

load('space_fiber1.mat'); %originele tissue matrix
% T=space; clear space; 
A=T; clear T;
% space=permute(A,[3,2,1]); %permuten
% space=A(:,:,400);
% clear A;
i=1;
for fiber_nr=[3 5 6 8 9 11 12 14];
    
    fiber_angle=(fiber_nr-2)*27.69;
    
    T=imrotate(A,fiber_angle); %roteren in xy vlak
    T(T==0)=3; %extra toegevoegde pixels gelijk stellen aan 3 = water
    [r,k,d]=size(T);
    center_x=round(r/2); 
    xmin = center_x-409; %matrix weer 719x719x719 maken
    xmax = center_x+409;
    center_y = round(k/2);
    ymin = center_y-409;
    ymax = center_y+409;
    T=T(xmin:xmax,ymin:ymax,:);
    
    subplot(3,3,i); imshow(T(:,:,400),[]);
    i=i+1;
%     T=permute(T,[3,2,1]); %terug permuten, zodat de rest van het script in de goede dimensies werkt
%     
%     T(T==3)=0;
%     T(T~=0)=1;
    
%     B=T; clear T;
    filename=strcat(['space_fiber', num2str(fiber_nr),'.mat']);
    save(filename,'T','-v7.3')
end