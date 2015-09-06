clc
clear all
close all

% Main script for the BOS application

%% LOAD THE IMAGES
im1=im2double(imread('Data/Air_ref.tif'));
im2=im2double(imread('Data/Water_ref.tif'));
im3=im2double(imread('Data/4layers.tif'));


%% LOAD THE PARAMETERS FILE
[Mconversion,Const,Lx,Lz,val_up,val_down,nx_pixel,ny_pixel,overlap_x,overlap_y]=Parameters();

%% CREATE THE CALIBRATION FIELD: Correlation air-water
[Calibration]=BOS_correlation_OpenPIV(im1,im2,nx_pixel,ny_pixel,overlap_x);

%%
figure
%imshow(im1)
hold on
quiver(Calibration.x,Calibration.y,Calibration.u,Calibration.v,'AutoScale','off');
axis equal

Magn_cal=sqrt(Calibration.u.^2+Calibration.v.^2);
figure
contour(Calibration.x,Calibration.y,Magn_cal,50)
axis equal
colorbar

%% APPLY THE REMAPPING 

[im3_remapped] = BOS_Remapping_Alex(Calibration,im3);   % Interp2 works if we fix an overlap of 0.5. Otherwise Matlab suggests to use the TriScatteredInterp 
imwrite(im3_remapped,'Remapped_4layers.tif')
%figure;imshowpair(im3,im3_remapped,'diff');

%% CORRELATION REFERENCE-REMAPPED

nx_pixel=32;
ny_pixel=32;
overlap_x=0.25;
[Displacement_POisson]=BOS_correlation_OpenPIV(im1,im3_remapped,nx_pixel,ny_pixel,overlap_x);

% Check the displacement corrected

figure
quiver(Displacement_POisson.x,Displacement_POisson.y,Displacement_POisson.u,Displacement_POisson.v,5);
axis equal


%save('Displacement_POisson1.mat','Displacement_POisson')
%% For comparison between corrected and not-corrected case

% Small displacement im2-im3. The interrogation area An have to be reduced
[Displ_notcorr]=BOS_correlation_OpenPIV(im2,im3,nx_pixel,ny_pixel,overlap_x);

% Check the displacement not-corrected
figure
quiver(Displ_notcorr.x,Displ_notcorr.y,Displ_notcorr.u,Displ_notcorr.v,5);
axis equal


%% POISSON INTEGRATION

%load('Displacement_POisson.mat');

[n2, xc, zc]=BOS_PoissonSolver(Displacement_POisson,Const,Lx,Lz);
[n2_nc, x_nc, z_nc]=BOS_PoissonSolver(Displ_notcorr,Const,Lx,Lz);

%% Gladstone-Dale conversion
[Dens,Dens_av]=Gladstone_Dale(n2,xc,zc);
[Dens_2,Dens_av_2]=Gladstone_Dale(n2_nc, x_nc, z_nc);

%% GRAPHICAL OUTPUT 

% ++++++++++++ Comparison Magnitudo +++++++++++++++++++++++++++++++++++++++

Magnitudo=sqrt(Displacement_POisson.u.^2+Displacement_POisson.v.^2);
Magnitudo_nc=sqrt(Displ_notcorr.u.^2+Displ_notcorr.v.^2);

figure
subplot(121)
contour(Displacement_POisson.x*Mconversion,Displacement_POisson.y*Mconversion,Magnitudo,20)
axis equal
h=colorbar;
xlabel('x [cm]')
ylabel('y [cm]')
axis equal
xlim([min(min(Displacement_POisson.x*Mconversion)), max(max(Displacement_POisson.x*Mconversion))])
ylim([min(min(Displacement_POisson.y*Mconversion)), max(max(Displacement_POisson.y*Mconversion))])
caxis([0 5])
title('Corrected')
set(gca,'Ydir','reverse')

subplot(122)
contour(Displ_notcorr.x*Mconversion,Displ_notcorr.y*Mconversion,Magnitudo_nc,20)
axis equal
h=colorbar;
xlabel('x [cm]')
ylabel('y [cm]')
axis equal
xlim([min(min(Displ_notcorr.x*Mconversion)), max(max(Displ_notcorr.x*Mconversion))])
ylim([min(min(Displ_notcorr.y*Mconversion)), max(max(Displ_notcorr.y*Mconversion))])
caxis([0 5])
title('Not-corrected')
set(gca,'Ydir','reverse')

% +++++++++++ Results: Corrected Magnitude, Density field, Density profiles +++++++++++++++++++++++++++++++++++++++

figure
subplot('position',[0.08 0.35 0.3 0.4]);  
contour(Displacement_POisson.x*Mconversion,Displacement_POisson.y*Mconversion,Magnitudo,20)
axis equal
h=colorbar;
xlabel('x [cm]')
ylabel('y [cm]')
xlim([min(min(Displacement_POisson.x*Mconversion)), max(max(Displacement_POisson.x*Mconversion))])
ylim([min(min(Displacement_POisson.y*Mconversion)), max(max(Displacement_POisson.y*Mconversion))])
% caxis([0 10])
title('Corrected')
set(gca,'Ydir','reverse')

subplot('position',[0.42 0.35 0.3 0.4]);  
pcolor(Dens.x*Mconversion, Dens.z*Mconversion, Dens.f')
shading flat
axis equal tight 
xlabel('x [cm]')
title('Density Corrected')
colorbar;
set(gca,'Ydir','reverse')

subplot('position',[0.77 0.4 0.20 0.3]); 
hold on
plot(Dens_av,Dens.z*Mconversion,'b','linewidth',1.5)
plot(Dens_av_2,Dens_2.z*Mconversion,'b--','linewidth',1.5)
xlabel('\rho [g/mL]')
xlim([0.99 1.3])
ylim([0 28])
legend('Corrected','Not-Corrected')
set(gca,'Ydir','reverse')

% ++++++++++++ Comparison POisson solutions +++++++++++++++++++++++++++++++


figure
subplot(121)
surf(xc,zc,n2')
xlabel('x [px]')
ylabel('y [px]')
zlabel('n')
title('Corrected')

subplot(122)
surf(x_nc,z_nc,n2_nc')
xlabel('x [px]')
ylabel('y [px]')
zlabel('n')
title('Not-corrected')

