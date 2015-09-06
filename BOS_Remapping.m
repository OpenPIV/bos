function new_image = BOS_Remapping(Displ,im1)
% IMAGE_REMAPPED = BOS_REMAPPING(CALIBRATION,MAGE)
% remaps the second image according to the values in CALIBRATION
%

% remapping has two steps:
% 1. interpolate the dx,dy fields to each pixel
% 2. warp the image according to the interpolated dx,dy
[m,n] = size(im1); % note that m is vertical, n is horizontal
[X,Y] = meshgrid(0:n-1,0:m-1); 


vi = interp2(Displ.x',Displ.y',Displ.v',X,Y,'linear');
ui = interp2(Displ.x',Displ.y',Displ.u',X,Y,'linear');

% 
% xi=linspace(1,n,n);
% yi=linspace(1,m,m);
% % debugging, to be removed
% figure, quiver(Displ.x',Displ.y',Displ.u',Displ.v','AutoScale','off');
% hold on
% quiver(xi(1:33:end,1:33:end),yi(1:33:end,1:33:end),ui(1:33:end,1:33:end),vi(1:33:end,1:33:end),'r','AutoScale','off');
% axis equal


% Replace the NaN using zeros
ui(isnan(ui)) = 0;
vi(isnan(vi)) = 0;

% Rotate the image and then apply the remapping (Lilly's change)
im1=flipud(im1);

% rotate apply the imwarp and then we rotate back the image. Like this Im
% sure the reference is in the same position.
%new_image = (imwarp(im1,ui,vi,true)); 

new_image = flipud(imwarp(im1,ui,vi,true)); 
