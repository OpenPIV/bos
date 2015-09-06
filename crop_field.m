function [Displ]=crop_field(Displacement_POisson,Lx,Lz)

Displ=Displacement_POisson;
Minimum=min(min(Displ.y));
Displ.y=Displ.y-abs(Minimum);
Magnitude=sqrt(Displ.u.^2+Displ.v.^2);
[a,b]=size(Displ.x);

% Crop The figure because the remapping algorithm
% creates the external frames (lack in the data)
Dxx=Lx/a;
Dyy=Lz/b;
nx_pixels_crop=700;  %250
ny_pixels_crop=350;  %200
Lx=Lx-nx_pixels_crop;
Lz=Lz-ny_pixels_crop;
Dx_pixels=round(nx_pixels_crop/Dxx);   % number of pixels to crop in th eright and left side of the image
Dy_pixels=round(ny_pixels_crop/Dyy);  % number of pixels to crop in th eright and left side of the image
Displ.x=Displ.x(Dx_pixels:end-Dx_pixels,Dy_pixels:end-Dy_pixels);
Displ.y=Displ.y(Dx_pixels:end-Dx_pixels,Dy_pixels:end-Dy_pixels);
Displ.u=(Displ.u(Dx_pixels:end-Dx_pixels,Dy_pixels:end-Dy_pixels)); % check the image imshow(Displ{1,1}.u)
Displ.v=(Displ.v(Dx_pixels:end-Dx_pixels,Dy_pixels:end-Dy_pixels));
Magnitude_crop=sqrt(Displ.u.^2+Displ.v.^2);

Displ;