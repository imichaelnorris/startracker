%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer Vision Final Project
% Matthew Walmer - mwalmer3
%
% Star Blob - Creates a gaussian pattern to represent a blurred star. The
% center position of the gaussian need not be integer valued.
%
% Inputs
%   x_star      The x-position of the star. Only the decimal part matters
%   y_star      The y-position of the star. Only the decimal part matters
%   Fx          Brightness value for the star blob. Scales the gaussian
%               Note that in matlab's double images, anything above 1 is
%               just white.
%
% Output
%   B           A matrix storing the star blob pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = Star_Blob(x_star,y_star,Fx)
sig = 1;
mat_size = ceil(2*(3*sig + 1));

x_c = (mat_size/2) + x_star - floor(x_star);
y_c = (mat_size/2) + y_star - floor(y_star);

x = 1:1:mat_size;
y = x;
[X,Y] = meshgrid(x,y);

B = Fx / (2*pi*sig^2) * exp( -((X-x_c).^2 + (Y-y_c).^2)/(2*sig^2) );
end