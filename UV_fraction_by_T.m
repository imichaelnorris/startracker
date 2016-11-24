% This function estimates the ratio of a star's visible spectrum brightness
% to it's UV spectrum brightness by integrating the black body radiation
% formula. The variable is T, the star's surface temperature in Kelvin.
% This is estimated based on the star's spectral classification.
function uv_vis_frac = UV_fraction_by_T(T)
h = 6.626e-34;
c = 2.997e8;
kb = 1.38e-23;

lamda = 1:800;
lamda = lamda * 10^-9;

bb_rad = (2.*pi.*h.*c.^2) ./ ( lamda.^5 .* (exp( h.*c./(lamda.*kb.*T) )-1) );

%Integrate Radiation over Visible Spectrum (400-800)
int_vis = sum(bb_rad(400:800));

%Integrate Radiation over UV (10-400)
int_uv = sum(bb_rad(10:400));

uv_vis_frac = int_uv/int_vis;
end