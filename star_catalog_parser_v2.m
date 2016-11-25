clear all
close all
clc

load('hygdata_v3.mat')

Mv  = 5; % Apparent visual Magnitude threshold as in - A Brightness-Referenced Star Identification Algorithm for APS Star Trackers
ind = find(mag<= Mv);

star_id   = id(ind);        % Mapping to original table
position  = [x(ind) y(ind) z(ind)]./sqrt(x(ind).^2 + y(ind).^2 + z(ind).^2); %location in unit sphere
mag_Ap    = mag(ind);       % Apparent magnitude M_v
type      = {};
type      = spect(ind);     % Star Type
ind_table = 1:length(type); % Star Catalog (SC) table index

star_class = blanks(length(type));
for i = 1:length(type)
    if(~isempty(type{i}))
        star_class(i) = char(type{i}(1));
    else
        star_class(i) = [];
    end
end

% Selecting stars with high UV content
index_type_F = strfind(star_class, 'F')';
index_type_A = strfind(star_class, 'A')';
index_type_B = strfind(star_class, 'B')';
index_type_O = strfind(star_class, 'O')';


ind_uv          = sort([index_type_F;index_type_A;index_type_B;index_type_O]); %SC index of UV stars
ind_vis         = ind_table;
ind_vis(ind_uv) = []; % SC index of only visible stars

save('parsed_star_catalog.mat','position','mag_Ap','star_id','ind_table','ind_uv','ind_vis');

figure;
subplot(2,2,[1 2])
plot3(position(ind_vis,1),position(ind_vis,2),position(ind_vis,3),'r.'); grid on; hold on; axis equal
title('Star Catalog M_v < 5')
plot3(position(ind_uv,1),position(ind_uv,2),position(ind_uv,3),'b.');
subplot(2,2,3)
plot3(position(ind_vis,1),position(ind_vis,2),position(ind_vis,3),'r.'); grid on; axis equal
title('Only Visible Spectrum')
subplot(2,2,4)
plot3(position(ind_uv,1),position(ind_uv,2),position(ind_uv,3),'b.'); grid on; axis equal
title('Visible and UV Spectrum')
set(findall(gcf,'type','text'),'fontsize',18)

