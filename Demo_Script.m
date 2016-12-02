%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer Vision Final Project
% Matthew Walmer - mwalmer3
%
% Demo Script on how to use:
%   star_catalog_parser.m (or star_catalog_parser_full.m)
%   Sim_Star_Image.m
%   Star_Detector.m
%
% Other dependencies
%   hygdata_v3.mat
%   Star_Blob.m
%   UV_fraction_by_T.m
%   uf_root.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('parsed_star_db.mat', 'file')
    star_catalog_parser
    clear;
end
if ~exist('parsed_star_db_full.mat', 'file')
    star_catalog_parser_full
    clear;
end

%Simulate Star Image.
h_dim = 720;        % Horizontal image dimension
v_dim = 480;        % Vertical image dimension
h_fov = 20;         % Horizontal field of view in degrees
Sim_Star_Image('001', h_dim, v_dim, h_fov, true, true, false, false);

%%
[BV, C, UV] = Star_Detector('001_image_vis.jpg','001_image_uv.jpg',...
   '001_cp.mat', true, '001_truth.mat');
