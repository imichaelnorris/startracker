clear all
clc

load('parsed_star_catalog.mat')

% Body vector pyramid with angles smaller than 20deg
% star_1 = 1; uv_1 = ismember(star_1,ind_uv);
% star_2 = 3; uv_2 = ismember(star_2,ind_uv); 
% star_3 = 4; uv_3 = ismember(star_3,ind_uv);
% star_4 = 5; uv_4 = ismember(star_4,ind_uv);

star_1 = 265; uv_1 = ismember(star_1,ind_uv);
star_2 = 290; uv_2 = ismember(star_2,ind_uv); 
star_3 = 353; uv_3 = ismember(star_3,ind_uv);
star_4 = 361; uv_4 = ismember(star_4,ind_uv);

% Star ID settings
use_uv   = 0;   % Use UV information
strng_uv = 1;   % Use Strong UV threshold
sigma = 100e-6; % Image centroiding error

body_vectors = [position(star_1,:)
                position(star_2,:)
                position(star_3,:)
                position(star_4,:)];

uv_content = [uv_1; uv_2; uv_3; uv_4];

indices = star_ID(body_vectors,uv_content,sigma,use_uv,strng_uv)