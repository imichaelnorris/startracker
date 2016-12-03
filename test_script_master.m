%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer Vision Final Project
% Phillip Rivera - priver14
%
% test_script_master - master function to test the star_ID pyramid based
% star identification function.
%
% Internal Functions:
% true_star_ID - Obtains the real index of the test image using the
% intertial vector and camera rotation matrix used to generate the image.
%


function test_script_master
clc
close all


num_of_tests = 100; % Total Number of images tested for statistical analysis
score = zeros(num_of_tests,1); % Final percentage of good matches
max_attempts = 200;
load('parsed_star_catalog.mat')

for test_num = 1:num_of_tests
    
    
    rng(test_num); % Control Random Number Generation for monte carlo analysis
    
    close all
    disp(['Test Number: ' num2str(test_num)])
    
    %%% Simulate Star Image.
    Sim_Star_Image
    
    %%% Star Detection on Image
    [BV, C, UV] = Star_Detector('001_image_vis.jpg','001_image_uv.jpg',...
        '001_cp.mat', true, '001_truth.mat');
    
    %%% Determines True Star Indecies
    [~,star_ID_body] = true_star_ID(position,im_list_3,R,BV);
    
    %%% Call Star-ID
    
    
    % Star ID settings
    use_uv   = 1;   % Use UV information
    strng_uv = 1;   % Use Strong UV threshold
    sigma = 50e-6; % Image centroiding error
    star_ID_results = [];
    attempt  = 1;  % Pyramid Permutation Counter
    n        = size(BV,1); % number of identified stars
    
    while(isempty(star_ID_results) && attempt < max_attempts)
        
        disp(['Attempt: ' num2str(attempt)])
        index_selection = randsample(size(BV,1),4); 
        
        % Calling Star ID 
        star_ID_results = star_ID(BV(index_selection,:),UV(index_selection),sigma,use_uv,strng_uv); 
        match = ismember(star_ID_results,star_ID_body(index_selection));
        
        attempt = attempt + 1;
    end
    
    if(isempty(match))
        match = 0;
    end
    
    score(test_num) = sum(double(match))/4;
end

disp(['Score :' num2str(sum(score)/num_of_tests*100) '%'])

function [star_ID_inertial,star_ID_body] = true_star_ID(position,im_list_3,R,BV)

n = length(im_list_3);
body_vector_map = zeros(n,3);
star_ID_inertial = zeros(n,1);
for i = 1:n
    [~,id,~] = intersect(position,im_list_3(i,:),'rows');
    star_ID_inertial(i) = id;
    body_vector_map(i,:) = R * im_list_3(i,:)';
end

ind_inertial_image = zeros(size(BV,1),1);
for i = 1:size(BV,1)
    dist = zeros(n,1);
    for j = 1:n
        
        dist(j) = sqrt(sum((BV(i,:) - body_vector_map(j,:)).^2));
        
    end
    
    [~,ind_smallest_dist] = min(dist);
    
    
    ind_inertial_image(i) = ind_smallest_dist;
end


star_ID_body = star_ID_inertial(ind_inertial_image);


