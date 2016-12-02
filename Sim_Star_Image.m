%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer Vision Final Project
% Matthew Walmer - mwalmer3
%
% Star Image Simulator - Creates simulated star images for a near-earth
% camera with a random orientation.
%
% Inputs:
%   prefix      filename prefix for saved files
%   h_dim       horizonatal size of image (pixels)
%   v_dim       vertical size of image (pixels)
%   h_fov       horizonatl field of view (degrees)
%   fullDB      If true, full star database is used to generate images
%   addNoise    If true, gaussian noise is added to the images
%   showIm      If true, the final images will be displayed at the end
%   supSave     If true, file saving is suppressed
%   R_i_c       Rotation matrix from the inertial frame to the camera frame
%               If not given, a random rotation matrix is generated
%
% Outputs:
%  (PREFIX)_image_vis.jpg     Visible spectrum star image
%  (PREFIX)_image_uv.jpg      UV spectrum star image
%  (PREFIX)_cp.mat            File containing camera parameters fx,fy,ox,oy
%  (PREFIX)_truth.mat         File containing ground truth data about image:
%     im_list_ID              hygdata_v3 star ID for each star drawn in the
%                             image.
%     im_list                 Sub-pixel positions of stars in image in 2D pixel space
%     im_list_mag             Corresponding apparent magnitudes of stars in image
%     im_list_class           Corresponding classes of stars in image (using numerical coding)
%     im_list_3               3D points of included stars in inertial frame
%     im_list_3_cf            3D points in camera reference frame
%     R_c_i                   Rotation matrix that moves vectors from the
%                             camera reference frame to the inertial frame
%     fullDB                  Boolean value. True if full database was used to generate image.     
%     im_list_bright          Sub-pixel positions of just bright stars (magnitude <= Mvt)
%                             in image in 2D pixel space.
%     im_list_b_class         Corresponding classes of bright stars in image
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [I_vis, I_uv] = Sim_Star_Image(prefix, h_dim, v_dim, h_fov, fullDB, addNoise, showIm, supSave, R_i_c)

% default settings
if nargin < 5
    fullDB = true; % default to full database
end
if nargin < 6
    addNoise = true; % default to adding noise
end
if nargin < 7
    showIm = false; % default to not show images
end
if nargin < 8
    supSave = false; % default to saving files
end

%%%%%% OTHER PARAMS %%%%%%
Mvt = 5;            % Threshold for apparent magnitude of stars to include
                    % in im_list_bright if full database is used
Fx0 = 100;          % Relative brightness value - simulates exposure level
% the goal is to fine-tune Fx0 so that only stars with 
% apparent magnitude <= Mvt stand out signifigantly against the 
% background noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% STAR DATABASE %%%%%%
if fullDB
    load('parsed_star_db_full.mat');
else
    load('parsed_star_db.mat');
end
% A parsed database can be created using 'star_catalog_parser.m', which
% removes darker stars, or using 'star_catalog_parser_full.', which does
% not remove any stars.
%
% star database should have the following elements:
% -position     3D positions of stars on the unit sphere, as an Nx3 matrix
% -mag_Ap       Aparent magnitude of stars as an Nx1 matrix
% -spec_class   Spectral Classification of stars as an Nx1 matrix using the
%               following coding: [O B A F G K M ?] = [1 2 3 4 5 6 7 8],
%               where '?' means unclassified
% -star_id      ID of each star in the hygdata_v3 database
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute additional camera parameters
fy = h_dim / ( 2 * tan( h_fov * pi / 360 ));
fx = fy;
ox = h_dim/2;
oy = v_dim/2;

% Estimate ratio of UV Brightness to Visible Spectrum Brightness based on
% star temperature:
class_temp = [35000, 20000, 9000, 6500, 5500, 4500, 3000, 3000]; %O B A F G K M ?
class_UV_ratio = zeros(size(class_temp));
for i = 1:size(class_temp,2)
    class_UV_ratio(i) = UV_fraction_by_T(class_temp(i));
end

%%
%Generate Random Rotation Matrix
if nargin < 9
    % if no input rotation matrix is specified, a random one is generated
    tx = rand*2*pi;
    Rx = [1,0,0;0,cos(tx),-sin(tx);0,sin(tx),cos(tx)];
    ty = rand*2*pi;
    Ry = [cos(ty),0,-sin(ty);0,1,0;sin(ty),0,cos(ty)];
    tz = rand*2*pi;
    Rz = [cos(tz),-sin(tz),0;sin(tz),cos(tz),0;0,0,1];
    R_i_c = Rx * Ry * Rz;   %rotation from inertial frame to camera frame
end
R_c_i = R_i_c'; %rotation from camera frame to inertial frame

%%
%Generate Star Points in (u,v) pixel space:
num_stars = size(position,1);

%Store data about what stars appear in the image for ground truth file
num_in_frame = 0;
in_frame = zeros(num_stars, 2); % Sub-pixel positions of stars in image in 2D pixel space
in_frame_ID = zeros(num_stars, 1);

%Process all stars in database and check which appear in frame
for i = 1:num_stars
    p_w = position(i,:)';
    p_c = R_i_c * p_w;
    if p_c(3) > 0 % on correct side of image plane
        u = fx * p_c(1)/p_c(3) + ox;
        v = fy * p_c(2)/p_c(3) + oy;
        u_r = round(u);
        v_r = round(v);
        
        %Check if star appears in image. If so, add to im_list
        if u_r > 0 && u_r <= h_dim && v_r > 0 && v_r <= v_dim
            num_in_frame = num_in_frame + 1;
            in_frame(num_in_frame,:) = [u,v];
            in_frame_ID(num_in_frame) = i;
        end
    end
end
%truncate off zeros
in_frame = in_frame(1:num_in_frame,:);
in_frame_ID = in_frame_ID(1:num_in_frame,:);

%%
%Create Visible Spectrum and UV images
I_vis = zeros(v_dim,h_dim); % vis
I_uv = zeros(v_dim,h_dim); % uv

%add noise (optional)
if addNoise
    I_vis = imnoise(I_vis,'gaussian',0,0.005);
    I_uv = imnoise(I_uv,'gaussian',0,0.005);
end

%Store data about what stars are rendered in the image for the ground truth file
num_in_image = 0;
im_list_ID = zeros(num_in_frame, 1);       % hygdata_v3 star ID for each star drawn in the image.
im_list = zeros(num_in_frame, 2);          % Sub-pixel positions of stars in image in 2D pixel space
im_list_mag = zeros(num_in_frame, 1);      % Corresponding apparent magnitudes of stars in image
im_list_class = zeros(num_in_frame, 1);    % Corresponding classes of stars in image (using numerical coding)
im_list_3 = zeros(num_in_frame, 3);        % 3D points of included stars in inertial frame
im_list_3_cf = zeros(num_in_frame, 3);     % 3D points in camera reference frame

%If using full database, keep a seperate list of stars that appear in the
%image that have apparent magnitude less than Mvt (high brightness stars)
num_bright = 0;
im_list_bright = zeros(num_in_frame, 2);   % Sub-pixel positions of bright stars in image in 2D pixel space
im_list_b_class = zeros(num_in_frame, 1);   % Corresponding classes of bright stars in image

for i = 1:num_in_frame
    c = in_frame(i,1);
    r = in_frame(i,2);
    ID = in_frame_ID(i); %database ID of star

    Fx1 = Fx0 * 10^( mag_Ap(ID) / (-2.5) );
    B1 = Star_Blob(c,r,Fx1); % create star blob pattern for vis-image

    %Ratio of UV brightness to Visible-Spectrum brightness depends on
    %the star's class.
    Fx2 = Fx1 * class_UV_ratio( spec_class(ID) );
    B2 = Star_Blob(c,r,Fx2); % create star blob pattern for uv-image

    mat_size = size(B1,1);
    try % try adding blob patterns into images - this may fail if the star is close to an edge
        I_vis(round(floor(r)+1-mat_size/2):round(floor(r)+mat_size/2), round(floor(c)+1-mat_size/2):round(floor(c)+mat_size/2)) = ...
           I_vis(round(floor(r)+1-mat_size/2):round(floor(r)+mat_size/2), round(floor(c)+1-mat_size/2):round(floor(c)+mat_size/2)) + B1;
        I_uv(round(floor(r)+1-mat_size/2):round(floor(r)+mat_size/2), round(floor(c)+1-mat_size/2):round(floor(c)+mat_size/2)) = ...
           I_uv(round(floor(r)+1-mat_size/2):round(floor(r)+mat_size/2), round(floor(c)+1-mat_size/2):round(floor(c)+mat_size/2)) + B2;            

        %if the star is drawn sucessfully, add it's information to the
        %ground truth file. If not, an exception will be thrown and this
        %code will be skipped.
        num_in_image = num_in_image + 1;
        im_list_ID(num_in_image) = star_id(ID);
        im_list(num_in_image,:) = [c,r];
        im_list_mag(num_in_image) = mag_Ap(ID);
        im_list_class(num_in_image) = spec_class(ID);
        im_list_3(num_in_image,:) = position(ID,:);
        im_list_3_cf(num_in_image,:) = (R_i_c * position(ID,:)')';
        if mag_Ap(ID) <= Mvt
            %Bright Star
            num_bright = num_bright + 1;
            im_list_bright(num_bright,:) = [c,r];
            im_list_b_class(num_bright) = spec_class(ID);
        end
    catch
        %If a star is too close to the edge, it triggers an out of
        %bounds error. In this case we just skip the star
        continue
    end
end

%truncate off zeros
im_list_ID = im_list_ID(1:num_in_image,:);
im_list = im_list(1:num_in_image,:);
im_list_mag = im_list_mag(1:num_in_image,:);
im_list_class = im_list_class(1:num_in_image,:);
im_list_3 = im_list_3(1:num_in_image,:);
im_list_3_cf = im_list_3_cf(1:num_in_image,:);
im_list_bright = im_list_bright(1:num_bright,:);
im_list_b_class = im_list_b_class(1:num_bright,:);

%%
%save files (if not suppressed)
if ~supSave
    name = strcat(prefix, '_image_vis.jpg');
    imwrite(I_vis, name, 'jpg')
    name = strcat(prefix, '_image_uv.jpg');
    imwrite(I_uv, name, 'jpg')
    name = strcat(prefix, '_truth.mat');
    save(name,'im_list_ID','im_list','im_list_class','im_list_mag',...
        'im_list_3','im_list_3_cf','R_c_i','R_i_c','fullDB','im_list_bright',...
        'im_list_b_class');
    name = strcat(prefix, '_cp.mat');
    save(name,'fx','fy','ox','oy');
end

%%
% display simulated images (optional)
if showIm
    figure;
    imshow(I_vis)
    figure;
    imshow(I_uv)
    % display vis-spectrum image with overlay
    figure;
    imshow(I_vis)
    hold on
    plot(im_list(:,1),im_list(:,2),'ro')
    % display vis-spectrum image with overlay for only bright stars
    if fullDB
        figure;
        imshow(I_vis)
        hold on
        plot(im_list_bright(:,1),im_list_bright(:,2),'ro')
    end
end

end
