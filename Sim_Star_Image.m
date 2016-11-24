%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer Vision Final Project
% Matthew Walmer - mwalmer3
%
% Star Image Simulator - Creates simulated star images for a near-earth
% camera with a random orientation.
%
% Outputs:
%  (PREFIX)_image_vis.jpg     Visible spectrum star image
%  (PREFIX)_image_uv.jpg      UV spectrum star image
%  (PREFIX)_cp.mat            File containing camera parameters fx,fy,ox,oy
%  (PREFIX)_truth.mat         File containing ground truth data about image:
%     im_list                 Sub-pixel positions of stars in image in 2D pixel space
%     im_list_mag             Corresponding apparent magnitudes of stars in image
%     im_list_class           Corresponding classes of stars in image (using numerical coding)
%     im_list_3               3D points of included stars in inertial frame
%     im_list_3_cf            3D points in camera reference frame
%     R_c_i                   Rotation matrix that moves vectors from the
%                             camera reference frame to the inertial frame
%  Where (PREFIX) is a customizable variable
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% FILE NAME PREFIX %%%%%%
% varies names of output files
prefix = '001';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% PARAMS %%%%%%
u_max = 720;        % Horizontal image dimension
v_max = 480;        % Vertical image dimension
h_fov = 30;         % Horizontal field of view in degrees
Fx0 = 200;          % Relative brightness value - simulates exposure level
showIm = false;     % Display images at the end
addNoise = true;   % Add gaussian noise to simulated images
%%%%%%%%%%%%%%%%%%%%


%%%%%% STAR DATABASE %%%%%%
% Specify parsed star database file name:

star_db_file = 'parsed_star_db.mat';
%star_db_file = 'parsed_star_db_full.mat';

% A parsed database can be created using 'star_catalog_parser.m', which
% removes darker stars, or using 'star_catalog_parser_full.', which does
% not remove any stars.
%
% star database should have the following elements:
% -position     3D positions of stars on the unit sphere, as an Nx3 matrix
% -mag_Ap       Aparent magnitude of stars as an Nx1 matrix
% -type         Spectral Classification of stars as an Nx1 matrix using the
%               following coding: [O B A F G K M ?] = [1 2 3 4 5 6 7 8],
%               where '?' means unclassified
% -star_id      ID of each star in the hygdata_v3 database
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute additional camera parameters
fy = v_max / ( 2 * tan( h_fov * pi / 360 ));
fx = fy;
o_x = u_max/2;
o_y = v_max/2;

load(star_db_file);
num_stars = size(position,1);

% Estimate ratio of UV Brightness to Visible Spectrum Brightness based on
% star temperature:
class_temp = [35000, 20000, 9000, 6500, 5500, 4500, 3000, 3000]; %O B A F G K M ?
class_UV_ratio = zeros(size(class_temp));
for i = 1:size(class_temp,2)
    class_UV_ratio(i) = UV_fraction_by_T(class_temp(i));
end

%%
%Generate Random Rotation Matrix
tx = rand*2*pi;
Rx = [1,0,0;0,cos(tx),-sin(tx);0,sin(tx),cos(tx)];
ty = rand*2*pi;
Ry = [cos(ty),0,-sin(ty);0,1,0;sin(ty),0,cos(ty)];
tz = rand*2*pi;
Rz = [cos(tz),-sin(tz),0;sin(tz),cos(tz),0;0,0,1];
R = Rx * Ry * Rz;   %rotation from inertial frame to camera frame
R_c_i = R';         %rotation from camera frame to inertial frame

%%
%Generate Star Points in (u,v) pixel space:

num_in_image = 0;

%Store data about what stars appear in the image for ground truth file
im_list = zeros(num_stars, 2);          % Sub-pixel positions of stars in image in 2D pixel space
im_list_mag = zeros(num_stars, 1);      % Corresponding apparent magnitudes of stars in image
im_list_class = zeros(num_stars, 1);    % Corresponding classes of stars in image (using numerical coding)
im_list_3 = zeros(num_stars, 3);        % 3D points of included stars in inertial frame
im_list_3_cf = zeros(num_stars, 3);     % 3D points in camera reference frame

%Process all stars in database and check which appear in the image
for i = 1:num_stars
    p_w = position(i,:)';
    p_c = R * p_w;
    if p_c(3) > 0 % on correct side of image plane
        u = fx * p_c(1)/p_c(3) + o_x;
        v = fy * p_c(2)/p_c(3) + o_y;
        u_r = round(u);
        v_r = round(v);
        
        %Check if star appears in image. If so, add to im_list
        if u_r > 0 && u_r <= u_max && v_r > 0 && v_r <= v_max
            num_in_image = num_in_image + 1;
            im_list(num_in_image,:) = [u,v];
            im_list_mag(num_in_image) = mag_Ap(i);
            im_list_class(num_in_image) = type(i);
            im_list_3(num_in_image,:) = position(i,:);
            im_list_3_cf(num_in_image,:) = (R * position(i,:)')';
        end
    end
end
%truncate off zeros
im_list = im_list(1:num_in_image,:);
im_list_mag = im_list_mag(1:num_in_image,:);
im_list_class = im_list_class(1:num_in_image,:);
im_list_3 = im_list_3(1:num_in_image,:);
im_list_3_cf = im_list_3_cf(1:num_in_image,:);

%%
%Create Visible Spectrum and UV images
I1 = zeros(v_max,u_max); % vis
I2 = zeros(v_max,u_max); % uv
i = 1;
while i <= size(im_list,1)
    if im_list_mag(i) < -26
        disp('WARNING: Sol in image')
    else
        c = im_list(i,1);
        r = im_list(i,2);
        
        Fx1 = Fx0 * 10^( im_list_mag(i) / (-2.5) );
        B1 = Star_Blob(c,r,Fx1); % create star blob pattern for vis-image

        %Ratio of UV brightness to Visible-Spectrum brightness depends on
        %the star's class.
        Fx2 = Fx1 * class_UV_ratio(im_list_class(i));
        B2 = Star_Blob(c,r,Fx2); % create star blob pattern for uv-image
        
        mat_size = size(B1,1);
        try % try adding blob patterns into images
            I1(round(floor(r)+1-mat_size/2):round(floor(r)+mat_size/2), round(floor(c)+1-mat_size/2):round(floor(c)+mat_size/2)) = ...
               I1(round(floor(r)+1-mat_size/2):round(floor(r)+mat_size/2), round(floor(c)+1-mat_size/2):round(floor(c)+mat_size/2)) + B1;
            I2(round(floor(r)+1-mat_size/2):round(floor(r)+mat_size/2), round(floor(c)+1-mat_size/2):round(floor(c)+mat_size/2)) = ...
               I2(round(floor(r)+1-mat_size/2):round(floor(r)+mat_size/2), round(floor(c)+1-mat_size/2):round(floor(c)+mat_size/2)) + B2;            
            i = i+1; %advance
        catch
            %If a star is too close to the edge, it triggers an out of
            %bounds error. If this is caught, don't draw the stat, and
            %remove the star from im_list ground truth
            im_list(i,:) = [];
            im_list_mag(i) = [];
            im_list_class(i) = [];
            im_list_3(i,:) = [];
            im_list_3_cf(i,:) = [];
            continue
        end
    end
end

%%
%add noise
if addNoise
    I1 = imnoise(I1,'gaussian',0,0.005);
    I2 = imnoise(I2,'gaussian',0,0.005);
end

%%
%save files
name = strcat(prefix, '_image_vis.jpg');
imwrite(I1, name, 'jpg')
name = strcat(prefix, '_image_uv.jpg');
imwrite(I2, name, 'jpg')
name = strcat(prefix, '_truth.mat');
save(name,'im_list','im_list_class','im_list_mag','im_list_3','im_list_3_cf','R_c_i');
name = strcat(prefix, '_cp.mat');
save(name,'fx','fy','o_x','o_y');

%%
% display simulated images (optional)
if showIm
    figure;
    imshow(I1)
    figure;
    imshow(I2)
    % display vis-spectrum image with overlay
    figure;
    imshow(I1)
    hold on
    plot(im_list(:,1),im_list(:,2),'ro')
end