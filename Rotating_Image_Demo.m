% DEMO: Draw star field for a slowly rotating camera

% NOTE: File saving in Sim_Star_Image is supressed to speed up this script.
% (This is controlled by one of the parameters)

R_i_c = eye(3); %initial camera pose
inc = 4/360 * 2 * pi; %using 4 degree increments
R_inc = [cos(inc),0,-sin(inc);0,1,0;sin(inc),0,cos(inc)]; %incremental rotation

h_dim = 720;        % Horizontal image dimension
v_dim = 480;        % Vertical image dimension
h_fov = 20;         % Horizontal field of view in degrees

%%

figure;
for i = 1:90
    [I1, I2] = Sim_Star_Image('001', h_dim, v_dim, h_fov, true, true, false, true, R_i_c);    
    R_i_c = R_i_c * R_inc;
    imshow(I1)
    pause(0.2)
end
