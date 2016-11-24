%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Star Catalog Parser Full
%
% This parser processes the hygdata_v3 database. It outputs a reduced star 
% catalog to 'parsed_star_db_full.mat' with the following data:
% -position     3D positions of stars on the unit sphere, as an Nx3 matrix
% -mag_Ap       Aparent magnitude of stars as an Nx1 matrix
% -type         Spectral Classification of stars as an Nx1 matrix using the
%               following coding: [O B A F G K M ?] = [1 2 3 4 5 6 7 8],
%               where '?' means unclassified
% -star_id      ID of each star in the hygdata_v3 database
%
% Unlike star_catalog_parser.m, the catalog output by this parser includes
% every star in the hygdata_v3 database.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('hygdata_v3.mat')

star_id   = id;        % Mapping to original table
position  = [x y z]./repmat(sqrt(x.^2 + y.^2 + z.^2),1,3); %location in unit sphere
mag_Ap    = mag;       % Apparent magnitude M_v

%Convert Star Classifications to Numerical Code:
%[O B A F G K M ?] = [1 2 3 4 5 6 7 8]
%class 8 is for anything else, and things without a known class. We'll
%assume these stars have no UV signature
type = zeros(size(mag_Ap));
for i = 1:size(type,1)
    if size(spect{i},2) > 0
        l = spect{i}(1);
        if l == 'O'
            type(i) = 1;
        elseif l == 'B'
            type(i) = 2;
        elseif l == 'A'
            type(i) = 3;
        elseif l == 'F'
            type(i) = 4;
        elseif l == 'G'
            type(i) = 5;
        elseif l == 'K'
            type(i) = 6;
        elseif l == 'M'
            type(i) = 7;
        else
            type(i) = 8;
        end
    else
        type(i) = 8;
    end
end
save('parsed_star_db_full.mat','star_id','position','mag_Ap','type');