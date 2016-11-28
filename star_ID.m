%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computer Vision Final Project
% Phillip Rivera - priver14
% 
% Star-ID  -  Identifies the index of four stars with the provided body
% vectors and uv content. It utilizies the k-vector index database to
% obtain indices for internal angles between all body vectors. It has the
% option of utilizing uv image content to further reduce this search by
% discarding database results that do not match the uv content of earch
% star-blob.
%
% This process is achieved in 2 parts.
% Part 1 - Triangle List Generation
% The first three body vectors are used to generate a triangle, with the
% first entry selected as the pivot star. The k-vector search is called
% based on the provided body angle measurement, sigma, uv content and uv
% setting per body angle. There are two methods for using UV information.
% The strong threshold assumes that the UV content per blob is match
% perfectly with the UV content of the star-catalog (no false positive or
% false negative errors assumed). It eliminates the k-vector index pairs
% that do not have the same UV content. The not-strong threshold assumes
% that at least one of the two blobs UV content is correct, thus it retains
% the index results that have at least the same ammount of UV content as
% the blobs.
%
% Part 2 - Pyramid Identification
% For the list of triangles generated in Part 1, the fourth star is used to
% select which of the provided triangles might be a good match. In this
% portion, the benefit of having a small set of potential matches per blob
% is advantegios. For each blob in the triangle, a potential match is found
% with the forth blob. Only the result where all four stars are identified
% is outputed.
%
%
% Inputs
%   body_vectors:   4x3     vector matrix with a body vector for each blob 
%                           per row
%   uv_content:     4x1     binary vector 1: uv content, 0: no uv content 
%   sigma:          double  estimate of the standard deviation for the
%                           centroiding error
%   use_uv:         bool    use UV information to reduce k-vector result
%                           size 
%   strng_uv        bool    use strong UV thresholding if use_uv is true
%
% Outputs
% indeces:  nx4     matrix of identified indeces in the star catalog. n
%                   depends on the possible number of matches

function indices = star_ID(body_vectors,uv_content,sigma,use_uv,strng_uv)

%% Part 1 - Triangle List Generation
% Identify star catalog indices for the first three stars. Convention for
% index list per blob/star is 
% a -> star 1 (pivot)
% b -> star 2
% c -> star 3


% uv_triangle = sum of binary uv information for the first 3 stars
% v_triangle  = inner angles between vectors for the first 3 stars
% Index order is 1 -> star1 and star2
%                2 -> star1 and star3
%                3 -> star2 and star3

uv_triangle(1) = uv_content(1) + uv_content(2);
uv_triangle(2) = uv_content(1) + uv_content(3);
uv_triangle(3) = uv_content(2) + uv_content(3);

v_triangle(1) = acos(body_vectors(1,:)*body_vectors(2,:)');
v_triangle(2) = acos(body_vectors(1,:)*body_vectors(3,:)');
v_triangle(3) = acos(body_vectors(2,:)*body_vectors(3,:)');

% uv_pyramid = sum of binary uv information for the first 3 stars with the
%               4th star
% v_pyramid  = inner angles between vectors for the first 3 stars with the
%               4rth star
% Index order is 1 -> star1 and star4
%                2 -> star2 and star4
%                3 -> star3 and star4

uv_pyramid(1) = uv_content(1) + uv_content(4);
uv_pyramid(2) = uv_content(2) + uv_content(4);
uv_pyramid(3) = uv_content(3) + uv_content(4);

v_pyramid(1) = acos(body_vectors(1,:)*body_vectors(4,:)');
v_pyramid(2) = acos(body_vectors(2,:)*body_vectors(4,:)');
v_pyramid(3) = acos(body_vectors(3,:)*body_vectors(4,:)');

% Generate star pair hypothesis based on measured internal angles, UV
% content information and settings.
star_pairs = k_vector_search(v_triangle,sigma,uv_triangle,use_uv,strng_uv);

% Parse star pair data into two indices sets I and J for each angle in the
% triangle
I12 = star_pairs{1}(:,1); J12 = star_pairs{1}(:,2);
I13 = star_pairs{2}(:,1); J13 = star_pairs{2}(:,2);
I23 = star_pairs{3}(:,1); J23 = star_pairs{3}(:,2);

% Using the star pivot approach, select the potential indices for star 2 (b)
% and star 3 (c)
[~,i12,~] = intersect(I12,union(I13,J13));
[~,j12,~] = intersect(J12,union(I13,J13));
[~,i13,~] = intersect(I13,union(I12,J12));
[~,j13,~] = intersect(J13,union(I12,J12));

b = union(I12(j12),J12(i12));
c = union(I13(j13),J13(i13));

% Reduce the potential star 2 (b) indices by finding the intersection with the
% second angle, star 2 to star 3. Also provide potential indices for star 1
% (a)
b = intersect(b,union(I23,J23));
[~,~,j12b] = intersect(b,I12(j12));
[~,~,i12b] = intersect(b,J12(i12));
I12_resized = I12(i12);
J12_resized = J12(j12);

ab = union(I12_resized(i12b),J12_resized(j12b));

% Reduce the potential star 3 (c) indices by finding the intersection with the
% second angle, star 2 to star 3. Also provide potential indices for star 1
% (a)
c = intersect(c,union(I23,J23));
[~,~,j13c] = intersect(c,I13(j13));
[~,~,i13c] = intersect(c,J13(i13));
I13_resized = I13(i13);
J13_resized = J13(j13);
ac = union(I13_resized(i13c),J13_resized(j13c));

% Final list of indecies for star 1
a = intersect(ab,ac);

% Final list of indecies for star 2. Reduce list of star 2 with further
% reduced set of star 1 list.
[~,~,j12a] = intersect(a,I12(i12));
[~,~,i12a] = intersect(a,J12(j12));
J12_resized = J12(i12);
I12_resized = I12(j12);
b = intersect(b,union(J12_resized(j12a),I12_resized(i12a)));

% Final list of indecies for star 3. Reduce list of star 3 with further
% reduced set of star 1 list.
[~,~,j13a] = intersect(a,I13(i13));
[~,~,i13a] = intersect(a,J13(j13));
J13_resized = J13(i13);
I13_resized = I13(j13);
c = intersect(c,union(J13_resized(j13a),I13_resized(i13a)));

%% Part 2 - Pyramid Identification
% Finalize star-ID process by finding a match with the fourth star.

num_tgls = length(a)*length(b)*length(c); % number of possible triangles
star_triangle = zeros(num_tgls,3);        % list of possible triangles

nd = 1;                 % Counter for star-4 index
d = cell(num_tgls,1);   % List of possible indeces for star-4

for na = 1:length(a)
    for nb = 1:length(b)
        for nc = 1:length(c)
            
            star_triangle(nd,:) = [a(na) b(nb) c(nc)]; 
            
            d_star = [];      
           
            star_pairs = k_vector_search(v_pyramid,sigma,uv_pyramid,use_uv,strng_uv);
            
            for k = 1:3
                
                I4 = star_pairs{k}(:,1); J4 = star_pairs{k}(:,2);
                
                [~,~,i] = intersect(star_triangle(nd,k),I4);
                [~,~,j] = intersect(star_triangle(nd,k),J4);
                
                if(k >1)
                    d_star = intersect(d_star,union(I4(j),J4(i)));
                else
                    d_star = union(I4(j),J4(i));
                end
            end
            
            d{nd} = d_star;
            nd = nd+1;
            
        end
    end
end

% Find triangle hypothesis permutation that provided a match for the fourth
% star
ind_d = find(~cellfun('isempty', d));

% Output the indices of pyramid
indices = [star_triangle(ind_d,:) d{ind_d}];


function star_pairs_out = k_vector_search(v,sigma,uv,use_uv,strng_uv)

load('k_vector.mat')

star_pairs_out = cell(3,1);
for i = 1:3
    
    ang = v(i);               % k-vector query
    uncertanty = 6.4 * sigma; % angular uncertanty
    
    ya = ang - uncertanty;
    yb = ang + uncertanty;
    
    jb = max([floor((ya - q)/m) 1]);
    jt = ceil((yb - q)/m);
    
    kstart = k_vector(jb) + 1;
    kend   = k_vector(jt) + 1;
    
    index_pairs_mat = features(I(kstart:kend),2:3);
    
    if(use_uv)
        
        uv_mat = features(I(kstart:kend),4:5);
        uv_vect = sum(uv_mat,2);
        
        if(strng_uv)

            uv_ind_thrsh = find(uv_vect ~= uv(i));
        
        else
            uv_ind_thrsh = find(uv_vect < uv(i));
        end
        
        index_pairs_mat(uv_ind_thrsh,:) = []; 
    end
    
    star_pairs_out{i} = index_pairs_mat;
    
end



