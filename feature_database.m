clear all
close all

load('parsed_star_catalog.mat')
n = length(position);
features = zeros(n,7);      % Feature order:[Star1_ind Star2_ind Star3_ind central_angle Star1_isUV? Star2_isUV? Star3_isUV?]
features_binary = zeros(n,3);% Feature order:[Star1_ind central_angle any_star_UV?]
for i = 1:length(position)

    central_star = position(i,:);
    position_difference_angle = acos(position * central_star'); % internal angle between central star and all other vectors
    [~,nearest_ind] = sort(position_difference_angle); % nearest one its with itself, second and third nearest are what we want
    nearest_1 = nearest_ind(2); % Star 2 index
    nearest_2 = nearest_ind(3); % Star 3 index
    vector_1 = position(nearest_1,:) - central_star; 
    vector_2 = position(nearest_2,:) - central_star;
    central_angle = acos(dot(vector_1,vector_2)/(norm(vector_1)*norm(vector_2)));
    
    uv_bool = [isempty(intersect(ind_uv,i)),...
                isempty(intersect(ind_uv,nearest_1)),...
                 isempty(intersect(ind_uv,nearest_2))]; % Binary vector describing if Star 1, 2 or 3 are UV stars
    
    features(i,:) = [i nearest_1 nearest_2 central_angle uv_bool];
    features_binary(i,:) = [i central_angle sum(uv_bool)/max([sum(uv_bool),1])];
end

save('feature_database.mat','features','features_binary')


%% Angle Difference Error
figure;
plot(diff(sort(features(:,4)))*180/pi); grid on
title('Angular difference between nearest central angle')

[~,ang_ind] = sort(features(:,4));
hamming_distance = zeros(n-1);

for i = 1:n-1
    
    X = features(ang_ind(i),5:7);
    Y = features(ang_ind(i+1),5:7);
    hamming_distance(i) = pdist2(X,Y,'hamming');

end
figure;
plot(hamming_distance*3); grid on
title('Hamming Distance of UV 3-elemtn feature')


