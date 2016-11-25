clear all
clc

%% k-vector construction

% Following construction in: 
% Mortari, Daniele, et al. "The pyramid star identification technique." Navigation 51.3 (2004): 171-183.
% APPENDIX: The k-vector range searching technique

load('feature_database.mat')
y = features(:,4); % Central Angles
[s,I] = sort(y);   % s - sorted central angles, I - index map to original data
n = length(y);
ymin = min(y);
ymax = max(y);

zeta = 2.2204e-16 * max(abs([ymin ymax]));
m = (ymax - ymin + 2*zeta)/(n - 1);
q = ymin - m - zeta;
x = 1:n;
z = m*x + q;     % line with steeper slope than the data

k_vector = zeros(n-1,1);

j = 1;
for i = 2:n-1

    while(z(i)>=s(j))
        j = j + 1;
    end


    k_vector(i) = j; % k-Vector construction
    
end

%% k-vector example

ang = 55 * pi/180;        % k-vector query
uncertanty = .001 * pi/180; % angular uncertanty

ya = ang - uncertanty;
yb = ang + uncertanty;

jb = floor((ya - q)/m);
jt = ceil((yb - q)/m);

kstart = k_vector(jb);
kend   = k_vector(jt)+1;

search_int = [y(I(kstart)) y(I(kend))];

disp(['For central angle: ' num2str(ang*180/pi) ' search range: ' num2str(search_int(1)*180/pi) ' to ' num2str(search_int(2)*180/pi)] )
disp(['Search element range: ' num2str(kstart) ' to ' num2str(kend)])

%% Hamming Distance
for i = kstart:kend

    disp(['Feature ' num2str(i) ' UV Vect: ' num2str(features(I(i),5:7))])

end