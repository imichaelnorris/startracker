clear all
clc

%% k-vector construction

% Following construction in: 
% Mortari, Daniele, et al. "The pyramid star identification technique." Navigation 51.3 (2004): 171-183.
% APPENDIX: The k-vector range searching technique

load('feature_database_anglePairs.mat')
y = features(:,1); % Angles difference
[s,I] = sort(y);   % s - sorted angle difference, I - index map to original data
n = length(y);
ymin = min(y);
ymax = max(y);

zeta = 2.2204e-16 * max(abs([ymin ymax]));
m = (ymax - ymin + 2*zeta)/(n - 1);
q = ymin - m - zeta;
x = 0:n-1;
z = m*x + q;     % line with steeper slope than the data

k_vector = zeros(n,1);

j = 0;
for i = 2:n

    while(z(i)>=s(j + 1))
        j = j + 1;
    end


    k_vector(i) = j; % k-Vector construction
    
end

save('k_vector.mat','features','s','I','q','m','z','k_vector')