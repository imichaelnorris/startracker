%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Arun's Algorithm for 3D point registration
% Authors: Christopher Hunt and Matthew Walmer
% Date: October 5, 2015
%
% Modified by Matthew Walmer on November 22, 2016
% -adapted to compute only the optimal rotation of two point clouds
%  assuming both frames have their origin in the same place
%
% * Arun's Algorithm fix comes from a paper from Stanford assessing 
% different graphics algorithms. Can be found at this URL: 
% http://graphics.stanford.edu/~smr/ICP/comparison/eggert_comparison_mva97.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters: A -> Points in camera frame (where each row is one point [x,y,z])
%             B -> Corresponding points in innertial frame B
% Return: R -> Optimal rotation matrix to rotate points from frame A to
%              frame B. If vA is a point in A and vB is a point in B, then
%              vB = R*vA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R = Calc_Rotation_Mat( A, B )
    %check dimensions
    if ( all( size(A) == size(B) ) ~= 1 ) %if any dimension does not agree
        error( 'Error: Frame dimensions must be equal!' )
    end
    if ( size(A, 2) > size(A, 1) )
        error( 'Error: Insufficient points! Amount < DOF' )
    end
    
    [ U, ~, V ] = svd( A' * B ); %H = A' * B
    R = V * U'; %optimal rotation
    if ( det(R) < 0 ) %check edge case, if R is a reflection
        disp( 'Warning: R was calculated as a reflection. Recalculating...' );
        V(:, end) = -V(:, end);
        R = V * U';
    end
end