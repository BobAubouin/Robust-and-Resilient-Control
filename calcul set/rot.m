function R = rot(theta)
%This function create a roration matrix in dimension 2
% Input :   - theta : angle (rad)
%
% Output:   - R: Rotation matrix
%
% Bob Aubouin, 21/09/2021
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
end