function [V] = positions_from_edge_lengths(E)

% positions_from_edge_lengths(edge_lengths([0,0;1.21,0;0.3,0.37],[1,2,3]))
% V = positions_from_edge_lengths(edge_lengths([0,0;2.3,0;1.3,0.37;0,0;1,0;0.3,0.3],[1,2,3;4,5,6]))
% reshape(1:3*f,[3,f])'
% reshape(permute(V, [2,1,3]), [-1,2])

f = size(E,1);

V = zeros(f,3,2);

% V1 = zeros(f,2);

V2 = [E(:,3),zeros(f,1)];

V(:,2,:) = V2;

S = (E(:,1)+E(:,2)+E(:,3)) / 2;

% https://en.wikipedia.org/wiki/Heron%27s_formula

V3 = [(E(:,2).^2+E(:,3).^2-E(:,1).^2)./(2*E(:,3)), 2*sqrt(S.*(S-E(:,1)).*(S-E(:,2)).*(S-E(:,3)))./(E(:,3)) ];

V(:,3,:) = V3;

