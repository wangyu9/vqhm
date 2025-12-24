function [He1,He2,He3] = extrinsic_edge_hessians()

core = 2 * [1,0,-1,0;0,1,0,-1;-1,0,1,0;0,-1,0,1];

He3 = zeros(6,6);
He3(1:4,1:4) = core;

He1 = zeros(6,6);
He2 = zeros(6,6);

He1(3:6,3:6) = core;
He2([1,2,5,6],[1,2,5,6]) = core;


end
