function [KTs] = batch_grad_edgesq_over_vertex(VD, F)
% grad_edgesq_over_vertex

% extrinsic-to-intrinsic

assert(size(VD,2)==2);

f = size(F,1);

KTs = zeros(f, 6, 3);

v12 = VD(F(:,2),:) - VD(F(:,1),:);
v23 = VD(F(:,3),:) - VD(F(:,2),:);
v31 = VD(F(:,1),:) - VD(F(:,3),:);

KTs(:,1:2,2) = 2*v31; 
KTs(:,1:2,3) = 2*(-v12);

KTs(:,3:4,3) = 2*v12; 
KTs(:,3:4,1) = 2*(-v23);

KTs(:,5:6,1) = 2*v23; 
KTs(:,5:6,2) = 2*(-v31);

% double check
if false
KTs2 = zeros(f, 6, 3);

for i=1:f

    VTi = VD(F(i,:),:);
    KT = grad_edgesq_over_vertex(VTi);
    %KT: deodv, 6 x 3

    KTs2(i,:,:) = KT;

end 

sum(sum(sum(abs(KTs - KTs2))))
end