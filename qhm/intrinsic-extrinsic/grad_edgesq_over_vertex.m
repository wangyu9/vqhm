function [deodv] = grad_edgesq_over_vertex(VX)

% extrinsic-to-intrinsic

% deodv: 6 x 3, specifying how vertex velocity determines ``edge sq velocity''
% de^2/dv
    deodv = zeros(6,3);
    
    deodv(1:2,[1,2,3]) = core(VX([1,2,3],:));
    deodv(3:4,[2,3,1]) = core(VX([2,3,1],:));
    deodv(5:6,[3,1,2]) = core(VX([3,1,2],:));
    
    % gv: 6 x 1
%     gv = [...
%         value(ge([1,2,3]), VX([1,2,3],:));
%         value(ge([2,3,1]), VX([2,3,1],:));
%         value(ge([3,1,2]), VX([3,1,2],:));
%     ];
    
end

function [deodv1] = core(VX)

    % deodv1: 2 x 3
    
    assert(size(VX,2)==2);

    v12 = VX(2,:) - VX(1,:);
    v23 = VX(3,:) - VX(2,:);
    v31 = VX(1,:) - VX(3,:);

    rnr = @(x) sqrt(sum(x.^2,2));

    e3 = rnr(v12);
    e1 = rnr(v23);
    e2 = rnr(v31);
        
    deodv1 = 2*[zeros(2,1),v31', (-v12)'];
    
    % e2 is unchanged so no effect. 
                       
end
