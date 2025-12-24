function [dvode] = grad_vertex_over_edge(VX)

% intrinsic-to-extrinsic

% dvode: 3 x 6, specifying one possibility of how   
% ``edge velocity'' determines vertex velocity. 
% dv/de: not unique, only give one solution. 

dvode = ones(3,6);
dvode(3,[1,2,3,4,5,6]) = core(VX([1,2,3],:));
dvode(1,[3,4,5,6,1,2]) = core(VX([2,3,1],:));
dvode(2,[5,6,1,2,3,4]) = core(VX([3,1,2],:));


% ge = ones(1,3);
% 
% ge(1,3) = value(gv([1,2,3],:), VX([1,2,3],:));
% ge(1,1) = value(gv([2,3,1],:), VX([2,3,1],:));
% ge(1,2) = value(gv([3,1,2],:), VX([3,1,2],:));

end


function [dvode3] = core(VX)

% the convention is used that 1,2,3 is counter-clock-wise.
% but the code works also in the case of inversion, i.e. 1,2,3 becomes
% counter-clock-wise, in which case negative angles pop up. 

    assert(size(VX,2)==2);

    v13 = VX(3,:) - VX(1,:);
    v13R = [-v13(:,2), v13(:,1)]; % v13 rotated counter-clockwise
    
    v32 = VX(2,:) - VX(3,:);
    v32R = [-v32(:,2), v32(:,1)];
    
    v12 = VX(2,:) - VX(1,:);
    
    rn = @(x) x./sqrt(sum(x.^2,2));
    
    v13RN = rn(v13R);
    v32RN = rn(v32R);
    
    v12N = rn(v12);
    v32N = rn(v32);
    
    sinangle2 = v12N(:,2) * v32N(:,1) -  v12N(:,1) * v32N(:,2) ; % sinangle2 is negative if the triangle indexing is not counter-clockwise 
    % sinangle2 = sin (<021 - < 023)
    
    dvode3 = [0,0, v32RN / sinangle2, 0, 0]; % one possibility.
    
    % (v12N .* v32RN) is cos angle < 123. % wangyu has checked this with rigorous calculations. 
    
    % de/dv1 / (de/dl) = dl/dv1 = sin < 123
    % (de/dl) = de/dv1 / sin < 123
    % de/dv1 = g2 .* v32RN;
        
end

function [ge] = value(gv,VX)

% the convention is used that 1,2,3 is counter-clock-wise.
% but the code works also in the case of inversion, i.e. 1,2,3 becomes
% counter-clock-wise, in which case negative angles pop up. 

    assert(size(VX,2)==2);

    g1 = gv(1,:);
    g2 = gv(2,:);
    g3 = gv(3,:);

    ge = zeros(3,size(gv,2));

    v13 = VX(3,:) - VX(1,:);
    v13R = [-v13(:,2), v13(:,1)]; % v13 rotated counter-clockwise
    
    v32 = VX(2,:) - VX(3,:);
    v32R = [-v32(:,2), v32(:,1)];
    
    v12 = VX(2,:) - VX(1,:);
    
    rn = @(x) x./sqrt(sum(x.^2,2));
    
    v13RN = rn(v13R);
    v32RN = rn(v32R);
    
    v12N = rn(v12);
    v32N = rn(v32);
    
    sinangle2 = v12N(:,2) * v32N(:,1) -  v12N(:,1) * v32N(:,2) ; % sinangle2 is negative if the triangle indexing is not counter-clockwise 
    % sinangle2 = sin (<021 - < 023)
    
    ge = g2 .* v32RN / sinangle2; 
    % (v12N .* v32RN) is cos angle < 123. % wangyu has checked this with rigorous calculations. 
    
    % de/dv1 / (de/dl) = dl/dv1 = sin < 123
    % (de/dl) = de/dv1 / sin < 123
    % de/dv1 = g2 .* v32RN;
    
    % as a double check, ge can be obtain by the following: TODO
    
end