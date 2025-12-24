function [VT] = tutte_init(V,F,VTI,graph)


n = size(V,1);
f = size(F,1);


if graph
    fprintf('Init with uniform graph Laplacian\n');
    GI = intrinsic_grad_equilateral(n,F);
    Area = ones([size(F,1),1]);  
else
    GI = intrinsic_grad(V,F);
    Area = doublearea(V,F)/2;
end


MF = sparse(1:2*f,1:2*f,[Area;Area],2*f,2*f);

L = GI' * MF * GI;


[BE, ~, ~] = circular_laplacian(size(V,1),F);
B = sort(BE(:,1));

UB = find(~sparse(1,B,true,1,n));
UB = UB(:);


VT = VTI(:,1:2);

VT(UB,:) =  - L(UB,UB) \ ( L(UB,B) * VT(B,:) );


