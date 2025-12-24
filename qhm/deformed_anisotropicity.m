function [au_LD] = deformed_anisotropicity(pre,VD)
%
F = pre.F;

EL = pre.DEC.mesh.Elist;
GD = grad(VD,F);

MFD = spdiags(doublearea(VD,F)/2,0,pre.f,pre.f);

LD = GD' * blkdiag(MFD,MFD) * GD;

ne = size(EL,1);
au_LD = zeros(ne,1);
for i=1:ne
    au_LD(i,1) = - LD(EL(i,1),EL(i,2));
end
