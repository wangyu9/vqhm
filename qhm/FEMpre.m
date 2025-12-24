classdef FEMpre < handle
  properties
    dim
    V
    F
    known
    BC
    f
    G
    FA
    MF
    J
    v2f
    n
    unknown
    L
    Mass
    invMass
    old_BR
    old_Wu
    s_pdapdt_lmul
    aG
    iter
    DEC
    D
    star
    use_dir_edge
  end
  
  methods
      
        function pre = FEMpre(V,F,known,BC,s_pdapdt_lmul)

        % pre = struct();
        
        pre.use_dir_edge = false;

        pre.dim = size(F,2) - 1;

        pre.V = V(:,1:pre.dim);

        pre.F = F;
        pre.known = known;
        pre.BC = BC;

        pre.f = size(pre.F,1);
        pre.G = my_grad(pre.V,pre.F);

        if pre.dim==2
            pre.FA = doublearea(pre.V,pre.F)/2;
        else
            assert(pre.dim==3);
            pre.FA = volume(pre.V,pre.F);
        end
        if(min(pre.FA)<0)
           warning("The area FA has negative terms, made positive.\n"); 
           fprintf("%d elements have negative volume.\n", length(find(pre.FA<0)));
           pre.FA = abs(pre.FA);
        end
        pre.MF = spdiags(pre.FA,0,pre.f,pre.f);
        pre.J = vertex_facet_adjacency(pre.F);
        pre.v2f = pre.J./(sum(pre.J,2));
        pre.n = size(pre.V,1);
        pre.unknown = find(~sparse(1,pre.known,true,1,pre.n));

        if isempty(pre.BC)
           pre.BC = eye(length(pre.known)); 
        end

        pre.L = - cotmatrix(pre.V,pre.F);
        pre.Mass = massmatrix(pre.V,pre.F,'voronoi');
        if(min(diag(pre.Mass))<0)
           warning("The area Mass has negative terms, made positive.\n"); 
           fprintf("%d vertices have negative volume.\n", length(find(diag(pre.Mass)<0)));
           pre.Mass = abs(pre.Mass);
        end
        assert(nnz(pre.Mass)==pre.n);
        pre.invMass = sparse(1:pre.n,1:pre.n,1./diag(pre.Mass));

        f = pre.f;
        % assert(pre.dim==2);
        % pre.aG = [pre.G(f+1:2*f,:);pre.G(1:f,:)];
        
        %% only for the DEC and edge-based approach
        
        if size(pre.F,2)==3
        
        X = pre.V;
        if size(X,2)==2
           X = [X,zeros(size(X,1),1)]; 
        end
        mesh = getMeshData(X,pre.F,2,'none');
        pre.DEC = discreteExteriorCalculus(mesh);
        pre.DEC.mesh = mesh;
        pre.D = pre.DEC.d01;
        pre.star = pre.DEC.star11D;
        
        end
        %

        %%

        pre.old_BR = [];
        pre.old_Wu = [];

        pre.s_pdapdt_lmul = s_pdapdt_lmul;
        
        pre.iter = 0;
        end
  end
end
