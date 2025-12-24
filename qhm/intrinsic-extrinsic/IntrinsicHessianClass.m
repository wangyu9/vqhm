classdef IntrinsicHessianClass
   methods(Static)
       
        function gI=ComputeIntrinsicGradients(VD, F, V, g) % VD, deformed pose, V, restpose
            % 3 x f

            E0 = edge_lengths(V,F);
            
            n = size(V,1);
            f = size(F,1);
            
            E = edge_lengths(VD,F);
             
            % [He1,He2,He3] = extrinsic_edge_hessians();
            
%             IDX = [...
%                 obj.F(:,1),n+obj.F(:,1),...
%                 obj.F(:,2),n+obj.F(:,2),...
%                 obj.F(:,3),n+obj.F(:,3),...
%             ];

            %
            
            gI = zeros(3,f);
            
            for i=1:f
            
                a0 = E0(i,1).^2;
                a1 = E(i,1).^2;
                b0 = E0(i,2).^2;
                b1 = E(i,2).^2; 
                c0 = E0(i,3).^2; 
                c1 = E(i,3).^2;

                gi = g(a0,a1,b0,b1,c0,c1);
                % hi = h(a0,a1,b0,b1,c0,c1);

                idx = [...
                    F(i,1),n+F(i,1),...
                    F(i,2),n+F(i,2),...
                    F(i,3),n+F(i,3),...
                    ];

                gI(:,i) = gi;
                
            end
        end
        
        function HI=ComputeIntrinsicHessians(VD, F, V, h)
            % 3 x 3 x f
            
            %
            
            E0 = edge_lengths(V,F);
            
            n = size(V,1);
            f = size(F,1);
            
            E = edge_lengths(VD,F);
                        
            % [He1,He2,He3] = extrinsic_edge_hessians();
            
%             IDX = [...
%                 obj.F(:,1),n+obj.F(:,1),...
%                 obj.F(:,2),n+obj.F(:,2),...
%                 obj.F(:,3),n+obj.F(:,3),...
%             ];

            %
            
            HI = zeros(3,3,f);
            
            for i=1:f
            
                a0 = E0(i,1).^2;
                a1 = E(i,1).^2;
                b0 = E0(i,2).^2;
                b1 = E(i,2).^2; 
                c0 = E0(i,3).^2; 
                c1 = E(i,3).^2;

                % gi = g(a0,a1,b0,b1,c0,c1);
                hi = h(a0,a1,b0,b1,c0,c1);

                idx = [...
                    F(i,1),n+F(i,1),...
                    F(i,2),n+F(i,2),...
                    F(i,3),n+F(i,3),...
                    ];

                HI(:,:,i) = hi;
                
            end
        end
        
        function H=HessiansIntrinsic2Extrinsic(VD, F, V, Ai, HI, gI, assemble, project) % Ai: area of restpose.
            % 6 x 6 x f from 3 x 3 x f
                        
            %
            
            E0 = edge_lengths(V,F);
            
            n = size(V,1);
            f = size(F,1);
            
            % R = []
            E = edge_lengths(VD, F);
                        
            [He1,He2,He3] = extrinsic_edge_hessians();
            
            IDX = [...
                F(:,1),n+F(:,1),...
                F(:,2),n+F(:,2),...
                F(:,3),n+F(:,3),...
            ];
            %
            
            H_ele = zeros(6,6,f);
            
            for i=1:f
            
%                 a0 = E0(i,1).^2;
%                 a1 = E(i,1).^2;
%                 b0 = E0(i,2).^2;
%                 b1 = E(i,2).^2; 
%                 c0 = E0(i,3).^2; 
%                 c1 = E(i,3).^2;
% 
%                 gi = g(a0,a1,b0,b1,c0,c1);
%                 hi = h(a0,a1,b0,b1,c0,c1);

                hi = HI(:,:,i);
                gi = gI(:,i);

                idx = [...
                    F(i,1),n+F(i,1),...
                    F(i,2),n+F(i,2),...
                    F(i,3),n+F(i,3),...
                    ];

                VTi = VD(F(i,:),:);
                KT = grad_edgesq_over_vertex(VTi);
                %KT: deodv, 6 x 3

                switch project
                    case 0
                        Hi = ( KT * hi *KT' + He1 * gi(1) + He2 * gi(2) + He3 * gi(3) ) * Ai(i);
                    case 1
                        Hi = ( KT * hi *KT' + He1 * gi(1) + He2 * gi(2) + He3 * gi(3) ) * Ai(i);
                        [EV0,ED0] = eig((Hi+Hi')/2); 
                        Hi = EV0*max(ED0,0)*EV0';
                    case 2
                        Hi = ( KT * hi *KT' + He1 * max(gi(1),0) + He2 * max(gi(2),0) + He3 * max(gi(3),0) ) * Ai(i);
                    case 3
                        EHi = He1 * gi(1) + He2 * gi(2) + He3 * gi(3);
                        [EV,ED] = eig(EHi); 
                        % [EV2,ED2] = eig(hi); 
                        % Hi = ( KT * EV2*max(ED2,0)*EV2' *KT' + EV*max(ED,0)*EV' ) * obj.Ai(i);
                        Hi = ( KT * hi *KT' + EV*max(ED,0)*EV' ) * Ai(i);
                    otherwise
                        error(['unsupport project']);
                end
                                
                H_ele(:,:,i) = Hi;
            
            end
            
            if assemble
                H = sparse([],[],[],n*2,n*2);
                        
                for i=1:6
                   for j=1:6
                      EH = squeeze(H_ele(i,j,:));
                      H = H + sparse(IDX(:,i),IDX(:,j),EH,n*2,n*2);
                   end
                end 
            else
                H = H_ele;
            end
        end
        
        function gE=GradsIntrinsic2Extrinsic(VD, F, Area, gI)
            
            n = size(VD,1);
            f = size(F,1);
            
            dim = size(F,2) - 1;
            assert(dim==2);
            
            gE = zeros(n*2,1);
            
            IDX = [...
                F(:,1),n+F(:,1),...
                F(:,2),n+F(:,2),...
                F(:,3),n+F(:,3),...
            ];
        
            
            
            if false % old code
        
            for i=1:f
                
                VTi = VD(F(i,:),:);
                KT = grad_edgesq_over_vertex(VTi);
                %KT: deodv, 6 x 3
                
                gE(IDX(i,:),1) = gE(IDX(i,:),1) + KT * gI(:,i) * Area(i);
                
            end 
            
            else
                [KTs] = batch_grad_edgesq_over_vertex(VD, F);
                
                % KTs f x 6 x 3
                % gI' is f x 3 
                % gI  is 3 x f 
                
                % out: AreaKTgI is 6f or f x 6
                AgI = zeros(3,1,f);
                AgI(:,1,:) = gI;
                % AgI is 3x1xf
                
                if ~isunix

                                
                AreaKTgI = pagemtimes(permute(KTs,[2,3,1]), AgI);
                
                % sparseIDX is 2n x 6f 
                % IDX is f x 6
                II = IDX';
                II = II(:);
                VV = repmat(Area(:)', [6,1]);
                VV = VV(:);
                sparseIDX = sparse(II, 1:6*f, VV, 2*n, 6*f);
                
                gE = sparseIDX * AreaKTgI(:);
                
                else
                    % this is much slower.
                
                gE2 = zeros(n*2,1);
                for i=1:f
                
                    KT = squeeze(KTs(i,:,:));
                
                    gE2(IDX(i,:),1) = gE2(IDX(i,:),1) + KT * gI(:,i) * Area(i);
                
                end 
                
                % disp(norm(gE-gE2));
                gE = gE2;

                end
            end
            
        end
        
        function H=AssembleHessian(H_ele, n, F)
            
            dim = size(F,2) - 1;
            assert(dim==2);
            
            IDX = [...
                F(:,1),n+F(:,1),...
                F(:,2),n+F(:,2),...
                F(:,3),n+F(:,3),...
            ];
            
            for i=1:6
               for j=1:6
                  EH = squeeze(H_ele(i,j,:));
                  H = H + sparse(IDX(:,i),IDX(:,j),EH,n*2,n*2);
               end
            end 
        end

        function e=ValueWithWeights(VD, F, V, ef, W)
            
            e = 0;
            
            f = size(F,1);
            
            E0 = edge_lengths(V,F);
            
            n = size(V,1);
            f = size(F,1);

            assert(numel(W)==f);
            
            E = edge_lengths(VD,F);
             
            for i=1:f
            
                a0 = E0(i,1).^2;
                a1 = E(i,1).^2;
                b0 = E0(i,2).^2;
                b1 = E(i,2).^2; 
                c0 = E0(i,3).^2; 
                c1 = E(i,3).^2;

                ei = ef(a0,a1,b0,b1,c0,c1);
                
                e = e + ei * W(i);
                
            end
        end
        
        function gE=GradWithWeights(VD, F, V, g, W)
            Area = doublearea(V, F)/2;
            
            gI = IntrinsicHessianClass.ComputeIntrinsicGradients(VD, F, V, g);
            
            gE = IntrinsicHessianClass.GradsIntrinsic2Extrinsic(VD, F, Area .* W, gI);
            
        end

        function e=Value(VD, F, V, ef)
            
            e = 0;
            
            f = size(F,1);
            
            E0 = edge_lengths(V,F);
            
            n = size(V,1);
            f = size(F,1);
            
            E = edge_lengths(VD,F);
             
            for i=1:f
            
                a0 = E0(i,1).^2;
                a1 = E(i,1).^2;
                b0 = E0(i,2).^2;
                b1 = E(i,2).^2; 
                c0 = E0(i,3).^2; 
                c1 = E(i,3).^2;

                ei = ef(a0,a1,b0,b1,c0,c1);
                
                e = e + ei;
                
            end
        end
        
        function gE=Grad(VD, F, V, g)
            
            Area = doublearea(V, F)/2;
            
            gI = IntrinsicHessianClass.ComputeIntrinsicGradients(VD, F, V, g);
            
            gE = IntrinsicHessianClass.GradsIntrinsic2Extrinsic(VD, F, Area, gI);
            
        end
        
        function H=ProjectedHessian(VD, F, V, g, h)
            % [g,h] = intrinsic_grad_hessian(energy);
            
            gI = IntrinsicHessianClass.ComputeIntrinsicGradients(VD, F, V, g);
            HI = IntrinsicHessianClass.ComputeIntrinsicHessians(VD, F, V, h);
                        
            dim = size(F,2) - 1;
            assert(dim==2);
            Ai = doublearea(V, F)/2;
            
            project = 1;
            assemble = 1;
            H = IntrinsicHessianClass.HessiansIntrinsic2Extrinsic(VD, F, V, Ai, HI, gI, assemble, project);
            
        end
   end
end