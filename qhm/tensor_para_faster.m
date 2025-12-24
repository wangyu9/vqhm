function [au,paupat] = tensor_para_faster(at,para_type)


% ff = @(rr) (1 - 1./(1+log(1+rr)));
% gg = @(rr) 1./(1+log(1+rr)).^2 ./ (1+rr);
% 
% ff = @(rr) tanh(rr);
% gg = @(rr) 1 - tanh(rr).^2;

% case 'complex-plane-det1'
%     ff = @(rr) (1 - 1./(1+log(1+rr)));
%     gg = @(rr) 1./(1+log(1+rr)).^2 ./ (1+rr);
%     res_sym_grad_complex_plane_det1;
% case 'complex-plane-det1-tanh'
%     ff = @(rr) tanh(rr);
%     gg = @(rr) 1 - tanh(rr).^2;
%     res_sym_grad_complex_plane_det1;

% s_at2au(at); 
% s_pdapdt_lmul = tp.s_pdapdt_lmul;


w1 = at(:,1);
w2 = at(:,2);

% fr1 = @(w1,w2) w1 ./ sqrt(w1.^2+w2.^2) .* ff(sqrt(w1.^2+w2.^2));
% fr2 = @(w1,w2) w2 ./ sqrt(w1.^2+w2.^2) .* ff(sqrt(w1.^2+w2.^2));

w_rr_sq = w1.^2+w2.^2;
w_rr = sqrt(w_rr_sq);
w_rr_cb = w_rr .* w_rr_sq;

nw1 = w1 ./ w_rr;
nw2 = w2 ./ w_rr;


% fff = ff( w_rr );
% ggg = gg( w_rr );

switch para_type
    case 'complex-plane-det1'
        fff = (1 - 1./(1+log(1+w_rr)));
        ggg = 1./(1+log(1+w_rr)).^2 ./ (1+w_rr);
        
    case 'complex-plane-det1-tanh'
        fff = tanh(w_rr);
        ggg = 1 - tanh(w_rr).^2;

    otherwise
        error('Error: Unimplemented!\n');
end

% fr1 = w1 ./ w_rr .* fff;
% fr2 = w2 ./ w_rr .* fff;

% let r = sqrt(x^2+y^2)
% ?( r ) / ?x = x/r
% ?( r ) / ?y = y/r
% ?( x/r ) / ?x = (r - x * x/r) / r^2 = (r^2-x^2)/r^3 = y^2/r^3 
% ?( x/r ) / ?y = ( - x * y/r) / r^2 = -x*y/r^3 
% ?( y/r ) / ?x = ( - y * x/r) / r^2 = -x*y/r^3 
% ?( y/r ) / ?y = (r - y * y/r) / r^2 = (r^2-y^2)/r^3 = x^2/r^3 

% let g(r) = f'(r)
% ?( x/r f(r)) / ?x = ( (f+xg)r - xf * x/r) / r^2 = (r^2(f+xg)-x^2f)/r^3 = (y^2f + r^2xg )/r^3 
% or 
% ?( x/r f(r)) / ?x = f (y^2/r^3 ) + g (x/r) * (x/r)
% ?( x/r f(r)) / ?y = f (-x*y/r^3 ) + g (x/r) * (y/r)
% ?( y/r f(r)) / ?x = f (-x*y/r^3 ) + g (y/r) * (x/r)
% ?( y/r f(r)) / ?y = f (x^2/r^3 ) + g (y/r) * (y/r)




% pfr11 = @(w1,w2) ff(sqrt(w1.^2+w2.^2)) .* (w2.^2 ./ (sqrt(w1.^2+w2.^2).^3)) ...
%             + gg(sqrt(w1.^2+w2.^2)) .* w1.^2 ./ (w1.^2+w2.^2); 
% 
% pfr12 = @(w1,w2) ff(sqrt(w1.^2+w2.^2)) .* (-w1.*w2 ./ (sqrt(w1.^2+w2.^2).^3)) ...
%             + gg(sqrt(w1.^2+w2.^2)) .* w1 .* w2 ./ (w1.^2+w2.^2); 


% pfr11 = fff .* (w2.^2 ./ w_rr_cb) + ggg .* w1.^2 ./ w_rr_sq; 
% pfr12 = fff .* (-w1.*w2 ./ w_rr_cb ) + ggg .* w1 .* w2 ./ w_rr_sq; 


pfr11 = fff .* (nw2.^2 ./ w_rr) + ggg .* nw1.^2; 
pfr12 = fff .* (-nw1.*nw2 ./ w_rr ) + ggg .* nw1 .* nw2; 


pfr21 = pfr12;
% we do not distinguish 12 and 21, which are same. 
        
% pfr22 = @(w1,w2) ff(sqrt(w1.^2+w2.^2)) .* (w1.^2 ./ (sqrt(w1.^2+w2.^2).^3)) ...
%             + gg(sqrt(w1.^2+w2.^2)) .* w2.^2 ./ (w1.^2+w2.^2); 

% pfr22 = fff .* (w1.^2 ./ w_rr_cb ) + ggg .* w2.^2 ./ w_rr_sq; 

pfr22 = fff .* (nw1.^2 ./ w_rr ) + ggg .* nw2.^2; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


p1 = w1 ./ w_rr .* fff;
p2 = w2 ./ w_rr .* fff;


rr_sq = p1.^2+p2.^2;
rr = sqrt(rr_sq);

% pmA1=...
%     @(p1,p2,p3)-((p1-1.0).^2+p2.^2)./(p1.^2+p2.^2-1.0);
% 
% pmA2=...
%     @(p1,p2,p3)(p2.*2.0)./(p1.^2+p2.^2-1.0);
% 
% pmA3=...
%     @(p1,p2,p3)-((p1+1.0).^2+p2.^2)./(p1.^2+p2.^2-1.0);


mA1=...
    -((p1-1.0).^2+p2.^2)./(rr_sq-1.0);

mA2=...
    (p2.*2.0)./(rr_sq-1.0);

mA3=...
    -((p1+1.0).^2+p2.^2)./(rr_sq-1.0);

au = [mA1, mA2, mA3];

% ps_dAdP11=...
%     @(p1,p2,p3)-(p1.*2.0-2.0)./(p1.^2+p2.^2-1.0)+p1.*((p1-1.0).^2+p2.^2).*1.0./(p1.^2+p2.^2-1.0).^2.*2.0;
% 
% ps_dAdP12=...
%     @(p1,p2,p3)(p2.*-2.0)./(p1.^2+p2.^2-1.0)+p2.*((p1-1.0).^2+p2.^2).*1.0./(p1.^2+p2.^2-1.0).^2.*2.0;
% 
% ps_dAdP13=...
%     @(p1,p2,p3)0.0;

ps_dAdP11=...
    -(p1.*2.0-2.0)./(rr_sq-1.0)+p1.*((p1-1.0).^2+p2.^2).*1.0./(rr_sq-1.0).^2.*2.0;

ps_dAdP12=...
    (p2.*-2.0)./(rr_sq-1.0)+p2.*((p1-1.0).^2+p2.^2).*1.0./(rr_sq-1.0).^2.*2.0;

ps_dAdP13=...
    0.0;


% ps_dAdP21=...
%     @(p1,p2,p3)p1.*p2.*1.0./(p1.^2+p2.^2-1.0).^2.*-4.0;
% 
% ps_dAdP22=...
%     @(p1,p2,p3)2.0./(p1.^2+p2.^2-1.0)-p2.^2.*1.0./(p1.^2+p2.^2-1.0).^2.*4.0;
% 
% ps_dAdP23=...
%     @(p1,p2,p3)0.0;


ps_dAdP21=...
    p1.*p2.*1.0./(rr_sq-1.0).^2.*-4.0;

ps_dAdP22=...
    2.0./(rr_sq-1.0)-p2.^2.*1.0./(rr_sq-1.0).^2.*4.0;

ps_dAdP23=...
    0.0;


% ps_dAdP31=...
%     @(p1,p2,p3)-(p1.*2.0+2.0)./(p1.^2+p2.^2-1.0)+p1.*((p1+1.0).^2+p2.^2).*1.0./(p1.^2+p2.^2-1.0).^2.*2.0;
% 
% ps_dAdP32=...
%     @(p1,p2,p3)(p2.*-2.0)./(p1.^2+p2.^2-1.0)+p2.*((p1+1.0).^2+p2.^2).*1.0./(p1.^2+p2.^2-1.0).^2.*2.0;
% 
% ps_dAdP33=...
%     @(p1,p2,p3)0.0;


ps_dAdP31=...
    -(p1.*2.0+2.0)./(rr_sq-1.0)+p1.*((p1+1.0).^2+p2.^2).*1.0./(rr_sq-1.0).^2.*2.0;

ps_dAdP32=...
    (p2.*-2.0)./(rr_sq-1.0)+p2.*((p1+1.0).^2+p2.^2).*1.0./(rr_sq-1.0).^2.*2.0;

ps_dAdP33=...
    0.0;

 
% @(p1,p2) p2.^2 ./ (sqrt(p1.^2+p2.^2).^3); 
% @(p1,p2) - p1.*p2 ./ (sqrt(p1.^2+p2.^2).^3); 
% @(p1,p2) p1.^2 ./ (sqrt(p1.^2+p2.^2).^3); 



% mA1=...
%     @(p1,p2,p3) pmA1(fr1(p1,p2), fr2(p1,p2), p3);
% 
% mA2=...
%     @(p1,p2,p3) pmA2(fr1(p1,p2), fr2(p1,p2), p3);
% 
% mA3=...
%     @(p1,p2,p3) pmA3(fr1(p1,p2), fr2(p1,p2), p3);
% 

np = 3;

paupat = struct(); 

paupat.s_dAdP1 = cell(1,np);
paupat.s_dAdP2 = cell(1,np);
paupat.s_dAdP3 = cell(1,np);

if true

% 
% 
% s_dAdP1{1}=...
%     @(p1,p2,p3) ps_dAdP11(fr1(p1,p2), fr2(p1,p2), p3) .* pfr11(p1,p2) ...
%               + ps_dAdP12(fr1(p1,p2), fr2(p1,p2), p3) .* pfr21(p1,p2);
%     
% s_dAdP1{2}=...
%     @(p1,p2,p3) ps_dAdP11(fr1(p1,p2), fr2(p1,p2), p3) .* pfr12(p1,p2) ...
%               + ps_dAdP12(fr1(p1,p2), fr2(p1,p2), p3) .* pfr22(p1,p2);
%           
% s_dAdP1{3}=...
%     @(p1,p2,p3) ps_dAdP13(p1,p2,p3);

s_dAdP1{1}=...
    ps_dAdP11 .* pfr11 + ps_dAdP12.* pfr21;
    
s_dAdP1{2}=...
    ps_dAdP11 .* pfr12 + ps_dAdP12.* pfr22; 
          
s_dAdP1{3}=...
    ps_dAdP13; 

%     
% s_dAdP2{1}=...
%     @(p1,p2,p3) ps_dAdP21(fr1(p1,p2), fr2(p1,p2), p3) .* pfr11(p1,p2) ...
%               + ps_dAdP22(fr1(p1,p2), fr2(p1,p2), p3) .* pfr21(p1,p2);
% 
% s_dAdP2{2}=...
%     @(p1,p2,p3) ps_dAdP21(fr1(p1,p2), fr2(p1,p2), p3) .* pfr12(p1,p2) ...
%               + ps_dAdP22(fr1(p1,p2), fr2(p1,p2), p3) .* pfr22(p1,p2);          
%           
% s_dAdP2{3}=...
%     @(p1,p2,p3) ps_dAdP23(p1,p2,p3);
% 


s_dAdP2{1}=...
    ps_dAdP21 .* pfr11 + ps_dAdP22 .* pfr21;

s_dAdP2{2}=...
    ps_dAdP21 .* pfr12 + ps_dAdP22 .* pfr22;       
          
s_dAdP2{3}=...
    ps_dAdP23;


% 
% s_dAdP3{1}=...
%     @(p1,p2,p3) ps_dAdP31(fr1(p1,p2), fr2(p1,p2), p3) .* pfr11(p1,p2) ...
%               + ps_dAdP32(fr1(p1,p2), fr2(p1,p2), p3) .* pfr21(p1,p2);
% 
% s_dAdP3{2}=...
%     @(p1,p2,p3) ps_dAdP31(fr1(p1,p2), fr2(p1,p2), p3) .* pfr12(p1,p2) ...
%               + ps_dAdP32(fr1(p1,p2), fr2(p1,p2), p3) .* pfr22(p1,p2);          
%           
% s_dAdP3{3}=...
%     @(p1,p2,p3) ps_dAdP33(p1,p2,p3);
% 


s_dAdP3{1}=...
    ps_dAdP31 .* pfr11 + ps_dAdP32 .* pfr21; 

s_dAdP3{2}=...
    ps_dAdP31 .* pfr12 + ps_dAdP32 .* pfr22;         
          
s_dAdP3{3}=...
    ps_dAdP33; 

end

paupat.s_pdapdt = @(gg,Area) eval_grad(gg,s_dAdP1, s_dAdP2, s_dAdP3, Area);


end


function [rg] = eval_grad(gg,s_1,s_2,s_3, Area)
 
    f = numel(gg) / 3;
    m = 3; 

    gg = reshape(gg, [f,3]);

    
    rg = zeros(f*m, 1);
    
    for ii=1:m
        rg(1+(ii-1)*f:ii*f) = ...
            Area .* ( s_1{ii} .* gg(:,1) + s_2{ii} .* gg(:,2) + s_3{ii} .* gg(:,3) );
    end
end