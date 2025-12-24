function [beltrami] = tensor_para_forward(au,para_type)


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

switch para_type
    case 'complex-plane-det1'
        fff = @(w_rr) (1 - 1./(1+log(1+w_rr)));
        
    case 'complex-plane-det1-tanh'
        fff = @(w_rr)  tanh(w_rr);

    otherwise
        error('Error: Unimplemented!\n');
end


w1 = au(:,1);
w2 = au(:,2);

f = size(au,1); 

% fr1 = @(w1,w2) w1 ./ sqrt(w1.^2+w2.^2) .* ff(sqrt(w1.^2+w2.^2));
% fr2 = @(w1,w2) w2 ./ sqrt(w1.^2+w2.^2) .* ff(sqrt(w1.^2+w2.^2));

w_rr_sq = w1.^2+w2.^2;
w_rr = sqrt(w_rr_sq);


w_rr_aa = fff(w_rr); 

nw1 = w1 ./ w_rr;
nw2 = w2 ./ w_rr;

beltrami = [w_rr_aa.*nw1,  w_rr_aa.*nw2]; 