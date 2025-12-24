function [xxx] = augmented_system2(H,ga,pre,L,extra,RU,RB,eq_rhs)

if isempty(eq_rhs)
   eq_rhs = zeros([0,1]);
end

n = pre.n;
nu = numel(pre.unknown(:));
na = numel(ga);

r = size(eq_rhs,1);
sp_zeros_like = @(AA) sparse([],[],[],size(AA,1),size(AA,2));
zeros_like = @(AA) zeros(size(AA));

%%
assert(pre.dim==2);
% known: 2, unknown: 1.
id2 = [pre.known(:)', pre.known(:)' + n];
id1 = [pre.unknown(:)', pre.unknown(:)' + n];
%%
% H = [L,sp_zeros_like(L);sp_zeros_like(L),L] + HC;
LAB = L(pre.unknown(:),:);
LAB2 = [LAB, sp_zeros_like(LAB); sp_zeros_like(LAB), LAB];

Oaa = speye(na);
Zan2 = sparse([],[],[],na,n*2);
Zu2u2 = sparse([],[],[],nu*2,nu*2);

Zn2 = zeros(n*2,1);
Zu2 = zeros(nu*2,1);

assert(pre.dim==2);
T = [extra.Tx,extra.Ty];
%%

lhs = [...
         H,-LAB2',Zan2';...
     -LAB2, Zu2u2,   T';...
      Zan2,     T,  1e-6*Oaa;...
];
 
rhs = ...
    [Zn2;Zu2;ga];

%%
r = lhs \ rhs;

xxx = r(end-na+1:end);


end