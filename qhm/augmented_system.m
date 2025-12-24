function [xxx] = augmented_system(H,dW,pre,L,LA,LB,RU,RB,eq_rhs)

if isempty(eq_rhs)
   RU = zeros(0,size(LA,1));
   RB = zeros(0,size(LB,1));
   eq_rhs = zeros([0,1]);
end

n = pre.n;
r = size(eq_rhs,1);
sp_zeros_like = @(AA) sparse([],[],[],size(AA,1),size(AA,2));
zeros_like = @(AA) zeros(size(AA));

%%
gg = dW(:);
assert(pre.dim==2);
% known: 2, unknown: 1.
id2 = [pre.known(:)', pre.known(:)' + n];
id1 = [pre.unknown(:)', pre.unknown(:)' + n];
%%
% H = [L,sp_zeros_like(L);sp_zeros_like(L),L] + HC;
H22 = H(id2,id2);
H11 = H(id1,id1);
H12 = H(id1,id2);
H21 = H(id2,id1);


LLA = [LA, sp_zeros_like(LA); sp_zeros_like(LA), LA];
LLB = [LB, sp_zeros_like(LB); sp_zeros_like(LB), LB];

g1 = gg(id1);
g2 = gg(id2);

ZZA = sp_zeros_like(LLA);
ZZB = sp_zeros_like(LLB);
Zg1 = zeros_like(g1);
ZRA = zeros([r, size(LLA,1)]);
ZRR = zeros([r,r]);
%%
if false

lhs = [...
     H22, H21, LLB', LLB', ZZB';...
     H12, ZZA,-LLA , ZZA , ZZA;...
     LLB,-LLA, ZZA , ZZA , ZZA;...
     LLB, ZZA, ZZA , ZZA , LLA;...
     ZZB, ZZA, ZZA,  LLA , H11;...
];
 
rhs = ...
    [g2;g1;Zg1;Zg1;Zg1];
end
%%

lhs = [...
     H22, H21, LLB',RB';...
     H12, H11,-LLA, RU';...
     LLB,-LLA, ZZA, ZRA';...
     RB,   RU, ZRA, ZRR;...
];
 
rhs = ...
    [g2;g1;Zg1;eq_rhs;];
%%
r = lhs \ rhs;

xxx = r(1:numel(g2));


end