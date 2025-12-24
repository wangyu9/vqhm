function [L_sym] = symbolic_lap(n,F) 


II = [F(:,1);F(:,2);F(:,3);F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1);];
JJ = [F(:,1);F(:,2);F(:,3);F(:,2);F(:,3);F(:,1);F(:,1);F(:,2);F(:,3);];

L_sym = sparse(II,JJ,1,n,n); 