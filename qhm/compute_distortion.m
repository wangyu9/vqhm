function compute_distortion(G,u,v)

f = size(G,1)/2;

n = numel(u);

assert(size(G,2)==n);


Gu = G * u;
Gv = G * v;

wf = zeros(f,1);

for j=1:f

    Jj = [Gu([j,j+f],1), Gv([j,j+f],2) ]'; % note the transpose


    
    wf(j) = 
end