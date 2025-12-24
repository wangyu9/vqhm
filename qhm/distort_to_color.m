function [wf] = distort_to_color(stats,encoder,args)


%%
ratio = stats.sigma_ratio;

assert(all(ratio>=1));


% encoder = 'log10';

%%

switch encoder
    case 'temp'
        
        base = 50;
        ratio(ratio> base) = base;
        
        ratio = log(ratio) / log(base);

    case 'log10'



        % log10 scale
        
        ratio(ratio>10) = 10;
        
        ratio = log10(ratio);

    case 'linear10'
        

        ratio = (ratio-1)/(10-1);

end

%%

wf = ratio; 