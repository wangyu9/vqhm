function [exjet] = my_colormap(map_type)

% remember to have: caxis([-0.2,1]);

if(strcmp(map_type,'weights-neg'))
    lambda1 = (1:140)'./140;
    lambda2 = 1 - lambda1;
    red = [1,0,0];
    green = [0,1,0];
    blue = [0,0,1];
    magenta = [1,0,1];

    map =  [...
        lambda1*blue+lambda2*magenta;... 
        lambda1*green+lambda2*blue;...
        lambda1*red+lambda2*green;...
        ];
    deepred_to_blue = jet(800);
    deepred_to_blue = deepred_to_blue(101:800,:);
    exjet = [...
        lambda1*blue+lambda2*magenta;... 
        deepred_to_blue;];
elseif(strcmp(map_type,'RdYlBu'))
    exjet = RdYlBu(800); 
    exjet = exjet(end:-1:1,:);
elseif strcmp(map_type,'heat')
	c = colormap('hot');
	exjet = c(end:-1:1,:);    
elseif strcmp(map_type,'halfHot')
	c = 0.77*colormap('hot');
	exjet = c;%c(end:-1:1,:);  
else 
    colormap(map_type);
    exjet = colormap;
    %error('Unknown colormap type');    
end



    

    

