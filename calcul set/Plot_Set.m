function Plot_Set(Z,P,c)
%Function to plot a set in each couple of dimension
% Input :   - Z is a set define with CORA
%           - P is a list of point
%           - c is the color of theset plot
%

if nargin == 1
    P = [];
    c = 'b';
elseif nargin == 2
    c = 'b';
end


n_z = Z.dimension;
nb_plot = n_z*(n_z-1)/2;
for i=1:nb_plot
    
     
    for j=1:n_z-1
        if i<=nb_plot-(n_z-j)*(n_z-1-j)/2
            x_2 = j;
            break
        end
    end
    x_1 = i -(nb_plot-(n_z-x_2)*(n_z-x_2+1)/2) + x_2;
    subplot(n_z,n_z,x_1+(x_2-1)*n_z)
    hold on
    plot(Z,[x_1,x_2],c);
    
    xlabel("Dim " + x_1)
    ylabel("Dim " + x_2)
    if ~isempty(P)
        plot(P(x_1),P(x_2),'ro')
    end
end
end