function [Li,se] = plot_JackKnife(arg,col,solid,TR)
%adapted from JackKnife written by Chandrasekaran used by Erie
%stuff to plot
if solid
    linest = '-';
else
    linest = '--';
end
y = nanmean(arg,1);

x = [1/(10/TR):1/(10/TR):length(y)/(10/TR)];
se = nanstd((arg))./sqrt(size(arg,1)-1);
L = y - se;
U = y + se;
Xcoords = [x x(end:-1:1)];
Ycoords = [U L(end:-1:1)];

% % Helen's additions !!
xlabel('Time (s)');
ylabel('Effect size');

%shadow color - add white noise (0.8 * diff to [1 1 1])
col2 = col;
col2(find(col<1)) = col2(find(col<1))+ 0.6*[1-col(find(col<1))];

% copied and adapted from JackKnife(x,y,se,col,col2);
hold on;
Pa = patch(Xcoords,Ycoords,col2);
set(Pa,'linestyle','none','linewidth',1);
Li = plot(x,y,'color',col,'linewidth',2,'LineStyle',linest);
alpha(0.5)


