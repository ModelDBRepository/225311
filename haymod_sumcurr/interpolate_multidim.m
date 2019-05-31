function yi = interpolate(x,y,xi)

Ndim = size(y,1);
if any(xi < x(1)) || any(xi > x(end))
    error('extrapolation needed');
end

yi = zeros(Ndim,length(xi));
for i=1:length(xi)
    ind1 = find(x<=xi(i),1,'last');
    if ind1 == length(x)
        ind1 = ind1-1;
    end
    ind2 = ind1 + 1;
    yi(:,i) = y(:,ind1) + (y(:,ind2) - y(:,ind1))*(xi(i)-x(ind1))/(x(ind2)-x(ind1));
end
