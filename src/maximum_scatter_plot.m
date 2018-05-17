function maximum_scatter_plot(X,fig,dim)

n = size(X);
if nargin == 2
    nx = n(1);
    ny = n(2);
    z = max(X,[],3);
    Z = zeros(nx,ny);
    Z(:,:) = z(:,:,1);
else
    if dim == 3
        nx = n(1);
        ny = n(2);
        z = max(X,[],3);
        Z = zeros(nx,ny);
        Z(:,:) = z(:,:,1);
    elseif dim == 1
        nx = n(2);
        ny = n(3);
        z = max(X,[],1);
        Z = zeros(nx,ny);
        Z(:,:) = z(1,:,:);
    else
        nx = n(1);
        ny = n(3);
        z = max(X,[],2);
        Z = zeros(nx,ny);
        Z(:,:) = z(:,1,:);
    end
end
x = zeros(nx,ny);
y = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        x(i,j) = i;
        y(i,j) = j;
    end
end
Z = double(Z);
Z = Z/max(max(Z))*100;
figure(fig)
Z = Z(:,ny:-1:1);
image(Z')
%imshow(Z')
colormap(parula(100))
%axis equal
axis tight
%colorbar