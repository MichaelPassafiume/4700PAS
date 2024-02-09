nx = 100;
ny = 100;
ni = 10000;
V = zeros(nx,ny);

SidesToZero = 1;
for k = 1:ni
    for i = 2:nx-1
        for j = 2:ny-1

            V(i, j) = (V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1))/4;
        end
    end
    
    if mod(k,50) == 0
       surf(V')
       pause(0.05)
    end
end

[Ex,Ey] = gradient(V);

figure
quiver(-Ey',-Ex',1)

