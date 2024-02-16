
nx = 50;
ny = 50;
V = zeros(nx,ny);
G = sparse(nx*ny,nx*ny);
k=0;
Inclusion = 0;

for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        nxm = j + (i-2)*ny;
        nxp = j+(i)*ny;
        nym = j-1+(i-1)*ny;
        nyp = j+1+(i-1)*ny;
        if i == 1
            G(n,n) = 1;
        elseif i == nx
            G(n,n) = 1;
        elseif j == 1
            G(n,n) = 1;
        elseif j == ny
            G(n,n) = 1;
        elseif (i > 10 && i < 20 && j > 10 && j < 20)
            G(n,n) = -2;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        else
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end


    end
end

figure 

figure('name', 'Matrix')
spy(G)

nmodes = 9;
[E,D] = eigs(G,nmodes,'SM');

figure('name', 'EigenValues')
plot(diag(D),'*');

np = ceil(sqrt(nmodes))
figure('name','Modes')
for k = 1:1:nmodes
    M = E(:,k);
    for i = 1:nx
        for j = 1:ny
            n = i + (j-1)*nx;
            V(i,j) = M(n);
        end
        subplot(np,np,k), surf(V)
        title(['EV =' num2str(D(k,k))])
    end
end
