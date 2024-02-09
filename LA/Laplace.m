nx = 100;
ny = 100;
ni = 7000;
V = zeros(nx,ny);
V(:,ny) = 0; % Top
V(:,1) = 0; % Bottom
V(1,:) = 1; % Left
V(nx,:) = 1; % Right
V_temp = V;


for k = 1:ni
    for i = 1:nx
        for j = 1:ny
            if i == 1 %Left
                V(i, j) = V(i+1, j);
            elseif i == nx %Right
                V(i, j) = V(i-1, j);
            elseif j == 1 %Bottom Set to one above
                V(i, j) = V(i, j+1);
            elseif j == ny %Top Set to one bellow
                V(i, j) = V(i, j-1);
            else
                V(i, j) = (V(i-1,j)+V(i+1,j)+V(i,j-1)+V(i,j+1))/4;
                
            end
            V(:,ny) = 0; %Reset Top
            V(:,1) = 0; %Reset Bottom
            V(1,:) =1; %Reset Left
            V(nx,:) = 1; %Reset Right
        end
    end
    
    if mod(k,50) == 0
       surf(V')
       pause(0.05)
    end
    V_temp = imboxfilt(V_temp,3);
    V_temp(1,:) =0; %Reset Left
    V_temp(nx,:) = 1; %Reset Right
end

[Ex,Ey] = gradient(V);

figure
quiver(-Ey',-Ex',1)
figure 
surf(V_temp')

