% This example shows how to calculate and plot both the
% fundamental TE and TM eigenmodes of an example 3-layer ridge
% waveguide using the full-vector eigenmode solver.  

% Refractive indices:
n1 = 3.34;          % Lower cladding
n2 = 3.44;          % Core
n3 = 1.00;          % Upper cladding (air)

% Layer heights:
h1 = 2.0;           % Lower cladding
h2 = 1.3;           % Core thickness
h3 = 0.5;           % Upper cladding

% Horizontal dimensions:
rh = 1.1;           % Ridge height
rw = 1.0;           % Ridge half-width
side = 1.5;         % Space on side

% Grid size:
dx = 0.0125*2.828;        % grid size (horizontal)
dy = 0.0125*2.828;        % grid size (vertical)

lambda = 1.55;      % vacuum wavelength
nmodes = 1;         % number of modes to compute

effective_array =[];


% First consider the fundamental TE mode:
i = 1;
c = (1-0.325)/10;
for index = 0.325:c:1
    [x,y,xc,yc,nx,ny,eps,edges] = waveguidemesh([n1,n2,n3],[h1,h2,h3], ...
                                            rh,index,side,dx,dy); 
    [Hx,Hy,neff] = wgmodes(lambda,n2,nmodes,dx,dy,eps,'000A');
    
    fprintf(1,'neff = %.6f\n',neff);
    effective_array = [effective_array, neff];
    figure(i);
    subplot(1,2,1);
    contourmode(x,y,Hx);
    title(['Hx (TE mode )', num2str(nmodes)]); xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end
    
    subplot(1,2,2);
    contourmode(x,y,Hy);
    title(['Hy (TE mode )',num2str(nmodes)]); xlabel('x'); ylabel('y'); 
    for v = edges, line(v{:}); end
    i = i + 1;
    % Next consider the fundamental TM mode
    % (same calculation, but with opposite symmetry)
    
    % [Hx,Hy,neff] = wgmodes(lambda,n2,mode_index,dx,dy,eps,'000S');
    % 
    % fprintf(1,'neff = %.6f\n',neff);
    % 
    % figure(mode_index);
    % subplot(2,2,3);
    % contourmode(x,y,Hx(:,:,mode_index));
    % title(['Hx (TM mode )',num2str(mode_index)]); xlabel('x'); ylabel('y'); 
    % for v = edges, line(v{:}); end
    % 
    % subplot(2,2,4);
    % contourmode(x,y,Hy(:,:,mode_index));
    % title(['Hy (TM mode )',num2str(mode_index)]); xlabel('x'); ylabel('y'); 
    % for v = edges, line(v{:}); end
end 
figure(i);
subplot(1,1,1);
plot(effective_array);
title("N-effective vs Iteration"); xlabel('Iteration'); ylabel('Neff');
