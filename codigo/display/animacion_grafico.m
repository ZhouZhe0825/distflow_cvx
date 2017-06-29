function [] = animacion_grafico(z3, xlimV, ylimV, zlimV, viewV, pauseTime, titleGraph)
% Z = peaks;
% surf(Z)
% axis tight
% set(gca,'nextplot','replacechildren','visible','off')
% f = getframe;
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,20) = 0;
% for k = 1:20
%   surf(cos(2*pi*k/20)*Z,Z)
%   f = getframe;
%   size(rgb2ind(f.cdata,map,'nodither'))
%   im(:,:,1,k) = rgb2ind(f.cdata,map,'nodither');
% end
% imwrite(im,map,'DancingPeaks.gif','DelayTime',0,'LoopCount',inf)

figure('Position', [0, 80, 1024, 768])

%// Generate data
Top = size(z3,3);
%// Create handles to access/modify data.
hSurf = surf(z3(:,:,1));
xlim(xlimV);
ylim(ylimV);
zlim(zlimV);
view(viewV);
title([titleGraph, ' - Iteracion ', num2str(1)]);

k = 1;

pause(10);

%// Set up name to create animated gif.
%// Just a loop
while k <= Top

    %// IMPORTANT part. Update the Z data
    Z = z3(:,:,k);
    set(hSurf,'ZData',Z);

    %// Set limits so the graph looks nice.
    xlim(xlimV);
    ylim(ylimV);
    zlim(zlimV);
    view(viewV);
    title([titleGraph, ' - Iteracion ', num2str(k)]);
    drawnow

%     %// Capture frame to write to gif.
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind,cm] = rgb2ind(im,256);
%     if k == 1;
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
% 
    pause(pauseTime)

    k = k+1;
end