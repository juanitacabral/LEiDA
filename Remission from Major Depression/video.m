

%% Set up 3D plot to record
%figure;clf;
% >>>> COPY HERE THE CODE OF THE PICTURE THAT WILL BE THE FIRST FRAME OF THE VIDEO
figure
colormap(jet)
subplot(1,1,1)
    plot_nodes_in_cortex_new(V);

    camlight('headlight');

% surf(peaks,'EdgeColor','none','FaceColor','interp','FaceLighting','phong')
%daspect([1,1,.3]);axis tight;

%% Set up recording parameters (optional), and record
OptionZ.FrameRate=20;OptionZ.Duration=10;OptionZ.Periodic=true;
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'FCstate',OptionZ)
