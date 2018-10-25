function plot_nodes_in_cortex_new (V)

figure
hold on

% % PLOT CORTEX
% 
cortex.path='MNI152_T1_2mm_brain_mask.nii';
cortex.pial=mapPial(cortex.path);
cortex.color=[0.9 0.9 0.9];
cortex.transparency=0.3; % To view only opaque cortex =1;
cortex.val=0.2;
redux=1;
sregion=smooth3(cortex.pial);
psregion=patch(isosurface(sregion,cortex.val,'verbose'), 'FaceColor', cortex.color, 'EdgeColor', 'none');
reducepatch(psregion,redux,'verbose');
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', cortex.transparency); %transparency

% PLOT NODES

% Max of V or -V is 1
V=V/max(abs(V));

% center origin
ori=[65 45 36];

load aal_cog.txt
MNI_coord=aal_cog/10;
clear aal_cog

% PLOT LINKS AS CYLINDERS

% %For illustration we choose C as the Structural connectivity matrix
C=V'*V;
n_pos=find(V>0);
n_neg=find(V<0);

% scaling strongest connection to 1
C=(C/max(abs(C(:)))); 
scale=5.5;
hold on
% Plot edges as scaled cylinders according to values of matrix C
for n=n_pos
    for p = n_pos
        if C(n,p) > 0.001
            % The cylinder's diameter is scaled proportionally to C
            cdia=.1;
            % plot cylinder
            c1=[scale*MNI_coord(n,2)+ori(1) scale*MNI_coord(n,1)+ori(2) scale*MNI_coord(n,3)+ori(3)];
            c2=[scale*MNI_coord(p,2)+ori(1) scale*MNI_coord(p,1)+ori(2) scale*MNI_coord(p,3)+ori(3)];
            [X,Y,Z]=cylinder1(c1,c2,cdia,10);
            surf(X,Y,Z,'FaceColor',[1 C(n,p) 0],'EdgeColor','none','FaceAlpha', .5);
        end
    end
end

%PLOT gray links only if there is no n_pos
if isempty(n_pos)
    for n=n_neg
        for p = n_neg
            if C(n,p) > 0.9
                %The cylinder's diameter is scaled proportionally to C
                cdia=.02;
                %plot cylinder
                c1=[scale*MNI_coord(n,2)+ori(1) scale*MNI_coord(n,1)+ori(2) scale*MNI_coord(n,3)+ori(3)];
                c2=[scale*MNI_coord(p,2)+ori(1) scale*MNI_coord(p,1)+ori(2) scale*MNI_coord(p,3)+ori(3)];
                [X,Y,Z]=cylinder1(c1,c2,cdia,10);
                surf(X,Y,Z,'FaceColor',[.4 .4 .4],'EdgeColor','none','FaceAlpha', .5);
            end
end
end
end
       
   

[x,y,z] = sphere;
a=2;

for n=1:length(V)
    if V(n)>0
        surf(x*a+scale*MNI_coord(n,2)+ori(1), y*a+scale*MNI_coord(n,1)+ori(2),z*a+scale*MNI_coord(n,3)+ori(3),'FaceColor',[1 V(n) 0],'EdgeColor','none','FaceAlpha',1);
    elseif V(n)<0
        surf(x*a+scale*MNI_coord(n,2)+ori(1), y*a+scale*MNI_coord(n,1)+ori(2),z*a+scale*MNI_coord(n,3)+ori(3),'FaceColor',[0 -V(n) 1],'EdgeColor','none','FaceAlpha',1);
    end
end


axis off;
axis equal


% -------------------------------------------------------
% Setting image properties - light, material, angle
% -------------------------------------------------------
set(gcf,'Renderer', 'OpenGL') % USE UNDER LINUX FOR TRANSPARENCY 
view(3); axis off;
daspect([1 1 1]);
pbaspect([1 1 1]);
%daspect([0.90  1.37  0.90])
set(gca,'CameraViewAngle', 6);
set(gca, 'Projection', 'orthographic')
set(gca, 'CameraTarget', [51 68 90])
%set(gca, 'CameraPosition', [51 1360 90]) % saggital
%set(gca, 'CameraPosition', [1020 1360 1800]) % free angle
%set(gca, 'CameraUpVector', [1 0 0])

%view([-90 60])
%view([90 -90]) % ventral
view([-90 90]) % top
%view([0 0]) % R sideways
% view([-180 0]) % L sideways
% view([45 20]) % perspective
%view([90 0]) % front

material dull; lighting phong;
camlight;
rotate3d;

end

function pial=mapPial(region)

VG=spm_vol(region(1,:));
pial=zeros(VG.dim(1:3)); 
for i=1:VG.dim(3),
  pial(:,:,i) = spm_slice_vol(VG,spm_matrix([0 0 i]),VG.dim(1:2),1);
end

end