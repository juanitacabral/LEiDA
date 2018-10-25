function plot_nodes_in_cortex (V)

hold on

% PLOT CORTEX

cortex.path='MNI152_T1_10mm_brain_mask.nii';
cortex.pial=mapPial(cortex.path);
cortex.color=[0.9 0.9 0.9];
cortex.transparency=0.3; % To view only opaque cortex =1;
cortex.val=0.3;
redux=0.5;
sregion=smooth3(cortex.pial);
psregion=patch(isosurface(sregion,cortex.val,'verbose'), 'FaceColor', cortex.color, 'EdgeColor', 'none');
reducepatch(psregion,redux,'verbose');
isonormals(sregion,psregion);
set(psregion,'FaceAlpha', cortex.transparency); %transparency

% PLOT NODES

V=V/max(abs(V));

% center origin
ori=[13.5 9.5 7.5];

load aal_cog.txt
MNI_coord=aal_cog/10;
clear aal_cog

[x,y,z] = sphere;

for n=1:length(V)
    if V(n)>0
        surf(x*V(n)+MNI_coord(n,2)+ori(1), y*V(n)+MNI_coord(n,1)+ori(2),z*V(n)+MNI_coord(n,3)+ori(3),'FaceColor','r','EdgeColor','none','FaceAlpha',0.5);
    elseif V(n)<0
        surf(x*V(n)+MNI_coord(n,2)+ori(1), y*V(n)+MNI_coord(n,1)+ori(2),z*V(n)+MNI_coord(n,3)+ori(3),'FaceColor','b','EdgeColor','none','FaceAlpha',0.5);
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