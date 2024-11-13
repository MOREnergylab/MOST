function Nemoh(projdir,plot_mesh_Nemoh)
fprintf('\n------ Starting NEMOH ----------- \n');
system(['preProc ',projdir]);

if plot_mesh_Nemoh
plot_total_mesh([projdir filesep 'mesh' filesep 'L12.dat'])
drawnow
end

fprintf('------ Solving BVPs ------------- \n');
system(['solver ',projdir]);
fprintf('------ Postprocessing results --- \n');
system(['postProc ',projdir]);
end
%% FUNCTIONS
function plot_total_mesh(filepath)
fid = fopen(filepath, 'rt');
coord_mat=[];
conn_mat=[];
zeroLineFound = false;
fgetl(fid);
tline=fgetl(fid);
while ischar(tline)
    numbers = str2num(tline);
    if all(numbers == 0)
        zeroLineFound = true;
    else
        if zeroLineFound
            conn_mat = [conn_mat; numbers];
        else
            coord_mat = [coord_mat; numbers];
        end
    end
    
    tline = fgetl(fid);
end
fclose(fid);

figure()
trimesh(conn_mat,coord_mat(:,2),coord_mat(:,3),coord_mat(:,4))
axis('equal')
title('Total Nemoh Mesh')
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

end
