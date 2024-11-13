function NemohRun_from_Salome(Platform,data)
%% Folder and file name input
folderMesh = 'Mesh'; 
folderResults = data.nemoh.results_folder_name; 
BodyNames = fields(Platform); 
nBodies=length(BodyNames);
RefP=zeros(nBodies,3);
for i=1:nBodies
    RefP(i,:)=Platform.(BodyNames{i}).COG'+[data.geom.locations(i,1:2) 0];
end
%% Conversion of the mesh file format from SALOME to Nemoh and rototranslation
nodes = zeros(1,nBodies);
panels = zeros(1,nBodies);
for ii=1:nBodies
result = dlmread([folderMesh,filesep,BodyNames{ii},'.dat']); %#ok<*DLMRD>
fid = fopen ([folderMesh,filesep,BodyNames{ii},'_Nemoh.dat'], 'w');
fprintf(fid, '%g %g\n',[2 0]);
for jj = 1:result(1,1)
   fprintf(fid,'%g %f %f %f\n',[result(jj+1,1) (Rz(data.geom.locations(ii,3))*result(jj+1,2:4)')'+[data.geom.locations(ii,1:2) 0]]);     
end
rows203 = find(result(:,2)==203); %index where the connectivity table refers to 3 nodes triangle reference
rows204 = find(result(:,2)==204); %index where the connectivity table refers to 4 nodes quadrangular reference
fprintf(fid, '%g %g %g %g\n',[0 0 0 0]);
for kk = 1:length(rows203)
   fprintf(fid,'%g %g %g %g\n',result(rows203(kk),3),result(rows203(kk),4),result(rows203(kk),5),result(rows203(kk),3));     
end
for kk = 1:length(rows204)
   fprintf(fid,'%g %g %g %g\n',result(rows204(kk),3),result(rows204(kk),4),result(rows204(kk),5),result(rows204(kk),6));     
end
fprintf(fid, '%g %g %g %g\n',[0 0 0 0]);
fclose(fid);
nodes(1,ii) = result(1,1);
panels(1,ii) = length(rows203)+length(rows204);
end
%% Write Nemoh.cal file
fid=fopen('Nemoh.cal','w');
fprintf(fid,'--- Environment ------------------------------------------------------------------------------------------------------------------ \n');
fprintf(fid,'%g				! RHO 			! KG/M**3 	! Fluid specific volume \n',data.env.rho);
fprintf(fid,'%g				! G			! M/S**2	! Gravity \n',data.env.g);
fprintf(fid,'%g                 ! DEPTH			! M		! Water depth\n',data.env.depth);
fprintf(fid,'0	0              ! XEFF YEFF		! M		! Wave measurement point\n');
fprintf(fid,'--- Description of floating bodies -----------------------------------------------------------------------------------------------\n');
fprintf(fid,'%g				! Number of bodies\n',nBodies);
for c=1:nBodies
    fprintf(fid,'--- Body %g -----------------------------------------------------------------------------------------------------------------------\n',c);
    fprintf(fid,'%s_Nemoh.dat ! Name of mesh file\n',BodyNames{c});
    fprintf(fid,'%g %g			! Number of points and number of panels 	\n',nodes(c),panels(c));
    fprintf(fid,'6				! Number of degrees of freedom\n');
    fprintf(fid,'1 1 0 0 0 0 0		! Surge\n');
    fprintf(fid,'1 0 1 0 0 0 0		! Sway\n');
    fprintf(fid,'1 0 0 1 0 0 0		! Heave\n');
    fprintf(fid,'2 1 0 0 %g %g %g		! Roll about a point\n',RefP(c,:));
    fprintf(fid,'2 0 1 0 %g %g %g		! Pitch about a point\n',RefP(c,:));
    fprintf(fid,'2 0 0 1 %g %g %g		! Yaw about a point\n',RefP(c,:));
    fprintf(fid,'6				! Number of resulting generalised forces\n');
    fprintf(fid,'1 1 0 0 0 0 0		! Force in x direction\n');
    fprintf(fid,'1 0 1 0 0 0 0		! Force in y direction\n');
    fprintf(fid,'1 0 0 1 0 0 0		! Force in z direction\n');
    fprintf(fid,'2 1 0 0 %g %g %g		! Moment force in x direction about a point\n',RefP(c,:));
    fprintf(fid,'2 0 1 0 %g %g %g		! Moment force in y direction about a point\n',RefP(c,:));
    fprintf(fid,'2 0 0 1 %g %g %g		! Moment force in z direction about a point\n',RefP(c,:));
    fprintf(fid,'0				! Number of lines of additional information \n');
end
fprintf(fid,'--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'%g %g	%g %g	! Output Type, Number of wave frequencies, Min, and Max (rad/s)\n',[1 data.nemoh.frequencies]);
fprintf(fid,'%g	%g %g		! Number of wave directions, Min and Max (degrees)\n',data.nemoh.wavedir);
fprintf(fid,'--- Post processing ---------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'0 0 0		! IRF 				! IRF calculation (0 for no calculation), time step and duration\n');
fprintf(fid,'0				! Show pressure\n');
fprintf(fid,'0 0 0		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n');
fprintf(fid,'0 0 0 0	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction\n');	
fprintf(fid,'0					! Response Amplitude Operator (RAO), 0 no calculation, 1 calculated\n');
fprintf(fid,'1					! output freq type, 1,2,3=[rad/s,Hz,s]\n');
fprintf(fid,'--- QTF----------------------------------------------------------------------------------------------------------------------------\n');
fprintf(fid,'0         				! QTF flag, 1 is calculated\n');
fclose(fid);
%% Manage Nemoh folders
system(['mkdir ',folderResults,filesep,'mesh']);
system(['mkdir ',folderResults,filesep,'results']);
for ii=1:nBodies
    movefile([folderMesh,filesep,char(BodyNames{ii}),'_Nemoh.dat'],folderResults)
end
movefile('Nemoh.cal',folderResults)
copyfile('input_solver.txt',folderResults)
%% Run Nemoh
Nemoh(['.',filesep,folderResults],data.flags.plot_mesh_Nemoh); 
end
%% FUNCTIONS
function [Rz] = Rz(psi)

Rz=[cos(psi)  -sin(psi)  0;
    sin(psi)  cos(psi)   0;
    0          0      1];

end