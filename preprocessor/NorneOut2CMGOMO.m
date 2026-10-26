
mrstModule add deckformat;

% Send to output folder NorneforCMG

outpth = fullfile('C:\Users\cigweo1\Desktop\mrst-2024a_Geomech_v3\examples\data\NorneforCMGOMO', 'NorneforCMGOMO') ;

% outpth = fullfile('C:\Users\cigweo1\Desktop\mrst-2024a_Geomech_v3\examples\data') ;

%% read deck for Norne

grdecl = fullfile(getDatasetPath('norne'), 'NORNE.GRDECL');

grdecl = readGRDECL(grdecl);
fullCartDims = grdecl.cartDims;
usys   = getUnitSystem('METRIC');
grdecl = convertInputUnits(grdecl, usys);

grdecl = cutGrdecl(grdecl, [10 25; 35 55; 1 22]); 
grdecl.ACTNUM = ones(size(grdecl.ACTNUM));

G = processGRDECL(grdecl);
G = G(1);
G = computeGeometry(G);

% figure;
% plotGrid(G);

%% Write out data for CMG
fid = fopen([outpth 'ZCORN0.GRID'], 'w'); 
fprintf(fid, '*ZCORN ** %d float\n',  numel(grdecl.ZCORN));
fprintf(fid, '%0.8f %0.8f %0.8f %0.8f %0.8f %0.8f\n', grdecl.ZCORN(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'COORD0.GRID'], 'w'); 
fprintf(fid, '*COORD ** %d float\n', numel(grdecl.COORD));
fprintf(fid, '%0.8f %0.8f %0.8f %0.8f %0.8f %0.8f\n', grdecl.COORD(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'ACTNUM.GRID'], 'w'); 
fprintf(fid, '*NULL *ALL ** %d float\n', numel(grdecl.ACTNUM));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.ACTNUM(:));fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'PERMI0.txt'], 'w'); 
fprintf(fid, '*PERMI *ALL ** %d float\n', numel(grdecl.PERMX));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.PERMX(:)*1.013e15);fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'PERMJ0.txt'], 'w'); 
fprintf(fid, '*PERMI *ALL ** %d float\n', numel(grdecl.PERMX));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.PERMY(:)*1.013e15);fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'PERMK0.txt'], 'w'); 
fprintf(fid, '*PERMI *ALL ** %d float\n', numel(grdecl.PERMX));
fprintf(fid, '%d %d %d %d %d %d\n', grdecl.PERMZ(:)*1.013e15);fprintf(fid, '\n');
fclose(fid);

fid = fopen([outpth 'POR0.txt'], 'w'); 
fprintf(fid, '*POR *ALL ** %d float\n', numel(grdecl.PORO));
fprintf(fid, '%d %d %d %d %d %d\n', max(grdecl.PORO(:), 0.1));fprintf(fid, '\n');
fclose(fid);










