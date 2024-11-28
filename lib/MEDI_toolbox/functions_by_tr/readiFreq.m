% filename = 'iField.bin';
% matrix_size = [288,384,96];
% fid = fopen(filename,'r');
% a = fread(fid,inf,'float');
% a= reshape(a,matrix_size);
% fclose(fid);

fid = fopen('a.bin','w');
fwrite(fid,a,'complex');
fclose(fid);

fid = fopen('a.bin','r');
b= fread(fid,inf);
fclose(fid);