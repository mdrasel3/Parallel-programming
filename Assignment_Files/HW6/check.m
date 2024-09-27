% test of matrix multiplication
% check.m

fid = fopen('C.txt','r'); % change filename if necessary
m = fscanf(fid,'%d\n',1);
k = fscanf(fid,'%d\n',1);
C = uint8(transpose(reshape(fscanf(fid,'%f\n',k*m),k,m)));
fclose(fid)
image(C);
colormap(gray(256));
set(gca,'DataAspectRatio',[1 1 1]);