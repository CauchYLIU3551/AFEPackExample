%% read the stiff matrix and mass matrix of the situation that mesh = Cube.2
%the dimension of the matrices are 26 x 26
function [A,M]=Get_Matrix(n)
  fid=fopen('stiffmatrix.txt','r');
  A=fscanf(fid,'%f',[n,n]);
  fclose(fid);

  fid2=fopen('massmatrix.txt','r');
  M=fscanf(fid2,'%f',[n,n]);
  fclose(fid2);
end