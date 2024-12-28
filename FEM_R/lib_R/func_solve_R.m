function [FEM] = func_solve_R(FEM)


diagx        = sparse(1:FEM.N,1:FEM.N,FEM.x);
Khat         = FEM.Kinihat + FEM.A1hat*diagx*FEM.A1hat' + FEM.A2hat*diagx*FEM.A2hat';
C            = chol(Khat(FEM.perm,FEM.perm));
FEM.Gtilde   =  C \ (C'\FEM.Ftilde(FEM.perm,:));
FEM.Gtilde   = FEM.Gtilde(FEM.iperm,:);     
FEM.Y        = FEM.Gtilde'*FEM.F;        
FEM.y        = FEM.Y(FEM.idxM);
   