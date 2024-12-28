function J = func_getJ_Green_R(FEM)


Nelec = length(FEM.Y);

idx_row = FEM.idx_row;
idx_col = FEM.idx_col;

GtA1=(FEM.Gtilde'*FEM.A1hat); 
GtA2=(FEM.Gtilde'*FEM.A2hat);
J = (GtA1(idx_row,:).*GtA1(idx_col,:) + GtA2(idx_row,:).*GtA2(idx_col,:));
