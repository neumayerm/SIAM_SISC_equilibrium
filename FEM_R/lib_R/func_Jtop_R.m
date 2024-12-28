function [Jtr] = func_Jtop_R(FEM,r)


    R = zeros(size(FEM.Y));
    R(FEM.idxM) = r;
    R = sparse(R);
   
    GtA1=(FEM.Gtilde'*FEM.A1hat);
    GtA2=(FEM.Gtilde'*FEM.A2hat);
    
    Jtr = sum((GtA1.*(R*GtA1) + GtA2.*(R*GtA2)));
    
    Jtr = Jtr(:);
    
