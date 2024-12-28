function J = func_getJ_Jtop_R(FEM)

    J = zeros(length(FEM.y),FEM.N);
    r0 = zeros(length(FEM.y),1);
    r = r0;
    
        GtA1=(FEM.Gtilde'*FEM.A1hat);
    GtA2=(FEM.Gtilde'*FEM.A2hat);
    
    for ii = 1:length(FEM.y)
        r(ii) = 1;
        
       
           R = zeros(size(FEM.Y));
           R(FEM.idxM) = r;
            R = sparse(R);

             Jtr = sum((GtA1.*(R*GtA1) + GtA2.*(R*GtA2)));
             
        J(ii,:) = Jtr';
        r = r0;
    end