function J = func_getJ_Jop_R(FEM)

J = zeros(length(FEM.y),FEM.N);
    x0 = zeros(FEM.N,1);
    Deltax = x0;
    
    for ii = 1:FEM.N
        %Deltax(ii)=1;
        
        %[DY,Dy] = func_Jop(FEM,Deltax);
        
           %DK = FEM.A1hat*(diagDx*FEM.A1hat') + FEM.A2hat*(diagDx*FEM.A2hat');      
           
           DK = FEM.A1hat(:,ii)*FEM.A1hat(:,ii)' + FEM.A2hat(:,ii)*FEM.A2hat(:,ii)';
            DY = (FEM.Gtilde')*(DK*FEM.Gtilde);
            Dy = DY(FEM.idxM);
        

        J(:,ii) = Dy;
        Deltax = x0;
    end