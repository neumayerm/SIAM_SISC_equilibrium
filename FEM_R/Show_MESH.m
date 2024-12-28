% Fast Numerical Techniques for Inverse Problems with Underlying Equilibrium Systems
% 
% Presentation of the FE mesh for the demo implementation
%  
%  Inputs
%   + Node         Node coordinates   
%   + Inzid        Inzidenz matrix 
%   + Material     Marker assigned to the finite elements: 
%       o  0  ROI 
%       o  50 Pipe 
%       o 150 Backside     
%   + BndMrk       Marker assigned to nodes
%       o 1:8 Electrodes  
%       o 100 Shield
%   + NrElec   
%
% EMS 2022
% Contact: neumayer@tugraz.at
clear all, close all, clc


% Show these electrodes
IDELEC = [1 2 3];

load MESH   % Data for sensor with 8 electrodes


% Plot Mesh
figure; hold on, set(gcf,'Color','White'), set(gca,'FontSize',16);
    for ii = 1:size(MESH.Inzid,1)
    
         inz = MESH.Inzid(ii,:);

        T = MESH.Node(inz,:);
        Te = [T;T(1,:)];

        if MESH.Material(ii) == 0  % Interior         
         h1=   fill(T(:,1),T(:,2),'w','Edgecolor','k');            
        end        
        if MESH.Material(ii) == 50   % Pipe 
         h2= fill(T(:,1),T(:,2),[1 1 1]*0.6);            
        end
        if MESH.Material(ii) == 150   % Pipe         
         h3=   fill(T(:,1),T(:,2),[1 1 1]*0.8);                     
        end        
    end   
    axis off 
    axis square

% Create legend
cnt = 1;
LEG{cnt} = 'ROI:  MESH.Material =  0'; cnt = cnt+1;
LEG{cnt} = 'Pipe: MESH.Material = 50'; cnt = cnt+1;
LEG{cnt} = 'Backside: MESH.Material = 150'; cnt = cnt+1;
LEG{cnt} = 'Shield: MESH.BndMrk = 100'; cnt = cnt+1;    
    
        
% Mark the shield
id = find(MESH.BndMrk==100); hS = plot(MESH.Node(id,1),MESH.Node(id,2),'+r','LineWidth',2);

% Mark the electrodes from IDELEC
hELEC = [];
for ii = 1:length(IDELEC)
    id = find(MESH.BndMrk==IDELEC(ii)); h4 = plot(MESH.Node(id,1),MESH.Node(id,2),'o','LineWidth',2);
    LEG{cnt} = sprintf('Electrode %i, MESH.BndMrk = %i',IDELEC(ii),IDELEC(ii)); cnt = cnt+1;
    hELEC = [hELEC,h4];
end

legend([h1,h2,h3,hS,hELEC],LEG,'Location','northeastoutside')

