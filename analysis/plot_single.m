%creates gradient map for P_Occlusion for single filament


%initialization

clear all

path = '../../PolymerData/';
A=dlmread(append(path,'single1_300.txt'));
NFil=1;
pOcc = [];
N_Array = 1:300;
N_vec = [];
iSite_vec = [];

%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
for iN=1:300
    N = N_Array(iN);
    N_All =NFil*N;
    for iy =1:N
        N_vec = [N_vec N]; %correct length list of filament lengths
        iSite_vec = [iSite_vec iy]; %all iSites used
        pOcc = [pOcc A(iN, 16 + 2*(N_All +1) + 7*(iy - 1))];
    end
end    

figure()
X(:,1) = transpose(N_vec);
X(:,2) = transpose(iSite_vec);
X(:,3) = transpose(pOcc);
disp(X)
scatter3(X(:,1),X(:,2),X(:,3),'.')

% plot(N_vec, pOcc)
% % 
% u=linspace(min(X(:,1)),max(X(:,1)),300);
% v=linspace(min(X(:,2)),max(X(:,2)),300);
% [U,V]=meshgrid(u,v);
% % F=scatteredInterpolant(X(:,1),X(:,2),X(:,3),'linear','none');
% contourf(U,V,X,100,'LineColor','none') 
% % c = colorbar;
% % caxis([0,1]);
% % c.Label.String = 'Probability of Occlusion' ;
% % xlabel('Total Length of FH1 domain')
% % ylabel('Binding Site Location Along Filament')