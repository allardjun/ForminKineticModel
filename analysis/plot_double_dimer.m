%% 
% graphs probability of occlusion for double and n-dimerized fh1 of all
% lengths and binding site locations
% graphs double minus dimerized probability of occlusion

% higher probability of occlusion corresponds to slower
% polymerization rate

%% initialization
err = 0;

path = '../../PolymerData/';
m2 = dlmread(append(path,'double_200.txt'));
m3 = dlmread(append(path,'dimer_122.txt'));

%% double
%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
NFil = 2;
p_occ2a_all = [];
p_occ2b_all = [];
length_array2 = 1:200;
length_vec2 = [];
iSite_vec2 = [];
err = 0;

%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
for iN=length_array2
    N_All =NFil*iN;
    for iy =1:iN
        p_a = m2(iN, 16 + 2*(N_All +1) + 7*(iy - 1));
        p_b = m2(iN, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil));
        if p_a > 1 || p_b >1
            err = err+1;
        end
        length_vec2 = [length_vec2 iN]; %list of filament lengths
        iSite_vec2 = [iSite_vec2 iy]; %all iSites used
        p_occ2a_all = [p_occ2a_all p_a];
        p_occ2b_all = [p_occ2b_all p_b];
    end
end

X2a = [length_vec2; iSite_vec2; p_occ2a_all]' ;
F2a = scatteredInterpolant(X2a(:,1),X2a(:,2),X2a(:,3),'linear','nearest');

X2b = [length_vec2; iSite_vec2; p_occ2b_all]';
F2b = scatteredInterpolant(X2b(:,1),X2b(:,2),X2b(:,3),'linear','nearest');

%% n-terminal dimerized

NFil=2;
p_occ3a_all = [];
p_occ3b_all = [];
length_array3 = 1:122;
length_vec3 = [];
iSite_vec3 = [];
err = 0;
err2 = 0;

%POcc1(iN, isite)=[occlusion probability of the iNth N value at isite]
for iN=length_array3
    N_All =NFil*iN;
    for iy =1:iN
        length_vec3 = [length_vec3 iN]; %list of filament lengths
        iSite_vec3 = [iSite_vec3 iy]; %all iSites used in percent of total filament length
        p_occ3a_all = [p_occ3a_all m3(iN, 16 + 2*(N_All +1) + 7*(iy - 1))];
        p_occ3b_all = [p_occ3b_all m3(iN, 16 + 2*(N_All +1) + 7*(iy - 1) + (6 + 9*iN + 2 + NFil + NFil))];
    end
end


X3a = [length_vec3' iSite_vec3' p_occ3a_all'];
X3b = [length_vec3' iSite_vec3' p_occ3b_all'];

F3a = scatteredInterpolant(X3a(:,1),X3a(:,2),X3a(:,3),'linear','nearest');
F3b = scatteredInterpolant(X3b(:,1),X3b(:,2),X3b(:,3),'linear','nearest');

%% scatter plots
% graphs scatter plot of p_occlusion for all lengths and iSites
% graphs comparison of p_occDouble - p_occDimer

%double
figure()
plot3(X2a(:,1),X2a(:,2),X2a(:,3),'.')
xlabel('length');
ylabel('iSite');
zlabel('p-occ');
title('double probability of occlusion');

%dimer
figure()
plot3(X3a(:,1),X3a(:,2),X3a(:,3),'.')
xlabel('length');
ylabel('iSite');
zlabel('p-occ');
title('dimerized probability of occlusion');

% comparison
szX3a = size(X3a);
n= szX3a(1);
X_comp = X2a(1:n,3) - X3a(:,3);
X_compa = []; X_compb = [];
for m = 1:n
    if X_comp(m) < 0
        X_compa = [X_compa; X3a(m,1) X3a(m,2) X_comp(m)];
    else
        X_compb = [X_compb; X3a(m,1),X3a(m,2),X_comp(m)];
    end
end

figure()

plot3(X_compa(:,1),X_compa(:,2),X_compa(:,3),'.','Color','b')
hold on
plot3(X_compb(:,1),X_compb(:,2),X_compb(:,3),'.','Color','r')
hold off

% X_compa has iSites where probability of occlusion is higher for double
% than NTD, i.e. actin polymerization rate is faster for NTD
% X_compb has iSites where p_occ is higher for NDT than double

legend('dimer faster','double faster'); 
xlabel('length')
ylabel('iSite')
zlabel('p_occ double - dimer')






