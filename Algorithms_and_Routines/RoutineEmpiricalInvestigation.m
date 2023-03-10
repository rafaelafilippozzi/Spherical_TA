% Routine for Empirical investigation of iteration complexity for the Spherical CHMP with p in convS
%Choice one
Extrainformation = 1; %On Extra information;
%Extrainformation = 0; %Off Extra information;

rng(11)
Num_rep=2;
Dim_array=[50];
Vtx_array=[500];
Eps_array=[0.00001,0.000025,0.00005,0.0001,0.00025,0.0005, 0.001];
time_spta = zeros(Num_rep,length(Eps_array));
time_sptawithheu = zeros(Num_rep,length(Eps_array));
time_spasfw = zeros(Num_rep,length(Eps_array));
time_spasfwwithheu = zeros(Num_rep,length(Eps_array));
aaMatrixinformation = [];
aaMatrixinformationasfw = [];
vDecisionsphe = [];
vDecisionheu = [];
vDecisionspheasfw = [];
vDecisionheuasfw = [];
vcountas = [];
vcountasheu = [];
iterationssemheu =  zeros(Num_rep,length(Eps_array));
iterationsheusubprob =  zeros(Num_rep,length(Eps_array));
iterationssemheuasfw =  zeros(Num_rep,length(Eps_array));
iterationsheusubprobasfw =  zeros(Num_rep,length(Eps_array));
miterations = zeros(size(Eps_array,2),1); 
miterationssemheu = zeros(size(Eps_array,2),1);
miterationsasfw = zeros(size(Eps_array,2),1); 
miterationssemheuasfw = zeros(size(Eps_array,2),1); 

Num_dim=Dim_array(1);
Num_vtx=Vtx_array(1);
    
    for kk=1:Num_rep
        data_mat = generateArandom(Num_dim, Num_vtx );     %each element of A={v1,...,vn} randomly generated according to a uniform distribution in the unit ball of Rm
       %Pdentro Fácil CASO A
%        p = max(2*ones(Num_dim,1)+ randn(Num_dim,1),0);
%        p = p/(Num_vtx*norm(p));              %p with a lot of chance of being inside the relative interior of the convA

      %Casos B,C,D
        aux = sum(data_mat,1);    
        [~, idx] = sort(aux,'descend');
        qidx = idx(1);
        lidx = idx(2);   

%        %Pdentro Difícil CASO B
          p = 0.5*data_mat(:,qidx) + 0.5*data_mat(:,lidx);
%        %Pfora Fácil CASO C
%          p = 1.5*(data_mat(:,qidx)+ data_mat(:,lidx))/2;
%        %Pfora Dificil CASO D
%        p = 1.01*(data_mat(:,qidx)+ data_mat(:,lidx))/2;        
        for jj=1:length(Eps_array)
            epsilon = Eps_array(jj);
            jj
            kk

            disp('spta nosso')
            tic; 
            [Decision,pk,a,iterations,~]= SPHERICALTAPLUSHEU(data_mat,p,epsilon,1);
            scale_spta_nosso=toc;
            vDecisionsphe = [vDecisionsphe; Decision];
            iterationssemheu(kk,jj) =  iterations;
            time_spta(kk,jj)= time_spta(kk,jj)+scale_spta_nosso;
            disp('end spta nosso')         
            
            disp('spta with heu')
            tic; 
           [Decisionheu,pkheu,a,iterationsheusubpro,aMatrixinformation]= SPHERICALTAPLUSHEU(data_mat,p,epsilon,1,1,1);           
            tasphericalnosso_end = toc;
            vDecisionheu = [vDecisionheu; Decisionheu];
            iterationsheusubprob(kk,jj) = iterationsheusubpro;
            if Extrainformation == 1
                aaMatrixinformation = [aaMatrixinformation; aMatrixinformation iterationsheusubpro*ones(size(aMatrixinformation,1),1)];
            end
            time_sptawithheu(kk,jj)= time_sptawithheu(kk,jj)+tasphericalnosso_end;
            disp('end spta with heu')
        end
    end     
miterationssemheu = mean(iterationssemheu); 
miterations = mean(iterationsheusubprob);
miterationssemheuasfw = mean(iterationssemheuasfw); 
miterationsasfw = mean(iterationsheusubprobasfw);

names=strings(size(Eps_array,2),1);
for i=1:length((Eps_array(1,:)))
    this_size=[num2str(Eps_array(i))];
    names(i)=this_size;
end
umsobepsilon = [];
umsobepsilon2 = [];
for i = 1:length((Eps_array(1,:)))
    umsobepsilon = [umsobepsilon 1/Eps_array(i)];
    umsobepsilon2 = [umsobepsilon2 1/(Eps_array(i)^2)];
end

    this_title=['Relação entre Iteração e Tolerância ',num2str(Dim_array(1)), 'x', num2str(Vtx_array(1))];
    hold on
    title(this_title);
    plot(log10(miterationssemheu'),'k-','LineWidth',1.5)
    plot(log10(miterations'),'k--','LineWidth',1.5)
    plot(log10(umsobepsilon'),'k:','LineWidth',1.5)
    plot(log10(umsobepsilon2'),'k-.','LineWidth',1.5)
    leg2 = legend('TA-Esférico','TA-EsféricoH','$1/\varepsilon$','$1/{\varepsilon^2}$','Interpreter','latex');
    set(gca,'xtick',1:length((Eps_array(1,:))),'xticklabel',names)
    xlabel ('Tolerâncias');
    ylabel (['$\log_{10}$(itera\c{c}\~oes)'],'Interpreter','latex');
    hold off;
