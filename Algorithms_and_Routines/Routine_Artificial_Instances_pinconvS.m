% Routine for Artificial Instances for the Spherical CHMP with p in convS
%Choice one
%Extrainformation = 1; %On Extra information;
Extrainformation = 0; %Off Extra information;

rng(111)
Num_rep=10;
Dim_array=[50,100,200,300];
Vtx_array=[500,1000,2000,3000];
Eps_array=[0.0000001];
time_spta=zeros(length(Dim_array),length(Eps_array));
time_sta_with_heu=zeros(length(Dim_array),length(Eps_array));
Matrizinformacao = [];
viterationsheusubprob=[];
viterationssemheu=[];
vDecisionheu = [];
vDecisionsphe = [];
iterrr=0;
for ii=1:length(Dim_array)
    Num_dim=Dim_array(ii);
    Num_vtx=Vtx_array(ii);
    jj = 1;
%    for jj=1:length(Eps_array)
        epsilon=Eps_array(jj);
        for kk=1:Num_rep
            ii
            jj
            kk
            iterrr=iterrr+1;
          
             
            A=Random_pts(Num_dim,Num_vtx,'unit ball');
            X=rand(Num_vtx,1);
            b=A*X;
            [M_con,N_var]=size(A);
            
            M=1200;

            tmp_mat=[A,zeros(M_con,1);ones(1,N_var),1];
            tmp_b=[-b;-M];

            data_mat=[tmp_mat,tmp_b;zeros(1,N_var+1),1];
            p=[zeros(Num_dim+1,1);1/(1+M)];

            disp('spta nosso')
            tic; 
            [Decision,pk,a,iterationssemheu]= SPHERICALTAPLUSHEU(data_mat,p,epsilon,1);
            scale_spta_nosso=toc;
            vDecisionsphe = [vDecisionsphe; Decision];
            viterationssemheu = [viterationssemheu; iterationssemheu];
            time_spta(ii,jj)= time_spta(ii,jj)+scale_spta_nosso;
            disp('end spta nosso')                         

            disp('Spherical with Heuristic')
            tic;
            [Decisionheu,pkheu,a,iterationsheusubpro,info]= SPHERICALTAPLUSHEU(data_mat,p,epsilon,1,1,1);
            tasphericalnosso_end = toc;
            vDecisionheu = [vDecisionheu; Decisionheu];
            viterationsheusubprob = [viterationsheusubprob; iterationsheusubpro];
            if Extrainformation == 1
                Matrizinformacao = [Matrizinformacao; info];
            end
            time_sta_with_heu(ii,jj)=time_sta_with_heu(ii,jj)+tasphericalnosso_end;
            disp('end Spherical with Heuristic')
            
        end      
%    end     
end

miterationsspta = [];
miterationssptaheu = [];
for i = 1:(length(Dim_array))
    miterationsspta = [miterationsspta; mean(viterationssemheu(10*i-9:10*i))];
    miterationssptaheu = [miterationssptaheu; mean(viterationsheusubprob(10*i-9:10*i))];
end
        

names=strings(4,1);
for i=1:length(Dim_array)
    this_size=[num2str(Dim_array(i)), 'x', num2str(Vtx_array(i))];
    names(i)=this_size;
end
        
    i=1;
    this_title=['$0 \in conv(\mathcal{S})$, ','$\varepsilon = 10^{-7}$'];
    figure(i);
    hold on
    title(this_title,'Interpreter','latex');
    plot(log10(time_spta(:,i)),'k-','DisplayName','TA-Esférico','LineWidth',1.5)
    plot(log10(time_sta_with_heu(:,i)),'k--','DisplayName','TA-EsféricoH','LineWidth',1.5)
    legend('show','Location','northwest')%,'Orientation','horizontal')
    set(gca,'xtick',1:length(Dim_array),'xticklabel',names)
    xlabel ('Dimensões');
      ylabel (['$\log_{10}$(tempo em segundos)'],'Interpreter','latex');
    hold off;
