function [Decision,pk,a,iterations, Matrixinformation]= SPHERICALTAPLUSHEU(data_mat,p,tol,gamma,heur_on,extrainformation)
%% Triangle Algorithm for spherical problem with option to have heuristic
%
% Syntax: 
%       [Decision,pk,a,iterations, Matrixinformation]= SPHERICALTAPLUSHEU(data_mat,p,tol,gamma,heur_on,extrainformation)
%
% Input: 
%         data_mat: set of n points in R^m (with or not norm(A(:,j))=1 for
%         all j)
%          p: query point
%         tol: stop criterion to find a p_epsilon-solution
%         gamma: constant for epsilon/M property
%        heur_on = 1 with subproblem otherwise the original SPHERICAL TA 
% Output: 
%         Decision: 1 when p \in conv(A)  
%         Decision: 0 when p \notin conv(A)
%         pk: witness or point that satisfies the epsilon solution
%         a: coeficients (Aa=pk)
%         iterations: number of external iterations
%         Matrixinformation: [iteration that the epsilon property is not satisfied,count of iterations that the epsilon property is not satisfied, number of times it satisfies the tol-property, \norm{p_t}, \norm{p_k}, \norm{p_t}- \norm{p_k}  1+gamma*tol iterations]; 

%% Initialization
    Matrixinformation = []; 
    countsatisfazepsilonsolution = 0;
    knaosatisfazepsilon = 0;
    A =[];
    for i= 1:size(data_mat,2)
         v = (data_mat(:,i)-p)/norm((data_mat(:,i)-p));
         A = [A v];
    end

    [~, n] = size(A);
    i=1;
    iterations = 0;
    pk = A(:,i);
    a = zeros(n,1); 
    a(i)= 1;
    normpk = norm(pk,2);
    if nargin<=5
       extrainformation = 0; 
        if nargin <=4
           heur_on = 0;
        end
    end
%% loop principal    
while (normpk > tol)               %while pk it's not a p_epsilon-solution
    prod = A'*pk;
    [~, idxj] = min(prod);       %Choose vj \in \argmin {vi^T(pk - p): v_i \in A/{pk}}
    j = idxj(1);
    vj = A(:,j);
    normpk2 = norm(pk,2)^2; 
    diff = vj - pk;
    dist = norm(diff, 2);  
    
    if 2*vj'*pk <= normpk2           % Simple pivot characterization     
        if (dist^2<1+gamma*tol) && (heur_on==1)            % eps/M - property does not hold
            [idxa,~] = find(abs(a) > 1e-8);% | [1:size(a,1)]'== j); 
            Atil = A(:,idxa);  
            distold = dist; 
            [pk, diff, dist, asub] = solverheuristic(Atil,vj,pk,dist,tol);
            if extrainformation == 1
               knaosatisfazepsilon = knaosatisfazepsilon +1; 
               diferencadists = dist-distold;
               if dist^2>=1+gamma*tol
                  countsatisfazepsilonsolution = countsatisfazepsilonsolution+1;                 
               end
               Matrixinformation = [Matrixinformation; iterations knaosatisfazepsilon countsatisfazepsilonsolution dist distold diferencadists 1+gamma*tol]; 
            end
            a = zeros(n,1);
            a(idxa) = asub;
        end 

        alpha = min(1,(pk'*(-diff))/(dist^2));  %Linear search
        pk = (1-alpha)*pk + alpha*vj;            %new iteration
        normpk = norm(pk,2);
   
       if (iterations > 1/(tol^2))
           fprintf('Maximum number of iterations reached!\n')
           Decision = -1;
           return
       end
        
        if (alpha == 1)
            a = 0*a;
            a(j) = 1;
        else
            a = (1-alpha)*a;
            a(j) = a(j) + alpha;
        end 
            
    else 
        v = -pk;                               
        normv = norm(v)^2/2;
        Decision = 0;
        % fprintf('Spherical TA: The point p is OUTSIDE of conv(A)\n')
        return
    end
    
    iterations=iterations+1;
end

Decision = 1;
%        fprintf('Spherical TA: The point p is INSIDE of conv(A)\n')
return;
