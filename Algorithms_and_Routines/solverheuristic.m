function [pk, diff, dist,asub] = solverheuristic(Stil,vj,pk, distold, tol)
%% Heuristic for solve min (Stilx)^Ts_j; e^Tx=1 x\geq 0 \norm{Stilx}^2 \leq \norm{p_k}^2
%
% Syntax: 
%       [pk, diff, dist,asub] = solverheuristic(Stil,vj,pk, distold, tol)
%
% Input: 
%         Stil: is a matrix with the L elements of S that are active in the combination convex of pk (the current iterate)
%          vj: pivot for pk
%          pk: the current iterate
%         distold: norm(vj-pk)
%         tol: tolerance for epsilon solution
%
% Output: 
%         pk: new pk that satisfies the epsilon property
%         diff: pk-vj
%         dist: \norm{pk-vj}
%         asub: coeficients (Aa=new pk)

    maxit = 1e+6;
    c = Stil'*vj; 
    [~,idxq] = min(c);
    xl=zeros(size(Stil,2),1);
    xl(idxq)=1;
    npkold = pk'*pk;  % r^2
    dist = distold;
    delta = 100; %1/sqrt(npkold); % 1/npkold;
    pk = Stil*xl;
    npk= pk'*pk;
    rquad = npk-npkold;
    S= Stil'*Stil;
    S= (S+ S')/2;
    optm = 0;
    kin = 0;
    eta = 1/sqrt(npkold); % 1/npkold;
    ee = ones(size(Stil,2),1)';
    lb = zeros(size(Stil,2),1);
    x = xl;
       
    while abs(rquad) > tol || optm > tol % rquad > npkold
          Q = S; %eta*S; 
          opts = optimoptions('quadprog','OptimalityTolerance',tol/eta);
          [x,~,~,OUTPUT,~] = quadprog(Q,c/eta,[],[],ee,1,lb,[],[],opts);
          pk = Stil*x;
          optm = eta*OUTPUT.firstorderopt;
          asub = x;     
       npk = pk'*pk;
       if abs(npk-npkold) > 0.8*abs(rquad)
           delta = 10*delta;
       end
       rquad = npk-npkold; 
       eta =  max(eta + 0.5*delta*rquad,0);
       kin = kin + 1;
    end
    
    diff = vj-pk;
    dist = norm(pk-vj,2);
end
