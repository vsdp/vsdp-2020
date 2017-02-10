function [infeas,x,y] = vsdpinfeas(A,b,c,K,choose,x0,y0,~,opts)
% VSDPINFEAS Infeasibility check for semidefinite-quadratic-linear programming.
%
%   [isinfeas,X,Y] = VSDPINFEAS(A,b,c,K,choose,x0,y0)
%      The block-diagonal format (A,b,c,K) is explained in 'mysdps.m'.
%
%         'choose'    If the character is 'p', primal infeasibility should be
%                     verified.  If 'd', dual infeasibility should be verified.
%
%         'x0'        If 'choose == d', provide an approximate solution
%                     computed by 'mysdps'.  Otherwise, this value is ignored.
%
%         'y0'        If 'choose == p', provide an approximate solution
%                     computed by 'mysdps'.  Otherwise, this value is ignored.
%
%      The output is:
%
%         'isinfeas'  Returns 1 if the primal or -1 if the dual problem is
%                     proved to be infeasible and 0 if infeasibility cannot be
%                     verified.
%
%         'x'         Contains a rigorous certificate (improving ray) of dual
%                     infeasibility, if it is not equal to NaN.
%
%         'y'         Is a rigorous certificate (improving ray) of primal
%                     infeasibility, if it is not equal to NaN.
%
%   VSDPINFEAS(...,z0,opts) optionally provide an approximate solution 'z0',
%      that is ignored and an option structure, explained in 'vsdpinit'.
%
%   Example:
%
%       %  min <[0 0; 0 0], X>
%       % s.t. <[1 0; 0 0], X> = [e];
%       %      <[0 1; 1 d], X> = [1];
%       %                   X >= 0
%
%       e = -0.01;  % Infeasible, because X(1,1) <= 0!
%       d = 0.1;
%       A1 = [1 0; 0 0];
%       A2 = [0 1; 1 d];
%       b = [e; 1];
%       c = zeros(4,1);
%       K.s = 2;
%       A = [A1(:), A2(:)];  % vectorize
%       [objt,x0,y0,z0,info] = mysdps(A,b,c,K);
%       choose = 'p';
%       [isinfeas, x, y] = vsdpinfeas(A,b,c,K,choose,x0,y0);
%
%   See also mysdps, vsdpinit.

% Copyright 2004-2012 Christian Jansson (jansson@tuhh.de)

narginchk(7,9);
if (nargin < 9)
  opts = [];
end

VSDP_OPTIONS = vsdpinit(opts);

[A,Arad,b,brad,c,crad,K,x0,y0,~,IF] = import_vsdp(A,b,c,K,x0,y0,[]);

% Save rounding mode
rnd = getround();
setround(0);

% Default output
infeas = 0;
x = NaN;
y = NaN;

for i = 1 % Allows break in switch statement
  switch (choose)
    
    case 'p' % Verify primal infeasibility
      
      if (isempty(y0) || any(isnan(y0)))
        error('VSDP:VSDPINFEAS', 'inapplicable approximation y0 provided');
      end
      
      % Step 1: Rigorous enclosure for z = -A*y
      if (K.f > 0)
        % Solve dual linear constraints rigorously
        y0 = vuls([],[],struct('mid',A(1:K.f,:),'rad',Arad(1:K.f,:)),...
          zeros(K.f,1),[],[],y0,[],opts);
        if ~isstruct(y0)
          break; % cannot verify infeasibility
        else
          y0rad = y0.rad;
          y0 = y0.mid;
        end
      else
        y0rad = sparse(size(y0,1),1);
      end
      % default rounding mode for verification
      setround(1);
      % z = -A*y  (with free variables)
      [z,zrad] = resmidrad(0,0,A,y0,0,0,Arad,y0rad);
      
      % Step 2: Compute positive cost
      alpha = prodsup(-b',y0,brad',y0rad);
      if (alpha >= 0)
        break; % cannot verify infeasibility
      end
      
      % Step 3: Verified lower bounds on cone eigenvalues
      dl = inf;
      % Bound for linear variables
      if (K.l > 0)
        ind = K.f+1:K.f+K.l;
        dl = -max(zrad(ind) - z(ind));
      end
      % Eigenvalue bound for second-order cone variables
      ind = K.f + K.l;
      for j = 1:length(K.q)
        ind = ind(end)+2:ind(end)+K.q(j);
        zj1 = zrad(ind(1)-1) - z(ind(1)-1);  % -inf(zq(1))
        zjn = abs(z(ind)) + zrad(ind);  % sup(abs(zq(2:end)))
        zjn = sqrtsup(zjn'*zjn);  % sup(||zq(2:end)||)
        dl = min(dl,-(zj1+zjn));  % inf(zq(1)-||zq(2:end)||)
      end
      if (dl < 0)
        break; % cannot verify infeasibility
      end
      % Eigenvalue bound for semidefinite cone
      blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
      for j = length(K.s):-1:1
        blks = blke - K.s(j)*(K.s(j)+1)/2 + 1;
        lambdaMin = bnd4sd(struct('mid',z(blks:blke),...
          'rad',zrad(blks:blke)),1,VSDP_OPTIONS.FULL_EIGS_ENCLOSURE);
        dl = min(dl,lambdaMin);
        if (dl < 0)
          break; % cannot verify infeasibility
        end
        blke = blks - 1;
      end
      
      % Step 4: Final feasibility check, compute result
      if (dl >= 0)
        infeas = 1;
        if any(y0rad)
          y = midrad(full(y0),full(y0rad));
        else
          y = y0;
        end
      end
      
      
    case 'd' % Verify dual infeasibility
      
      if (isempty(x0) || any(isnan(x0)))
        error('VSDP:VSDPINFEAS', 'inapplicable approximation x0 provided');
      end
      
      % Step 1: Check c'*x > 0
      [alpha, alpharad] = spdotK(c,x0,3);
      setround(1);  % default rounding mode for verification
      alpharad = alpharad + crad'*x0;
      if (alpharad >= -alpha)
        break; % cannot verify infeasibility
      end
      
      % Step 2: Inclusion for x such that A*x = 0
      vx = vuls([],[],struct('mid',A,'rad',Arad),0*b,[],[],x0,[],opts);
      if ~isstruct(vx)
        break; % cannot verify infeasibility
      elseif (prodsup(vx.mid',c,full(vx.rad)',crad) >= 0)  % if c'*x < 0
        A = cat(2,A,c);
        Arad = cat(2,Arad,crad);
        % find inclusion such that c'*x = alpha
        vx = vuls([],[],struct('mid',A,'rad',Arad),...
          struct('mid',[b; alpha],'rad',[brad; 0]),[],[],x0,[],opts);
        if ~isstruct(vx)
          break; % cannot verify infeasibility
        end
      end
      x0rad = sparse(vx.rad);
      x0 = full(vx.mid);
      
      % Step 3: Verified lower bounds on cone eigenvalues
      lb = inf;
      % Bound for linear variables
      if (K.l > 0)
        ind = K.f+1:K.f+K.l;
        lb = -max(x0rad(ind) - x0(ind));
      end
      % Eigenvalue bound for second-order cone variables
      ind = K.f + K.l;
      for j = 1:length(K.q)
        ind = ind(end)+2:ind(end)+K.q(j);
        xj1 = x0rad(ind(1)-1) - x0(ind(1)-1);  % -inf(xq(1))
        xjn = abs(x0(ind)) + x0rad(ind);  % sup(abs(xq(2:end)))
        xjn = sqrtsup(xjn'*xjn);  % sup(||xq(2:end)||)
        lb = min(lb,-(xj1+xjn));  % inf(xq(1)-||xq(2:end)||)
      end
      if (lb < 0)
        break; % cannot verify infeasibility
      end
      % Eigenvalue bound for semidefinite cone
      blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
      for j = length(K.s):-1:1
        dind = cumsum(1:K.s(j));  % index for diagonal entries
        blk3 = dind(end);
        vx = struct('mid',x0(blke-blk3+1:blke) / 2, ...
          'rad',x0rad(blke-blk3+1:blke) / 2);
        vx.mid(dind) = 2 * vx.mid(dind);  % regard mu=2 for x
        vx.rad = vx.rad + vx.rad.*sparse(dind,1,1,blk3,1);
        lambdaMin = bnd4sd(vx,1,VSDP_OPTIONS.FULL_EIGS_ENCLOSURE);
        lb = min(lb,lambdaMin);
        if (lb < 0)
          break; % cannot verify infeasibility
        end
        blke = blke - blk3;
      end
      
      % Step 4: Final feasibility check, compute result
      if (lb >= 0)
        infeas = -1;
        x = export_vsdp(IF,K,midrad(full(x0),full(x0rad)));
      end
      
    otherwise
      error('VSDP:VSDPINFEAS', '"choose" must be ''p'' or ''d''');
  end
end

if (VSDP_OPTIONS.VERBOSE_OUTPUT && (infeas == 0))
  str = {'primal', 'dual'};
  str = str{(choose == 'd') + 1};
  fprintf(1, 'VSDPINFEAS: could not verify %s infeasibility\n', str);
end

setround(rnd);  % restore rounding mode

end
