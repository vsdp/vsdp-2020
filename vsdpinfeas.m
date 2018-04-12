function [infeas,x,y] = vsdpinfeas(At,b,c,K,choose,x0,y0,~,opts)
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

[At,b,c,K,x0,y0,~,imported_fmt] = import_vsdp(At,b,c,K,x0,y0,[]);

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
        y0 = vuls([], [], At(1:K.f,:), zeros(K.f,1), [], [], y0, [], opts);
        if (~isintval (y0))
          break; % cannot verify infeasibility
        end
      end
      % default rounding mode for verification
      setround(1);
      
      % Compute rigorosly `z = -A' * y0` with free variables.
      z = -At * intval (y0);
      
      % Step 2: Compute positive cost
      if (inf_ (-b' * y0) >= 0)
        break; % cannot verify infeasibility
      end
      
      % Step 3: Verified lower bounds for each cone.
      dl = -inf;
      
      % LP cones
      if (K.l > 0)
        dl = max (dl, inf_ (z(K.f + (1:K.l))));
      end
      
      % Second-order cones
      ind = K.f + K.l;
      for j = 1:length(K.q)
        ind = (ind(end) + 2):(ind(end) + K.q(j));
        % dl = inf (z(1) - ||z(2:end)||)
        dl = min (dl, inf_ (z(ind(1)-1) - norm (intval (z(ind)))));
      end
      if (dl < 0)
        break; % cannot verify infeasibility
      end
      
      % SDP cones
      blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
      for j = length(K.s):-1:1
        blks = blke - K.s(j)*(K.s(j)+1)/2 + 1;
        lambda = veigsym (z(blks:blke));
        if (any (lambda > 0))
          break; % cannot verify infeasibility
        end
        blke = blks - 1;
      end

      % Step 4: Final feasibility check, compute result
      if (dl >= 0)
        infeas = 1;
        y = y0;
      end


    case 'd' % Verify dual infeasibility

      if (isempty(x0) || any(isnan(x0)))
        error('VSDP:VSDPINFEAS', 'inapplicable approximation x0 provided');
      end

      % Step 1: Check c'*x > 0
      setround(1);  % default rounding mode for verification
      if (inf_ (c' * x0) >= 0)
        break; % cannot verify infeasibility
      end

      % Step 2: Inclusion for x such that A*x = 0
      vx = vuls([], [], At, zeros (size (b)), [], [], x0, [], opts);
      if (~isintval (vx))
        break; % cannot verify infeasibility
      elseif ((c' * vx) >= 0)
        At = cat(2,At,c);
        % find inclusion such that c'*x = alpha
        vx = vuls([], [], At, [b; alpha], [], [], x0, [], opts);
        if (~isintval (vx))
          break; % cannot verify infeasibility
        end
      end
      x0rad = sparse(vx.rad);
      x0 = full(vx.mid);

      % Step 3: Verified lower bounds on cone eigenvalues
      lb = inf;
      
      % LP cones
      if (K.l > 0)
        ind = K.f+1:K.f+K.l;
        lb = -max(x0rad(ind) - x0(ind));
      end
      
      % Second-order cones
      ind = K.f + K.l;
      for j = 1:length(K.q)
        ind = (ind(end) + 2):(ind(end) + K.q(j));
        % dl = inf (z(1) - ||z(2:end)||)
        lb = min (lb, inf_ (z(ind(1)-1) - norm (intval (z(ind)))));
      end
      if (lb < 0)
        break; % cannot verify infeasibility
      end
      
      % SDP cones
      blke = K.f + K.l + sum(K.q) + sum(K.s.*(K.s+1))/2;
      for j = length(K.s):-1:1
        dind = cumsum(1:K.s(j));  % index for diagonal entries
        blk3 = dind(end);
        vx = x0(blke-blk3+1:blke) / 2;
        vx(dind) = 2 * vx(dind);  % regard mu=2 for x
        vxrad = x0rad(blke-blk3+1:blke) / 2;
        if (~iszero(vxrad))
          vxrad = vxrad + vxrad.*sparse(dind,1,1,blk3,1);
          vx = midrad(vx, vxrad);
        end
        lambdaMin = bnd4sd(vx, VSDP_OPTIONS.FULL_EIGS_ENCLOSURE);
        lb = min(lb,lambdaMin);
        if (lb < 0)
          break; % cannot verify infeasibility
        end
        blke = blke - blk3;
      end

      % Step 4: Final feasibility check, compute result
      if (lb >= 0)
        infeas = -1;
        x = export_vsdp(imported_fmt,K,midrad(full(x0),full(x0rad)));
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
