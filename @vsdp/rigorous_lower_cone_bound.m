function [lb, n, vidx, sdp_matrix] = rigorous_lower_cone_bound (obj, x, mu, weight_sdp)
% RIGOROUS_LOWER_CONE_BOUND  Helper function for rigorous lower bounds.
%
%   lb = obj.rigorous_lower_cone_bound(x, mu, weight_sdp)  Compute a rigorous
%      lower bound 'lb' of the conic variable '[x]' in (K.f, K.l, K.q, K.s),
%      where the semidefinite (SDP) cone is scaled by 'mu', see vsdp.svec.  If
%      'do_weight' is true, then for the SDP cone the minimal eigenvalue is
%      weighted by the number of negative eigenvalues.
%
%   [~, n, vidx, sdp_matrix] = obj.rigorous_lower_cone_bound()  With no input
%      arguments, except for the VSDP object, helper structures are computed.
%      Those are:
%
%        - 'n'  The number of required bounds: K.f + K.l + length(K.q)
%                                                        + length(K.s)
%
%        - 'vidx'  Index vector of the cone diagonal elements.  For this vector
%                  there is:
%
%                    length(vidx) = obj.n     % Like 'obj.c', 'obj.x', ...
%                       sum(vidx) = K.f + K.l + length(K.q) + sum(K.s)
%
%        - 'sdp_matrix'  Index matrix to translate the SDP part of 'lb' to the
%                        length 'sum(K.s)' because each SDP lower bound 'lb(j)'
%                        is applied to the whole diagonal 'eye(K.s(j))'.
%
%                        [1    ]
%                        [  1  ]
%                        [  1  ]   [lb_1]
%           sdp_matrix = [  1  ] * [lb_2],  for  obj.K.s = [1, 3, 2]
%                        [    1]   [lb_3]
%                        [    1]
%
%   See also vsdp, vsdp.rigorous_lower_bound, vsdp.rigorous_upper_bound.
%

% Copyright 2004-2018 Christian Jansson (jansson@tuhh.de)

n  = obj.K.f + obj.K.l + length (obj.K.q) + length (obj.K.s);

% Compute helper structures.
if (nargin == 1)
  % Index vector for perturbation.
  vidx = vsdp.sindex (obj);            % Get only diagonal entries of SDP cones.
  vidx(1:(obj.K.f + obj.K.l)) = true;  % Copy free and linear part directly.
  if (~isempty (obj.K.q))
    % In case of second-order cones, only the first element is perturbed.
    vidx(obj.K.idx.q(:,1)) = true;
  end
  
  % Index matrix to translate the SDP part of 'lb'.
  if (isempty (obj.K.s))
    sdp_matrix = [];
  else
    N = sum (obj.K.s);
    cols = zeros (N, 1);
    cols(cumsum ([1; obj.K.s(1:end-1)])) = 1;
    sdp_matrix = sparse (1:N, cumsum (cols), 1);
  end
  lb = [];  % Just to set a value.
else  % Compute cone lower bound.
  narginchk (4, 4);
  x  = vsdp_indexable (x, obj);
  lb = vsdp_indexable (zeros (n, 1), obj);
  
  % Rigorous lower bounds 'lb' on '[x]' for each cone.
  %
  % Free variables: just 0, already set above.
  % LP cone variables.
  lb.l = inf_ (x.l);
  
  % Second-order cone (SOCP) variables.
  offset = obj.K.f + obj.K.l;
  for j = 1:length (obj.K.q)
    xq = x.q(j);
    lb(j + offset) = inf_ (xq(1) - norm (xq(2:end)));
  end
  
  % SDP cone variables.
  offset = offset + length(obj.K.q);
  for j = 1:length(obj.K.s)
    E_ = inf_ (vsdp.verify_eigsym (vsdp.smat ([], x.s(j), mu)));
    if (weight_sdp)
      lb(j + offset) = min (E_) * sum (E_ < 0);
    else
      lb(j + offset) = min (E_);
    end
  end
  lb = lb.value;
end
