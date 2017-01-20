function R = prodsup(Amid,Bmid,Arad,Brad)
% helper function to compute a rigorous upper bound (supremum) for
%       R = A*B
%
%   R = prodsup(Amid,Bmid)
%   R = prodsup(Amid,Bmid,Arad,Brad)
%
% this function is used within VSDP (for efficiency reasons)
% Important!    setround(1) is assumed

% sup(Amid*Bmid)
if ~isreal(Amid) && ~isreal(Bmid)
  R = real(Amid)*Bmid + imag(Amid)*(1j*Bmid);
else
  R = Amid*Bmid;
end

% matrices are point-matrices   ->  finish!
if nargin<4
  return;
end

% initial interval radius
Rrad = sparse(size(R,1),size(R,2));

if ~isempty(find(Arad,1))  % regard Arad
  Rrad = Rrad + Arad*(abs(Bmid)+Brad);
end

if ~isempty(find(Brad,1))  % regard Brad
  Rrad = Rrad + abs(Amid)*Brad;
end

R = R + Rrad;
