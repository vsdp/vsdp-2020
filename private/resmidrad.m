function [R,Rrad] = resmidrad(Amid,Bmid,Cmid,Dmid,Arad,Brad,Crad,Drad)
% helper function to compute a rigorous midpoint radius enclosure for
%       R = A*B-C*D
%
%   [R,Rrad] = resmidrad(Amid,Bmid,Cmid,Dmid)
%   [R,Rrad] = resmidrad(Amid,Bmid,Cmid,Dmid,Arad,Brad,Crad,Drad)
%
% this function is used within VSDP (for efficiency reasons)
% Important!    setround(1) is assumed


%% calculate inclusion for point-matrices
if isreal(Amid) && isreal(Bmid) && isreal(Cmid) && isreal(Dmid)
  % sup(A*B-C*D)
  R = Amid*Bmid + Cmid*(-Dmid);
  % -inf(A*B-C*D)  - Rrad is used as placeholder
  Rrad = Amid*(-Bmid) + Cmid*Dmid;
  % 0.5*(sup(A*B-C*D)+inf(A*B-C*D)) = mid(A*B-C*D)
  R = 0.5*(R-Rrad);
  % Rrad = mid(A*B-C*D) - inf(A*B-C*D)
  Rrad = R + Rrad;
else
  % extract imaginary part of A
  if isreal(Amid)
    iAmid = sparse(size(Amid,1),size(Amid,2));
  else
    iAmid = imag(Amid);
  end
  % extract imaginary part of C
  if isreal(Cmid)
    iCmid = sparse(size(Cmid,1),size(Cmid,2));
  else
    iCmid = imag(Cmid);
  end
  
  % sup(A*B-C*D)
  R = real(Amid)*Bmid + iAmid*(1j*Bmid) + ...
    (-real(Cmid))*Dmid + iCmid*(-1j*Dmid);
  % -inf(A*B-C*D)  - Rrad is used as placeholder
  Rrad = (-real(Amid))*Bmid + iAmid*(-1j*Bmid) + ...
    real(Cmid)*Dmid + iCmid*(1j*Dmid);
  % 0.5*(sup(A*B-C*D)+inf(A*B-C*D)) = mid(A*B-C*D)
  R = 0.5*(R-Rrad);
  % Rrad = abs(mid(A*B-C*D)-inf(A*B-C*D))
  Rrad = abs(R+Rrad);
end

% matrices are point-matrices   ->  finish!
if nargin<8
  return;
end


%% extend interval radius
if ~isempty(find(Arad,1))  % regard Arad
  Rrad = Rrad + Arad*(abs(Bmid)+Brad);
end

if ~isempty(find(Brad,1))  % regard Brad
  Rrad = Rrad + abs(Amid)*Brad;
end

if ~isempty(find(Crad,1))  % regard Crad
  Rrad = Rrad + Crad*(abs(Dmid)+Drad);
end

if ~isempty(find(Drad,1))  % regard Drad
  Rrad = Rrad + abs(Cmid)*Drad;
end
