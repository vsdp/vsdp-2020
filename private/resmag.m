function R = resmag(Amid,Bmid,Cmid,Dmid,Arad,Brad,Crad,Drad)
% helper function to compute a rigorous upper bound for
%       R = magnitude(A*B-C*D),
%
%   R = resmag(Amid,Bmid,Cmid,Dmid)
%   R = resmag(Amid,Bmid,Cmid,Dmid,Arad,Brad,Crad,Drad)
%
% this function is used within VSDP (for efficiency reasons)
% Important!    setround(1) is assumed

if (isreal(Amid) && isreal(Bmid) && isreal(Cmid) && isreal(Dmid))
    % max(sup(A*B-C*D),-inf(A*B-C*D)) = abs(A*B-C*D)
    R = Amid*Bmid + Cmid*(-Dmid);
    R = max(Amid*(-Bmid) + Cmid*Dmid, R);
else
    % extract imaginary part of A
    if (isreal(Amid))
        iAmid = sparse(size(Amid,1),size(Amid,2));
    else
        iAmid = imag(Amid);
    end
    % extract imaginary part of C
    if (isreal(Cmid))
        iCmid = sparse(size(Cmid,1),size(Cmid,2));
    else
        iCmid = imag(Cmid);
    end

    % sup(A*B-C*D)  - Rrad is used as placeholder
    Rrad = real(Amid)*Bmid + iAmid*(1j*Bmid) + ...
        (-real(Cmid))*Dmid + iCmid*(-1j*Dmid);
    % -inf(A*B-C*D)
    R = (-real(Amid))*Bmid + iAmid*(-1j*Bmid) + ...
        real(Cmid)*Dmid + iCmid*(1j*Dmid);
    % mag(A*B-C*D)
    R = abs( max(real(Rrad),real(R)) + 1j*max(imag(Rrad),imag(R)) );
end

% matrices are point-matrices   ->  finish!
if (nargin<8)
    return;
end

% initial interval radius
Rrad = sparse(size(R,1),size(R,2));

if (any(any(Arad)))  % regard Arad
    Rrad = Rrad + Arad*(abs(Bmid)+Brad);
end

if (any(any(Brad)))  % regard Brad
    Rrad = Rrad + abs(Amid)*Brad;
end

if (any(any(Crad)))  % regard Crad
    Rrad = Rrad + Crad*(abs(Dmid)+Drad);
end

if (any(any(Drad)))  % regard Drad
    Rrad = Rrad + abs(Cmid)*Drad;
end

R = R + Rrad;