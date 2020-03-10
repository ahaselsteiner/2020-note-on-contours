function [xcont,ycont,plevel,Pmax]=hdr2D(x,y,fxy,alpha)

% calculate included probability as a function of density level
dx=x(2)-x(1);
dy=y(2)-y(1);
psorted=flip(sort(fxy(:)));
P=cumsum(psorted)*dx*dy;
psorted=flip(psorted);
P=flip(P);

% check on accuracy of simple numerical integration
Pmax=max(P);

% find density level corresponding to correct exclusion probability
for i=1:length(alpha)
    ind=find(P<1-alpha(i),1,'first');
    plevel(i)=psorted(ind);
    
    C=contourc(x,y,fxy,[plevel plevel]);
    n=C(2,1);
    xcont{i,1}=C(1,2:n+1);
    ycont{i,1}=C(2,2:n+1);
end

end