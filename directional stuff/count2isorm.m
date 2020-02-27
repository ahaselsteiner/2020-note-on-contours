function [xcont,ycont]=count2iform(x,y,count,P)

% conditional CDF of y given x
CDFygx=cumsum(count,1);
m=repmat(CDFygx(end,:),length(y),1);
CDFygx=CDFygx./m;
bad=sum(isnan(CDFygx),1)>0;
i0=find(bad==0,1,'first');
i1=find(bad(i0:end)==1,1,'first')+i0-2;
if isempty(i1)
    i1=length(x);
end
CDFygx=CDFygx(:,i0:i1);
xi=x(i0:i1);

% CDF of x
CDFx=sum(cumsum(count,2),1);
CDFx=CDFx/max(CDFx);

% contour points in normal space
r=sqrt(chi2inv(P,2));
theta=(0:1:360)*pi/180;
ux=r*cos(theta);
uy=r*sin(theta);
Px=normcdf(ux);
Py=normcdf(uy);

% contour points on original margins
[CDFx,inds]=unique(CDFx);
xcont=interp1(CDFx,x(inds),Px,'pchip');

ycont=0*xcont;
for i=1:length(uy)
%     disp(i)
    if xcont(i)<min(xi)
        condCDFy=CDFygx(:,1);
    else
        condCDFy=interp2(xi,y,CDFygx,xcont(i),y);
    end
%     [condCDFy,inds]=unique(condCDFy);
%     yi=y(inds);
%     if any(isnan(yi)) || any(isnan(condCDFy)) || any(isinf(yi)) || any(isinf(condCDFy))
%         error('got')
%     end
%     ycont(i)=interp1(condCDFy,yi,Py(i),'pchip');
    i0=find(condCDFy<Py(i),1,'last');
    i1=find(condCDFy>Py(i),1,'first');
    if isempty(i0) 
        i0=1;
    end
    if isempty(i1)
        disp('err')
    end
    if i1==i0
        ycont(i)=y(i0);
    else
        ycont(i)=interp1(condCDFy(i0:i1),y(i0:i1),Py(i));
    end
end