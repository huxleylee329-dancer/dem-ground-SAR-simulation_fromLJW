function [image2_quantify,SlaveReg,valid,info] = regis_subpixel(Master,Slave,nr,nc)
assert(isequal(size(Master),size(Slave),[nr,nc]));
assert(all(isfinite(Master(:))) && all(isfinite(Slave(:))));
bs=32; stride=16; lim=2; border=4;
assert(nr>=bs && nc>=bs,'Image is smaller than one block.');
rStarts=unique([1:stride:nr-bs+1,nr-bs+1]);
cStarts=unique([1:stride:nc-bs+1,nc-bs+1]);
[Xq,Yq]=meshgrid(1+border:bs-border);
peakPower=max(abs(Master(:)).^2);
points=zeros(numel(rStarts)*numel(cStarts),5); count=0;
options=optimset('Display','off','TolX',1e-3,'TolFun',1e-7,'MaxIter',150,'MaxFunEvals',300);
for ir=1:numel(rStarts)
 for ic=1:numel(cStarts)
  rr=rStarts(ir)+(0:bs-1); cc=cStarts(ic)+(0:bs-1);
  patchM=Master(rr,cc); patchS=Slave(rr,cc);
  A=abs(patchM(1+border:bs-border,1+border:bs-border));
  if mean(A(:).^2)<1e-3*peakPower,continue;end
  a=A(:)-mean(A(:)); na=norm(a);
  if na<0.05*norm(A(:)) || na==0,continue;end
  cost=@(s) patch_cost(s,patchS,Xq,Yq,a,na,lim);
  best=Inf; start=[0,0];
  for dr=-lim:lim
   for dc=-lim:lim
    f=cost([dr,dc]);
    if f<best,best=f;start=[dr,dc];end
   end
  end
  [s,f,flag]=fminsearch(cost,start,options);
  score=-f;
  if flag<=0 || score<0.6 || any(abs(s)>lim-0.05),continue;end
  count=count+1;
  points(count,:)=[mean(rr),mean(cc),s(1),s(2),score];
 end
end
points=points(1:count,:);
assert(count>=12,'Fewer than 12 reliable blocks: increase useful scene coverage or inspect image matching.');
mr=(points(:,1)-(nr+1)/2)/nr;
nc0=(points(:,2)-(nc+1)/2)/nc;
F=[ones(count,1),mr,nc0,mr.^2,nc0.^2,mr.*nc0];
keep=true(count,1);
for it=1:4
 assert(sum(keep)>=12 && rank(F(keep,:))==6,'Matched block distribution cannot support a quadratic model.');
 assert(numel(unique(points(keep,1)))>=3 && numel(unique(points(keep,2)))>=3,'Need matches in at least three rows and three columns.');
 w=points(keep,5);
 coef=(F(keep,:).*w)\(points(keep,3:4).*w);
 residual=sqrt(sum((F*coef-points(:,3:4)).^2,2));
 med=median(residual(keep));
 robustScale=1.4826*median(abs(residual(keep)-med));
 candidate=keep & residual<=med+max(0.1,3*robustScale);
 if isequal(candidate,keep),break;end
 keep=candidate;
end
assert(sum(keep)>=12 && rank(F(keep,:))==6,'Too few consistent matches after rejection.');
w=points(keep,5);
coef=(F(keep,:).*w)\(points(keep,3:4).*w);
[C,R]=meshgrid(1:nc,1:nr);
m=(R(:)-(nr+1)/2)/nr; n=(C(:)-(nc+1)/2)/nc;
G=[ones(numel(m),1),m,n,m.^2,n.^2,m.*n];
shift=G*coef;
dRow=reshape(shift(:,1),nr,nc); dCol=reshape(shift(:,2),nr,nc);
SlaveReg=interp2(Slave,C+dCol,R+dRow,'cubic',NaN);
hull=convhull(points(keep,2),points(keep,1));
accepted=points(keep,:);
inside=inpolygon(C,R,accepted(hull,2),accepted(hull,1));
valid=isfinite(SlaveReg) & inside & abs(dRow)<=lim & abs(dCol)<=lim;
SlaveReg(~valid)=NaN;
amp=abs(SlaveReg);amp(~valid)=0;
image2_quantify=20*log10(amp/max(max(amp(:)),eps)+eps);
info.points=points;info.keep=keep;info.dRow=dRow;info.dCol=dCol;
info.accepted=sum(keep);info.candidates=count;info.blocks=numel(rStarts)*numel(cStarts);
end
function f=patch_cost(s,S,X,Y,a,na,lim)
if any(abs(s)>lim),f=1e3+sum(s.^2);return;end
z=interp2(S,X+s(2),Y+s(1),'cubic');
b=abs(z(:)); b=b-mean(b); nb=norm(b);
if nb==0 || any(~isfinite(b)),f=1e3;else,f=-real(a'*b)/(na*nb);end
end
