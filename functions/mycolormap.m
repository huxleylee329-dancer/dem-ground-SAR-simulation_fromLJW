function pic=mycolormap(amp,phase,nr,naz)
  % Function to overlay phase image on the amplitude image
  
amplow=prctile(amp(:),5); %amp(:)中有5%的数据小于amplow
amphi=prctile(amp(:),95); %amp(:)中有95%的数据小于amphi
amp(amp<amplow)=amplow;
amp(amp>amphi)=amphi;
scale=(amp-amplow)/(amphi-amplow);
phlow=prctile(phase(:),1);
phhi=prctile(phase(:),99);

%creat a color table
colormap jet;
map=colormap;
stack=max(phase,phlow); %phase中小于phlow的替换成phlow
stack=min(phase,phhi);
colorstack=round((stack-phlow)/(phhi-phlow)*64);
colorstack=max(colorstack,1);
colorstack=min(colorstack,64);

for k=1:naz
    for kk=1:nr
        red(k,kk)=map(colorstack(k,kk),1);
        green(k,kk)=map(colorstack(k,kk),2);
        blue(k,kk)=map(colorstack(k,kk),3);
    end
end
pic(:,:,1)=red.*scale;
pic(:,:,2)=green.*scale;
pic(:,:,3)=blue.*scale;

end