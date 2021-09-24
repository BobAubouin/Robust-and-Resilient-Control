function Blin = bilin2lin(Bu,Bv,Bxu,Bvu,Xmean,Vmean)

Bbilinx = [];
Bbilinv = [];
for i=1:size(Bu,2)
    Bbilinx = cat(2,Bbilinx,Bxu(:,:,i)*Xmean);
    Bbilinv = cat(2,Bbilinv,Bvu(:,:,i)*Vmean);
end
Blin = [Bu + Bbilinx + Bbilinv, Bv];
end