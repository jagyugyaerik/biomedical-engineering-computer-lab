function dist=ImageDistance_L2(V1, V2)

%V1 should be the reference volume...

if size(V1)==size(V2)
    v1=double(V1(:));
    v2=double(V2(:));
    dist=sqrt(sum((v1-v2).^2)/sum(v1.^2))*100;
else 
    error ('Image dimensions are different');
    dist=NaN;
end
