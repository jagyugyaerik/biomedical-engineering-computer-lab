function dist=ImageDistance_CC(v1, v2)

% ImageDistance_CC	Calculates the CC distance of two images.
%
%      dist=ImageDistance_CC(v1, v2)
%
% Where
%       v1     is an image
%       v2     is another image (same size as v1)
%
% Returns the distance of the images is the correlation coefficient of the 
% two images.
%

if size(v1)==size(v2)
    v1=double(v1(:));
    v2=double(v2(:));
    R=corrcoef(v1,v2);
    if size(R) == [2 2]
	R=R(2);
    end
    dist=100*(1-abs(R));
else 
    error ('Image dimensions are different');
    dist=NaN;
end
