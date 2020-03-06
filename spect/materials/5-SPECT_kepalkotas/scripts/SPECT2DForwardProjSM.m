function p=SPECT2DForwardProjSM( img, A, thetaIDs)

if (size(img) ~= [A.imagesize, A.imagesize])
	error('Image dimensions are different from that defined in the system matrix!');
end
if (nargin < 3) 
	thetaIDs=1:size(A.theta,2);
end;

p=zeros( A.projsize, length(thetaIDs));
if (length(thetaIDs)==length(A.theta))
	p(:) = img(:)'*A.matrix;
elseif (length(thetaIDs)==1)
	p(:) = img(:)'*A.matrix(:, ((thetaIDs-1)*A.projsize+1):(thetaIDs*A.projsize));
else
	for t=1:length(thetaIDs)
        p(:,t) = SPECT2DForwardProjSM( img, A, thetaIDs(t) );
	end
end

end
