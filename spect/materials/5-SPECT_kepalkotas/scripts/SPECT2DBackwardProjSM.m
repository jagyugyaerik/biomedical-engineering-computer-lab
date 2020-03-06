function r=SPECT2DBackwardProjSM( sino, A, thetaIDs )

% SPECT2DBackwardProjSM	Calculates the backward projection of a sinogram.
%
%	r=SPECT2DBackwardProjSM( sino, A[, thetaIDs] )
%
% Where 
%	sino	    is a sinogram
%	A		    is a SPECT2D System matrix descriptor (see SPECT2DSystemMatrix)
%	thetaIDs	is a vector of projection IDs from the descriptor
%
% Example
%	img = SPECT2DBackwardProjSM( sino, S )
%		Creates a backprojection from a sinogram
%
%	img = SPECT2DBackwardProj( projection, S, 1 )
%		Creates a backprojection from the first projection image only
%
if (nargin < 3) 
  thetaIDs=1:length(A.theta);
end;
if (size(sino) ~= [A.projsize, length(thetaIDs)])
  error('Sinogram dimensions are different from that defined in the system matrix!');
end
r=zeros( A.imagesize, A.imagesize );

if (length(thetaIDs)==length(A.theta))
    r(:) = sino(:)'*A.matrix';
elseif (length(thetaIDs)==1)
    r(:) = sino(:)'*A.matrix(:,((thetaIDs-1)*A.projsize+1):(thetaIDs*A.projsize))';
else
    for i=1:length(thetaIDs)
        r = r+SPECT2DBackwardProjSM( sino(:,i), A, thetaIDs(i) );
    end
end

end
