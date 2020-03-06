function img=SPECT2DBackwardProj( projection, S, theta )

% SPECT2DBackwardProj	Calculates the backward projection of a sinogram.
%
%	img=SPECT2DBackwardProj( projection, S[, theta] )
%
% Where 
%	projection	is a projection series
%	S		is a SPECT2D System descriptor (see SPECT2DSystem)
%	theta		is a vector of projection IDs from the descriptor
%
% Example
%	img = SPECT2DBackwardProj( projection, S )
%		Creates a backprojection from a full series
%
%	img = SPECT2DBackwardProj( projection, S, 1 )
%		Creates a backprojection from the first projection image only
%

if (nargin < 3) 
	theta=1:size(S.theta,2);
end;

img=zeros( S.imagesize, S.imagesize );

imgsize= size(img);
offset = ceil((sqrt(2)-1)*max(imgsize)/2); 
newsize= imgsize+2*offset;


if ( max(newsize)*S.step > 2*S.ror )
        error('Invalid SPECT System: the whole image must be in the FOV!');
end


for t=1:size(theta,2)
	bigimg = zeros( newsize );

	proj = projection(:,t)';
	if (S.sigma_intr)
		pc = size(proj,2);
		gauss = 1/(sqrt(2*pi)*S.sigma_intr) * exp( -(((1:pc)-floor(pc/2)) * S.step ).^2 / 2 /S.sigma_intr^2 );
		gauss = gauss(find(gauss>1e-6));
		gc = size(gauss,2);
		sp = conv( proj, gauss);
% ezt le kell itt vágni?
%		proj = sp( (floor(gc/2)+1):(pc+floor(gc/2)) );
		proj = sp;
	end

	pc = floor(size(proj,2)/2);
	ac = floor(newsize(1) / 2);
	if (size(proj,2) > newsize(1))
		bigimg( 1, :)=proj( (pc-ac+1):(pc+ac) );
	else
		bigimg( 1, (ac-pc+1):(ac-pc+size(proj,2)))=proj( : );
	end
	for y=(2:newsize(2))
		bigimg(y,:)=bigimg(1,:);
	end

	if (S.psf_a || S.psf_b)
		pc = size(bigimg,2);
		for d=1:size(bigimg,1)
%			dist = S.ror + ((size(rot,1)-d+0.5) - size(rot,1)/2)*S.step;
			dist = S.ror + ((d+0.5) - size(bigimg,1)/2)*S.step;
			sigma = S.psf_a + S.psf_b * dist;
			gauss = 1/(sqrt(2*pi)*sigma) * exp( -(((1:pc)-(pc/2)) * S.step ).^2 / 2 /sigma^2 );
			gauss = gauss(find(gauss>1e-6));
			gc = size(gauss,2);
			sp = conv( bigimg(d,:), gauss);
			bigimg(d,:) = sp( (floor(gc/2)+1):(pc+floor(gc/2)) );
		end
	end

	rot = myimrotate( bigimg, S.theta(theta(t)), 'bilinear', 'crop' );
%	rot = rot * sum(bigimg(:))/sum(rot(:));

%	rot = myimrotate( bigimg, S.theta(theta(t)), 'bicubic', 'loose' );
%	rot = rot * (sum(bigimg(:))/sum(rot(:)));
%    ns = size(rot);
%    os = size(bigimg);
%    shift = ceil((ns-os)/2);
%    rot = rot( (1+shift):(os(1)+shift), (1+shift):(os(1)+shift) );
%
	img = img + rot( (1+offset):(imgsize(1)+offset), (1+offset):(imgsize(2)+offset) );
end

end
