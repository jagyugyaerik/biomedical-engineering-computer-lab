function noisy = SPECT2DAddNoise( proj, count )

% SPECT2DAddNoise Adds Poisson noise to a given projection series
%
%	noisy = SPECT2DAddNoise( proj, count )
%
% Where
%	proj is the original (noiseless) projection series
%	count is the mean total count of the generated projection series
%	noisy is the generated, noisy projection series
%
% Note: that this can take some time, when the statistics toolbox is not
%	installed.
%

if (nargin < 2)
	count=100000
end

noisy = proj / sum(proj(:)) * count;

if (exist('random')>=2)
	% if your Matlab has the random function from the statistics 
	% toolbox, then we can use that
	noisy = random('Poisson', noisy);
else
	% otherwise we use our own (very slow) implementation
	for i=1:prod(size(noisy))
		noisy(i) = poissonRand(noisy(i));
	end
end

end
