function r = poissonRand( lambda )

% poissonKnut	Creates a Poisson random number
%
%	r = poissonRand( lambda )
%
% Where
%	lambda 	is the expected value 
%	r 	is the generated Poisson random value
%
% This function is using either Knut's algorithm (lambda <= 30) or a
% regection method (otherwise).
%
	if (lambda > 30)
		r = poissonRejection( lambda );
	else
		r = poissonKnut( lambda );
	end
end
