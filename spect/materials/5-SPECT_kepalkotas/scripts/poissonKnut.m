function r = poissonKnut( lambda )

% poissonKnut	Creates a Poisson random number
%
%	r = poissonKnut( lambda )
%
% Where
%	lambda 	is the expected value 
%	r 	is the generated Poisson random value
%
% Note that this function uses Knuth's Poisson random generator, which is
%	slow. Very slow. 
%	(see http://en.wikipedia.org/wiki/Poisson_distribution)
%
	L = exp(-lambda);
	k = 0;
	p = 1;
	running=true;
	while (running)
		k=k+1;
		p=p*rand();
		running = (p > L);
	end;
	r = k-1;
end
