function r = poissonRejection( lambda )

% poissonRejection    Creates a Poisson random number
%
%	r = poissonRejection( lambda )
%
% Where
%	lambda 	is the expected value 
%	r 	is the generated Poisson random value
%
%
        c = 0.767 - 3.36/lambda;
        beta = pi/sqrt(3.0*lambda);
        alpha = beta*lambda;
        k = log(c) - lambda - log(beta);
        lhs = 1; rhs = 0;

        while (lhs > rhs)
            u = rand();
            x = (alpha - log((1.0 - u)/u))/beta;
            n = floor(x + 0.5);
            if (n < 0)
                continue;
            end
            v = rand();
            y = alpha - beta*x;
            lhs = y + log(v/(1.0 + exp(y))^2);
            rhs = k + n * log( lambda ) - gammaln( n + 1.0 );
        end
        r=n;
end
