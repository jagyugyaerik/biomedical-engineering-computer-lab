function subsetOrder = SPECTSubsetOrder(nSubsets)
	done = zeros([nSubsets 1]);
	subsetOrder = 1:(nSubsets);
        subsetOrder(1) = 1; % start with arbitrary subset
        done(1) = 1;
        for i=2:nSubsets
            % check all remaining subsets and take one with maximum distance
            % to previous subset
            ibest = 0;
            dbest = 0;
            j = 1;
            while (j < nSubsets)
                % next unprocessed subset
                while ((done(j) == 1) && (j < nSubsets))
                    j=j+1;
		end % while
		
                if (done(j) == 0) % we have found a yet unprocessed subset
                    % find distance to previous subset
                    prevSet = subsetOrder(i-1);
                    d = abs(j - prevSet);
                    if ((nSubsets - d) < d)
                        d = nSubsets - d;
                    end % if

                    if (d > dbest) % check if this is the largest distance found so far
                        dbest = d;
                        ibest = j;
                    end % if
	            j = j+1;
                end % if
            end % while

            subsetOrder(i) = ibest;
            done(ibest) = 1;

        end % for subsets
end % function
