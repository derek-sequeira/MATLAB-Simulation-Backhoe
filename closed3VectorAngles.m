function [thetaBC, thetaCA] = closed3VectorAngles(sideAB, sideBC, sideCA,...
    theta_AB, theta0_BC, theta0_CA)
%[thetaB, thetaC] = closed3VectorAngles(sideAB, sideBC, sideCA, theta_AB,
%   theta0_BC, theta0_CA)
% Given a 3 vector closed loop with a path A->B->C->A, where the magnitudes
% of 3  vectors are known and the angle of the vector AB is known, this 
% function will calculate  the other two angles BC & CA. All input and 
% output angles are taken with respect to the horizontal & measured in 
% radians.
% Input sideAB = magnitude of the vector AB whose angle is known
% Input sideBC = magnitude of the vector BC whose angle is to be solved
% Input sideCA = magnitude of the vector CA whose angle is to be solved
% Input theta_AB = angle for known vector AB (radians)
% Input theta0_BC = initial estimate of angle for vector BC (radians)
% Input theta0_CA = initial estimate of angle for vector CA (radians)
% Output theta_BC = calculated angle for vector BC (radians)
% Output theta_CA = calculated angle for vector CA (radians)

% Version 1: created 25/03/2020. Author: Derek Sequeira
% This MATLAB function calculates 2 unknown angles of a 3 vector system in
% a 2d plane by using closure equations. It solves for two equations.
% One stating that the sum of x-changes equals 0. Another stating that the
% sum of y-changes equals 0.
% If the iteration limit is reached the evaluated angles should not be
% considered.

%--------------------------------------------------------------------------

% Internal parameter: iterationMax = maximum number of Newton-Rhapson steps
% taken before the function gives up.
% Internal parameter: tolerence = minimum acceptable value for the 
% magnitude of the function evaluated at the estimated point

iterationMax = 20;
tolerence = 1e-9;
%--------------------------------------------------------------------------

% Check if inputs are valid
if(~isscalar(sideAB) || (~isreal(sideBC)) || (~isscalar(sideCA)) ...
   || (~isreal(sideCA)) || ~isscalar(theta_AB) || (~isreal(theta0_BC)) ||...
   ~isscalar(theta0_BC) || (~isreal(theta0_CA)))
     error('All input must be real number scalar values.')
end

if( (sideAB < 0) || (sideBC < 0) || (sideCA < 0))
    error('All magnitudes should be postive numbers')
end

%--------------------------------------------------------------------------

longestSide = max([sideAB sideBC sideCA]); % Normalize sides so that 
                                           % accuracy isnt affected by 
sideAB = sideAB/longestSide;               % sides lengths
sideBC = sideBC/longestSide; 
sideCA = sideCA/longestSide;

soln(1) = theta0_BC; %Set unknown angles to initial guesses
soln(2) = theta0_CA;

% Begin Newton-Raphson
for count = 0:iterationMax-1
    if count == iterationMax % if number of iterations is greater 
                             % iteration limit, report an error.
        error('Iteration limit reached. Iteration did not converge')
    end
    
    % Define output angles in terms of solution matrix
    thetaBC = soln(1);
    thetaCA = soln(2);
    
    % Calculate the values of the closure equations for X & Y respectively
    % at current value of thetaBC and thetaCA
    f(1) = sideAB*cos(theta_AB) + sideBC...
        *cos(thetaBC) + sideCA*cos(thetaCA);
    f(2) = sideAB*sin(theta_AB) + sideBC...
        *sin(thetaBC) + sideCA*sin(thetaCA);
    
    %Calculate average error between two functions
    avgerror = sqrt(f(1).^2+f(2).^2);
    
    if avgerror < tolerence % If the evaluated error is currently less than 
                            % our tolerence then the current values for
                            % thetaBC and thetaCA are accepted as the
                            % estimated solution.
        break
    end
    
    % Calculate the jacobian of the two equations with respect to the
    % unknown thetaBC and thetaCA
    Jac(1,1) = -1*sideBC*sin(thetaBC);
    Jac(1,2) = -1*sideCA*sin(thetaCA);
    Jac(2,1) = sideBC*cos(thetaBC);
    Jac(2,2) = sideCA*cos(thetaCA);
    
    %Compute matrix division of the Jacobian into the function
    del_theta = Jac\f';
    % Reassign solution for next iteration
    soln = soln - del_theta';
end
% Assign solutions to thetaBC & thetaCA 
thetaBC = soln(1);
thetaCA = soln(2);
% Ensure that angles are in the range 0-2pi
thetaBC = mod(thetaBC, 2*pi);
thetaCA = mod(thetaCA, 2*pi);

end



