function [SRC, DSRC] = stechmann_functions_general(T, F, D)
% FN_functions_general.m
% A simple function to define the source and Jacobian functions for the
% general form of the Stechmann equation in a form compatible
% with SPDE_IE_2D.  In general, the source term is
%
%    f = -1/T * q + F + D W',
%
% Where T, F, and D may be functions of space and time.  The functions
% defined here assume that the functions T, F, and D all take arguments
% (x,y,t) for vectors x and y and scalar t and return matrices with
% dimensions 1 by length(x) by length(y) (note the 1; it's important).

% Author: Cooper Brown
    
    % Function handles to be returned
    SRC = @(x,y,t,u,w) SRC_defn(x,y,t,u,w);
    DSRC = @(x,y,t,u,~) DSRC_defn(x,y,t,u);

    function f = SRC_defn(x, y, t, u, w)
        % Source term for the Stechmann equation
        f = -1./T(x,y,t) .* u + F(x,y,t) + D(x,y,t).*w;
    end

    function f = DSRC_defn(x, y, t, u)
        % Jacobian matrix for the Stechman equation
        f = -1./T(x,y,t);
    end
end