% A tool for formulating and solving Differential-Algebraic Equations
% where the system of equations is linear.

classdef LinDAETool < handle
	methods
		% Constructor
		function self = LinDAETool
			self.y = sym([]);
			self.dy = sym([]);
			self.eqns = sym([]);
		end

		% Adds a new state to the DAE. Needs a symbolic name for the state
		% as well as the value at time 0. Returns symbolic variable for the
		% state and its derivative.
		function [expr, dexpr] = addVar(self, name, init_value)
			expr = sym(name, 'real');
			dexpr = sym(['d' name], 'real');
			self.y(end+1, 1) = expr;
			self.dy(end+1, 1) = dexpr;
			self.y0(end+1, 1) = init_value;
		end

		% Adds a new equation to the DAE system.
		function addEqn(self, lhs, rhs)
			self.eqns(end+1, 1) = lhs - rhs;
		end

		% Solves the DAE for the given duration and caches the solution.
		function solve(self, duration)
			M = matlabFunction(jacobian(self.eqns, self.dy), 'vars', {sym([]), self.y});
			f = matlabFunction(-subs(self.eqns, self.dy, zeros(size(self.dy))), 'vars', {sym({}), self.y});
			self.soln = ode15s(f, [0 duration], self.y0, odeset('Mass', M, 'MassSingular', 'no'));
		end
	end

	properties
		y % State variables representation (symbolic)
		dy % Derivative of the state variables (symbolic)
		y0 % Value of y at time 0
		eqns % DAE equations. Will be set to 0 at runtime. Must be linear in dy
		soln % Solution to the DAE, as returned by MATLAB's ODE solver
	end
end
