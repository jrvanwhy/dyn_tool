% Class which represents dynamics in Standard Form, as described in:
% http://www.cds.caltech.edu/~murray/books/MLS/pdf/mls94-manipdyn_v1_2.pdf
%
% In short, the dynamics take the following form:
%
%     M(q) * ddq + C(q, dq) * dq + N(q, dq) = u
%
% where q is the position vector, dq is the velocity vector, ddq is the acceleration vector,
% and u is the input force vector.

classdef SFDyn < handle
	methods
		% Constructor; this builds up the system dynamics
		%
		% Parameters:
		%     Lagr  Lagrangian of the system
		%     pos   Position vector (vector of symbolic variables)
		%     vel   Velocity vector (vector of symbolic variables)
		function self = SFDyn(Lagr, pos, vel)
			% Jacobian of Lagr with respect to velocity, used for the mass matrix
			% and Coriolis matrix calculations
			D_vel_Lagr = jacobian(Lagr, vel);

			% Generate the mass matrix
			self.mass_matrix = jacobian(D_vel_Lagr, vel);

			% Generate the Coriolis matrix
			self.coriolis_matrix = jacobian(D_vel_Lagr, pos);

			% Generate the position-dependent (gravity) terms function
			self.pos_dep_terms = jacobian(-Lagr, pos).';

			% Store the position and velocity variables
			self.pos_vars = pos;
			self.vel_vars = vel;
		end

		% Generates an anonymous function that evaluates this system's dynamics.
		%
		% The function takes in the positions, velocities, and forces
		% and returns the acceleration vector.
		function accel_fcn = gen_accel_fcn(self)
			% Mass matrix function
			m_fcn = matlabFunction(self.mass_matrix, 'vars', {self.pos_vars, self.vel_vars});

			% Coriolis and position-dependent terms sum function
			rhs_fcn = matlabFunction(self.coriolis_matrix * self.vel_vars + self.pos_dep_terms, 'vars', {self.pos_vars, self.vel_vars});

			% Acceleration function. Solves numerically at runtime (this should generally
			% be more accurate and faster than a symbolic solution)
			accel_fcn = @(pos, vel, forces) m_fcn(pos, vel) \ (forces - rhs_fcn(pos, vel));
		end
	end

	properties
		mass_matrix
		coriolis_matrix
		pos_dep_terms
		pos_vars
		vel_vars
	end
end
