% Convenience class for assembling SFDyn objects.
% Allows for easy coordinate creation, energy specification,
% and creates the SFDyn object or acceleration function for the user.

classdef DynTool < handle
	methods (Access = public)
		% The constructor for this class; does some basic initialization
		%
		% Parameters: none
		%
		% Returns:
		%     this A handle to this object
		function this = DynTool
			this.lagrangian = sym(0);
			this.constraints = sym([]);
		end % LagrDyn

		% This function adds a new coordinate to this dynamics system.
		%
		% Parameters:
		%     name The symbolic variable name for this parameter
		%
		% Returns:
		%     expr  The symbolic expression for this coordinate
		%     dexpr The symbolic expression for this coordinate's derivative
		function [expr,dexpr] = addCoord(this, name)
			% Error checking
			if nargin < 2
				error('addCoord needs a name parameter.')
			end
			if ~ischar(name)
				error('Input "name" must be a string')
			end

			% Diagnostic message for the user
			disp(['Adding coordinate ' name])

			% Generate symbolic expressions
			expr   = sym(name,        'real');
			dexpr  = sym(['d'  name], 'real');

			% Append this coordinate to the coordinates list
			this.coords(end+1,1)   = expr;
			this.dcoords(end+1,1)  = dexpr;
		end % addCoord

		% Add a new kinetic energy term
		%
		% Parameters:
		%     expr Kinetic energy expression
		%
		% Returns:
		%     expr The kinetic energy expression
		function expr = addKE(this, expr)
			% Diagnostic message for the user
			disp(['Adding kinetic energy expression ' char(expr)])

			% Directly add the expression to the Lagrangian
			this.lagrangian = this.lagrangian + expr;
		end % addKE

		% Add a new potential energy term
		%
		% Parameters:
		%     expr Potential energy expression
		%
		% Returns:
		%     expr The potential energy expression
		function expr = addPE(this, expr)
			% Diagnostic message for the user
			disp(['Adding potential energy expression ' char(expr)])

			% Directly add (well, subtract) the expression from the Lagrangian
			this.lagrangian = this.lagrangian - expr;
		end % addPE

		% Adds a new constraint
		%
		% Parameters:
		%     expr  The expression for the constraint. Dynamically set to 0
		%
		% Returns:
		%     expr  The constraint expression.
		function expr = addConstraint(self, expr)
			% Diagnostic message for the user
			disp(['Adding constraint expression ' char(expr)])

			% Directly add the constraint to the symbolic constraint vector
			self.constraints(end+1) = expr;
		end

		% Generates and returns the acceleration function for a constrained system
		%
		% The returned function takes in the positions and velocities
		% and computes the accelerations
		function accel_fcn = genAccelFcn(self)
			% The linear systems of equations takes the form
			%     lhsMat * [ddx; lambda] = rhs

			% Various Jacobians which are incorporated into lhsMat and rhs
			dL_ddr = jacobian(self.lagrangian, self.dcoords);
			df_dq = jacobian(self.constraints, self.coords);
			dt_df_dq_qcomp = reshape(...
				jacobian(df_dq, self.coords) * self.dcoords...
			, numel(self.constraints), numel(self.coords));

			% Compute the left-hand side matrix, which includes both
			% the dynamics and the constraint equations.
			lhsMat = [ jacobian(dL_ddr, self.dcoords), df_dq.'
			           df_dq,           zeros(numel(self.constraints) * [1 1]) ];

			% Right hand side of the system of equations
			rhs = [ jacobian(self.lagrangian, self.coords).' - jacobian(dL_ddr, self.coords) * self.dcoords
			        -dt_df_dq_qcomp * self.dcoords ];

			% Generate anonymous functions to calculate lhsMat and rhs
			lhsMatFcn = matlabFunction(lhsMat, 'vars', {self.coords, self.dcoords});
			rhsFcn = matlabFunction(rhs, 'vars', {self.coords, self.dcoords});

			% Generate the anonymous acceleration function
			subsref_struct.type = '()';
			subsref_struct.subs = {1:numel(self.coords)};
			accel_fcn = @(pos, vel) subsref(lhsMatFcn(pos, vel) \ rhsFcn(pos, vel), subsref_struct);
		end
		% and returns the acceleration vector.

		% Generates and returns the SFDyn object.
		% Ignores constraints!!!
		function sfdyn = genSFDyn(self)
			% Pass our stored values to SFDyn's constructor, which does all the work
			sfdyn = SFDyn(self.lagrangian, self.coords, self.dcoords);
		end
	end % methods

	properties (GetAccess = public, SetAccess = private)
		% Coordinates and their derivatives
		coords@sym
		dcoords@sym

		% Lagrangian expression
		lagrangian@sym

		% Dynamic constraints. These are set equal to 0 when solving for the acceleration function
		constraints@sym
	end % properties
end % classdef
