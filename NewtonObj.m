% Object representing an object in 3D space using Newtonian mechanics

classdef NewtonObj < handle
	methods
		% Constructor.
		%
		% Parameters:
		%     name        Name for the object -- affects symbolic variable naming
		%     daetool     LinDAETool instance to embed this object's dynamics into
		%     init_pos    Initial position for this object
		%     init_orient Initial orientation for this object (symbolic quaternion)
		function self = NewtonObj(name, daetool, init_pos, init_orient)
			
		end
	end

	properties
		% Symbolic translation variables
		pos
		vel
		accel

		% Symbolic orientation variables
		orient
		ang_vel
		ang_accel

		% External force and torque expressions
		ext_force
		ext_torque
	end
end
