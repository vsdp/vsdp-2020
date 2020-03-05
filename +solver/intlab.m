classdef intlab < handle
  % INTLAB  Solver proxy class (not the actual solver!).
  %
  % For information about INTLAB, see:
  %
  %      http://www.ti3.tuhh.de/rump/intlab
  %
  %   See also vsdp.solve.
  %
  
  
  % Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)
  methods (Static)
    function spath = install (varargin)
      % Returns the path to the installed and usable solver.  Otherwise return
      % an empty array.  No error messages are thrown.
      %
      % By passing one or more arguments interactive installation actions
      % happen and, in case of failures, error messages are thrown.
      %
      
      sname          = 'intlab';
      is_available   = @() exist ('startintlab', 'file') == 2;
      get_path       = @() fileparts (which ('startintlab'));
      installer_file = 'startintlab.m';
      do_error       = true;  % Always error if INTLAB is missing!
      
      try
        spath = solver.registry.generic_install (sname, is_available, ...
          get_path, installer_file, do_error);
      catch
        error ('VSDP:SOLVER:INTLAB:missing', ...
          '%s.  %s\n\n\t%s\n\n%s.\n\n', ...
          'Cannot find interval toolbox "INTLAB"', ...
          'Get a recent version from', ...
          'http://www.ti3.tuhh.de/rump/intlab', ...
          'and run "startintlab" from the root directory');
      end
    end
  end
end
