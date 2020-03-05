classdef mosek < handle
  % MOSEK  Solver proxy class (not the actual solver!).
  %
  %   For more information on the MOSEK format, see:
  %
  %      [1] https://docs.mosek.com/8.1/toolbox/data-types.html.
  %
  %   See also vsdp.solve.
  %
  
  
  % Copyright 2004-2020 Christian Jansson (jansson@tuhh.de)
  methods (Static)
    function obj = solve (obj, sol_type)
      % SOLVE  Approximately solve conic problem instance with MOSEK.
      %
      %   See also vsdp.solve.
      %
      
      narginchk (1, 2);
      
      solver.mosek.install (true);                   % Show errors
      solver.registry.check_cones (obj, 'mosek', 1); % Show errors
      
      if (nargin == 1)
        sol_type = 'Approximate';
      end
      [A, b, c] = obj.get_midpoint_problem_data (sol_type);
      
      % Should initial solution guess be taken into account?
      if (obj.options.USE_INITIAL_GUESS)
        warning ('VSDP:solve_mosek:ignoreInitialGuess', ...
          ['solve_mosek: MOSEK does not support initial guesses (x0,y0,z0) ' ...
          'and proceeds without them.']);
      end
      
      % Should special solver options be taken into account?
      if (~isempty (obj.options.SOLVER_OPTIONS))
        pars = obj.options.SOLVER_OPTIONS;
      else
        pars = [];
      end
      
      % Adapt output verbosity.
      if (~obj.options.VERBOSE_OUTPUT)
        pars.fid = 0;
      end
      
      % Prepare data for solver.
      
      % Regard the free, linear, and second-order cone.
      idx = obj.K.f + obj.K.l + sum(obj.K.q);
      prob.a = sparse(A(1:idx,:)');
      prob.c = c(1:idx,1);
      
      % Define lower bounds for free, linear, and second-order cone variables.
      prob.blx = [-inf * ones(1, obj.K.f), zeros(1, obj.K.l), ...
        -inf * ones(1, sum (obj.K.q))];
      
      % If there are second-order cones.
      if (sum (obj.K.q) > 0)
        [~, res] = mosekopt('symbcon');
        prob.cones.type = repmat (res.symbcon.MSK_CT_QUAD, 1, length (obj.K.q));
        % The indices of all second-order cone variables.
        prob.cones.sub  = obj.K.idx.q(1,1):obj.K.idx.q(end,end);
        % The indices of all second-order cone starts.  Those are the indices
        % relative to 'prob.cones.sub', not the "global" indices.
        prob.cones.subptr = obj.K.idx.q(:,1)' - obj.K.idx.q(1,1) + 1;
      end
      
      % If there are semidefinite cones.
      if (sum (obj.K.s) > 0)
        prob.bardim = obj.K.s;
        [~, prob.barc.subj, ...
          prob.barc.subk, prob.barc.subl, prob.barc.val] = ...
          solver.mosek.to_mosek_fmt (obj, vsdp.smat (obj, c, 1));
        
        [prob.bara.subi, prob.bara.subj, ...
          prob.bara.subk, prob.bara.subl, prob.bara.val] = ...
          solver.mosek.to_mosek_fmt (obj, vsdp.smat (obj, A, 1));
      end
      
      % As the VSDP format support equality constraints only, we have to set
      % the lower and upper bound for 'A*x' equal to 'b'.
      prob.blc = b;
      prob.buc = b;
      
      % Adapt output verbosity.
      if (obj.options.VERBOSE_OUTPUT)
        echo_level = 3;
      else
        echo_level = 0;
      end
      
      % Call solver.
      tic;
      [r, res] = mosekopt ( ...
        sprintf ('minimize info echo(%d)', echo_level), prob);
      solver_info.elapsed_time = toc;
      
      % Store solution after normal termination.
      if (isscalar(r) && isnumeric(r))
        solver_info.name = 'mosek';
        if (r == 0)
          solver_info.termination = 'Normal termination';
        else
          solver_info.termination = 'Unknown';
        end
        % In case of a pure LP, prefer the basic solution 'bas' (simplex
        % optimizer).  Otherwise use the interior-point solution 'itr' if
        % available.
        if (isfield (res, 'sol') && isfield (res.sol, 'bas'))
          x = res.sol.bas.xx(:);
          y = res.sol.bas.y(:);
        elseif (isfield (res, 'sol') && isfield (res.sol, 'itr'))
          x = [res.sol.itr.xx(:); solver.mosek.to_vsdp_fmt(obj, ...
            res.sol.itr.barx)];
          y = res.sol.itr.y(:);
        else
          obj.add_solution (sol_type, [], [], [], nan(2,1), solver_info);
          warning ('VSDP:solve_mosek:solutionNotAvailable', ...
            'solve_mosek: No solution returned by MOSEK.');
          return;
        end
        z = vsdp.svec (obj, obj.c - obj.At*y, 1);
        f_objective = [obj.c'*x; obj.b'*y];
        obj.add_solution (sol_type, x, y, z, f_objective, solver_info);
      end
      
    end
    
    
    function [subi, subj, subk, subl, val] = to_mosek_fmt (obj, X)
      % TO_MOSEK_FMT Translate the VSDP quantities 'At' and 'c' to MOSEK format.
      %
      %   The output are four to five vectors of equal length representing data
      %   for the non-zero values of the VSDP quantities
      %
      %      subi - constraint index (only relevant for 'At')
      %      subj - cone       index
      %      subk - row        index
      %      subl - column     index
      %      val  - non-zero  values
      %
      %   See also vsdp, vsdp.sindex, mosek_idx.
      %
      
      
      [ltri, subk_all, subl_all] = cellfun(@solver.mosek.index, ...
        num2cell (obj.K.s), 'UniformOutput', false);
      ltri     = vertcat (ltri{:});   % lower triangular    indices (unshifted)
      subj_all = ones (size (ltri));  % all possible cone   indices (unshifted)
      subk_all = vertcat (subk_all{:});  % all possible row    indices
      subl_all = vertcat (subl_all{:});  % all possible column indices
      
      % Shift lower triangular 'ltri' and cone 'subj_all' indices.
      offset_n = 0;
      offset_N = 0;
      for i = 1:length (obj.K.s)
        N = obj.K.s(i);
        idx = (offset_n + 1):(offset_n + N*(N+1)/2);
        ltri(idx,1) = ltri(idx,1) + offset_N;
        subj_all(idx,1) = subj_all(idx,1) * i;
        offset_n = idx(end);
        offset_N = offset_N + N^2;
      end
      
      % Reduce 'X' to semidefinite part.
      X = X(obj.K.idx.s(1,1):end,:);
      
      % Get the non-zero entries including the indices.
      [all_idx, subi, val] = find (X(ltri,:));
      subi = subi';
      subj = subj_all(all_idx)';
      subk = subk_all(all_idx)';
      subl = subl_all(all_idx)';
      val = val';
    end
    
    
    function [ltri, subk, subl] = index (N)
      % INDEX  Gen. indices for a MOSEK lower triangular vectorized matrix.
      %
      %      ltri - lower triangular index of full Fortran vectorized matrix.
      %      subk - row    index
      %      subl - column index
      %
      %   Example for N = 3:
      %
      %             [1 4 7]        ltri = [1 2 3 5 6 9]
      %         A = [2 5 8]   ==>  subk = [1 2 3 2 3 3]
      %             [3 6 9]        subl = [1 1 1 2 2 3]
      %
      %   See also vsdp, vsdp.sindex, to_mosek_fmt.
      %
      
      n = N*(N+1)/2;
      diag_idx = cumsum ([1, N*ones(1, N-1) - (0:1:N-2)]);
      
      ltri = ones(n,1);
      ltri(diag_idx(2:end)) = 2:N;
      ltri = cumsum (ltri);
      
      subk = ones(n,1);
      subk(diag_idx(2:end)) = -N*ones(1, N-1) + (2:N);
      subk = cumsum (subk);
      
      subl = zeros(n,1);
      subl(diag_idx) = ones(N,1);
      subl = cumsum (subl);
    end
    
    
    function x = to_vsdp_fmt (obj, x)
      % TO_VSDP_FMT  Translate the MOSEK solution 'X' to VSDP format.
      %
      %   See also vsdp, vsdp.svec.
      %
      
      if (isempty (x))
        return;
      end
      
      % Offset to the semidefinite cones.
      sdp_offset = obj.K.idx.s(1,1) - 1;
      x = [zeros(sdp_offset, 1); x(:)];
      
      % Index vector for sorting.
      [~,~,~,vlidx] = vsdp.sindex (obj);
      
      S = warning ('off', 'VSDP:svec:justScale');
      x = vsdp.svec (obj, x([(1:sdp_offset)'; vlidx]), 2);
      warning (S);
      
      % Strip non-semidefinite portion.
      x(1:sdp_offset) = [];
    end
    
    
    function [f,l,q,s] = supported_cones ()
      f = true;  % free   variables.
      l = true;  % linear variables.
      q = true;  % second-order cones.
      s = true;  % semidefinite cones.
    end
    
    
    function spath = install (varargin)
      % Returns the path to the installed and usable solver.  Otherwise return
      % an empty array.  No error messages are thrown.
      %
      % By passing one or more arguments interactive installation actions
      % happen and, in case of failures, error messages are thrown.
      %
      
      sname          = 'mosek';
      is_available   = @() exist ('mosekopt', 'file') == 3;
      get_path       = @() fileparts (which ('mosekopt'));
      installer_file = [];
      do_error       = (nargin > 0);
      spath = solver.registry.generic_install (sname, is_available, ...
        get_path, installer_file, do_error);
    end
  end
end
