classdef vsdp_indexable < handle
  % VSDP_INDEXABLE  A quantity that can be seperated
  %   Detailed explanation goes here
  
  properties
    value
    vsdp_obj
  end
  
  methods
    function obj = vsdp_indexable (value, vsdp_obj)
      % VSDP_INDEXABLE  Constructor for an indexable object.
      %
      %   obj = vsdp_indexable(value, vsdp_obj)  where 'value' is an arbitrary
      %           quantity, that should be indexed and 'vsdp_obj' is a
      %           reference to a VSDP object.
      %
      narginchk (2, 2);
      if (~isa (vsdp_obj, 'vsdp'))
        error ('VSDP_INDEXABLE:vsdp_indexable:badInput', ...
          'vsdp_indexable: Second argument must be a VSDP object.');
      end
      obj.value = value;
      obj.vsdp_obj = vsdp_obj;
    end
    
    function subval = subsref(obj, S)
      % SUBSREF  Redefine subscripted reference for objects.
      %
      %   Detailed explanation goes here
      switch (S(1).type)
        case '.'
          if (strcmp (S(1).subs, 'value'))
            if (length (S) == 1)
              subval = obj.value;
            else
              subval = subsref (obj.value, S(2:end));
            end
            return;
          end
          idx = getfield (obj.vsdp_obj.K.idx, S(1).subs);
          if (length (S) == 1)
            % Like 'obj.f'
            subval = obj.value(idx(1):idx(end),:);
          elseif ((length (S) == 2) && (isequal (S(2).type, '()')) ...
              && isscalar (S(2).subs))
            cone_idx = S(2).subs{1};
            if (isnumeric (cone_idx) && isscalar (cone_idx))
              % Like 'obj.s(2)''
              subval = obj.value(idx(cone_idx,1):idx(cone_idx,2),:);
            else
              error ('VSDP_INDEXABLE:subsref:unsupportet', ...
                'subsref: Unsupported indexing.');
            end
          else
            error ('VSDP_INDEXABLE:subsref:unsupportet', ...
              'subsref: Unsupported indexing.');
          end
        case '()'  % Like obj(1:end)
          % Simply forward them to value
          subval = subsref (obj.value, S);
        otherwise
          error ('VSDP_INDEXABLE:subsref:unsupportet', ...
            'subsref: Unsupported indexing.');
      end
    end
    
    
    function obj = subsasgn (obj, S, varargin)
      % SUBASGN  Redefine subscripted assignment for objects.
      %
      %   Detailed explanation goes here
      narginchk (3, 3);
      val = varargin{1};
      switch (S(1).type)
        case '.'
          idx = getfield (obj.vsdp_obj.K.idx, S(1).subs);
          if (length (S) == 1)
            % Like 'obj.f'
            obj.value(idx(1):idx(end),:) = val;
          elseif ((length (S) == 2) && (isequal (S(2).type, '()')) ...
              && isscalar (S(2).subs))
            sidx = S(2).subs{:};
            if (isnumeric (sidx) &&  isscalar (sidx))
              % Like 'obj.s(2)'
              obj.value(idx(sidx,1):idx(sidx,2),:) = val;
            elseif (isnumeric (sidx) &&  isvector (sidx))
              % Like 'obj.s([2,3])' or 'obj.s(2:4)'
              obj.value(idx(sidx(1),1):idx(sidx(2),2),:) = val;
            elseif (isequal (sidx, ':'))
              % Like 'obj.s(:)'
              obj.value(idx(1):idx(end),:) = val;
            else
              error ('VSDP_INDEXABLE:subsref:unsupportet', ...
                'subsref: Unsupported indexing.');
            end
          else
            error ('VSDP_INDEXABLE:subsref:unsupportet', ...
              'subsref: Unsupported indexing.');
          end
        case '()'  % Like obj(1:end)
          % Simply forward them to value
          obj.value = subsasgn (obj.value, S, val);
        otherwise
          error ('VSDP_INDEXABLE:subsref:unsupportet', ...
            'subsref: Unsupported indexing.');
      end
    end
  end
end
