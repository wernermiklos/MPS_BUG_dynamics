function [Model] = MPO_setupLogical(Model,options)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  arguments
    Model struct;
    options.DELETE_INACTIVE logical = true; 
  end
  CL = Model.CL;
  for pos = 1:CL
    Model.MPOlogical{pos} = cellfun(@(x) ~isempty(x),Model.MPO{pos});
  end
  if ~options.DELETE_INACTIVE
    for pos = 1:CL
      [Model.MPOtasktable{pos}(:,1), Model.MPOtasktable{pos}(:,2)] = find(Model.MPOlogical{pos});
    end
    return;
  else
    something_erased = true;
    while something_erased
      something_erased = false;
      %left_to_right erase sweep
      for pos = 1:CL-1
        MPOlogical_tmp = Model.MPOlogical{pos};
        keepops = any(MPOlogical_tmp,1);
        if ~all(keepops)
          something_erased = true;
          MPOlogical_tmp = MPOlogical_tmp(:,keepops);
          Model.MPOlogical{pos} = MPOlogical_tmp;
          Model.MPO{pos} = Model.MPO{pos}(:,keepops);
          Model.MPO{pos+1} = Model.MPO{pos+1}(keepops,:);
          Model.MPOlogical{pos+1} = Model.MPOlogical{pos+1}(keepops,:);
        end
      end
      %right_to_left erase sweep
      for pos = CL:-1:2
        MPOlogical_tmp = Model.MPOlogical{pos};
        keepops = any(MPOlogical_tmp,2);
        if ~all(keepops)
          something_erased = true;
          MPOlogical_tmp = MPOlogical_tmp(keepops,:);
          Model.MPOlogical{pos} = MPOlogical_tmp;
          Model.MPO{pos} = Model.MPO{pos}(keepops,:);
          Model.MPO{pos-1} = Model.MPO{pos-1}(:,keepops);
          Model.MPOlogical{pos-1} = Model.MPOlogical{pos-1}(:,keepops);
        end
      end
    end
    for pos = 1:CL
      [Model.MPOtasktable{pos}(:,1), Model.MPOtasktable{pos}(:,2)] = find(Model.MPOlogical{pos});
    end
  end
end

