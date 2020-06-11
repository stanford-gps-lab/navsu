function plotOutput(obj,outputs,varargin)

% Just immediately call the other method :)
obj.plotOutputPpp(outputs,varargin{:})

residsData = [outputs.resids]';
residsFull = cat(1,residsData.resids);
measIdFull = cat(1,residsData.measId);
epochsFull = cat(1,residsData.epochs);

residsType = cat(1,measIdFull.TypeID);


end

