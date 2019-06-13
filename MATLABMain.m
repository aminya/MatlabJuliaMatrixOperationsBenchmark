function [  ] = MATLABMain( operationMode )

FALSE   = 0;
TRUE    = 1;
OFF     = 0;
ON      = 1;

OPERATION_MODE_PARTIAL  = 1; %<! For Testing (Runs Fast)
OPERATION_MODE_FULL     = 2;

if(exist('operationMode', 'var') == FALSE)
    operationMode = OPERATION_MODE_FULL;
end

mRunTime = MatlabMatrixBenchmark0001(operationMode);
mRunTime = MatlabMatrixBenchmark0002(operationMode);
mRunTime = MatlabMatrixBenchmark0003(operationMode);


end

