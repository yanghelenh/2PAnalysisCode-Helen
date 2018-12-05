% nestedEvalin.m
%
% Helper function that enables evalin function call using variables in the
%  function that calls nestedEvalin instead of the function that called the
%  function that contains evalin. I.e. uses the workspace of the present
%  function.
%
% Notes: so far, only tested when str is variable name, not when str is a
%  complex MATLAB expression. Should probably work, but depending on how
%  evalin outputs, may only return some outputs of the str expression.
%
% INPUTS:
%   str - the string for the expression to evaluate
%
% OUTPUTS:
%   val - the output of the evaluated expression
%
% CREATED: 12/5/18 HHY
% UPDATED: 12/5/18 HHY

function val = nestedEvalin(str)
    val = evalin('caller', str);
end