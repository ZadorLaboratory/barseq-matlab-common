function hsvcolmap=hsvcolmapQY(n,paramstruct)
% hsvcolorQY provides color and changes hue gradually.
%
%
%
% Inputs:
%   n      - The total number of colors you want to generate.
%
%   paramstruct - a Matlab structure of input parameters
%                    Use: "help struct" and "help datatypes" to
%                         learn about these.
%                    Create one, using commands of the form:
%
%       paramstruct = struct('field1',values1,...
%                            'field2',values2,...
%                            'field3',values3) ;
%
%                          where any of the following can be used,
%                          these are optional, misspecified values
%                          revert to defaults
%
%                    Version for easy copying and modification:
%    paramstruct = struct('',, ...
%                         '',, ...
%                         '',) ;
%
%    fields            values
%
%    saturation       A vector of reals numbers in the interval [0,1]. 
%                        0 a column vector: one saturation for each data
%                        point.
%                        1 (default) a column vector of 1s: all points' saturation to be 1.   
%                           (point's saturation being 1 means no saturation of white color and 
%                              point's saturation being 1 means white color. )
%
%    value              A vector of reals numbers in the interval [0,1].
%                         (value in hsv color means the brightness of the
%                         color: 0- darkest, 1- brightest)
%                        0 a column vector: one value for each data point
%                        1 (default) a column vector of 1s: all points' value to be 1. 
%
% Output:
%
%    hsvcolmap   - The transformed hsv colormap from magenta to red for the data.
%
% Copyright (c) Qunqun Yu 2014

% First, set saturation and value to defaults
saturation=ones(n,1);
value=ones(n,1);
%  Now update parameters as specified,
%  by parameter structure (if it is used)
%
if nargin > 1 ;   %  then paramstruct is an argument

  if isfield(paramstruct,'saturation') ;    %  then change to input value
    saturation = getfield(paramstruct,'saturation') ; 
  end ;

  if isfield(paramstruct,'value') ;    %  then change to input value
    value = getfield(paramstruct,'value') ; 
  end ;

end
hue=[0:(5/6)/(n-1):5/6];
 % Get equally spaced hue of the color
hsvcolmap = hsv2rgb([flipud(hue'), saturation, value] );

