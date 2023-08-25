classdef AutomaticityIndexHelper
    methods(Static)
        function b = findMatchingb(top, bottom, a, x)
            %given the top and bottom of a sigmoid, and desired x value to reach 0.9
            %of plateau, find the corresponding slope, a.
            plateauScale = 0.1;
            plateauShift = (top-bottom) * plateauScale;
            yPlateau = top - plateauShift;
            % a = (b - log(1/yMax -1)) / x;
            b = log((top-bottom)/(yPlateau - bottom)-1) + a*x;
        end
        
        function y = sigmoidValue(x, a, b, top, bottom)
            %find the ys for all given xs following a sigmoid function
            %a is the slope, b is the inflection point where function = 0.5
            %or halfway between top and bottom
            if nargin < 4
                top = 1; bottom = 0;
            end
            y = bottom + (top - bottom)./(1+exp(-a*x+a*b));
        end
        
        function x = inverseSigmoidValue(y, a, b, top, bottom)
            %given the top and bottom, slope and b for a sigmoid function,
            %find the x to reach a given y
            %a is the slope, b is the inflection point where function = 0.5
            %or halfway between top and bottom
            if nargin < 4
                top = 1; bottom = 0;
            end
            x = -(log((top-bottom) / (y - bottom) -1) - a*b)/a;
        end
        
        function val = generateRandValInRange(bottom, top)
            %generate random value in (botom, top), exclusive both ends
            val = (top-bottom).*rand(1,1) + bottom;
        end
        
        function y = linearValue(x, a, b)
            %find the ys for all given xs following y = ax+b
            y = a*x + b;
        end
    end
end