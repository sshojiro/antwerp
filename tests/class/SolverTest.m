classdef SolverTest < matlab.unittest.TestCase
    % SolverTest tests solutions to the quadratic equation
    % a*x^2 + b*x + c = 0
    
    methods (Test)
        function testRealSolution(testCase)
            actSolution = quadraticSolver(1,-3,2);
            % actual solution,
            % which is the implementation result
            expSolution = [2,1];% expected solution
            testCase.verifyEqual(actSolution,expSolution);
        end
        function testImaginarySolution(testCase)
            actSolution = quadraticSolver(1,2,10);% actual solution
            expSolution = [-1+3i, -1-3i];% expected solution
            testCase.verifyEqual(actSolution,expSolution);
        end
    end
    
end