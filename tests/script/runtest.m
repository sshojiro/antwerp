quadTests = matlab.unittest.TestSuite.fromClass(?SolverTest);

result = run(quadTests);
% Running SolverTest
% ..
% Done SolverTest
% __________

whos result
%   Name        Size            Bytes  Class                         Attributes
% 
%   result      1x2               250  matlab.unittest.TestResult              

result(1)
%   TestResult with properties:
% 
%           Name: 'SolverTest/testRealSolution'
%         Passed: 1
%         Failed: 0
%     Incomplete: 0
%       Duration: 0.6470
% 
% Totals:
%    1 Passed, 0 Failed, 0 Incomplete.
%    0.64695 seconds testing time.

rt = table(result)
% rt = 
% 
%                    Name                   Passed    Failed    Incomplete    Duration
%     __________________________________    ______    ______    __________    ________
% 
%     'SolverTest/testRealSolution'         true      false     false          0.64695
%     'SolverTest/testImaginarySolution'    true      false     false         0.080053

sortrows(rt,'Duration')
%                    Name                   Passed    Failed    Incomplete    Duration
%     __________________________________    ______    ______    __________    ________
% 
%     'SolverTest/testImaginarySolution'    true      false     false         0.080053
%     'SolverTest/testRealSolution'         true      false     false          0.64695

writetable(rt,'myTestResults.csv','QuoteStrings',true)
% create .csv