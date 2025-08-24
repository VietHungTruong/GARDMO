clear;
% Please make sure that GARDMO.m was input into PlatEMO V4.7
Algorithms = {@GARDMO};
AlgNames = {'GARDMO'};
NumRuns = 30;
SaveNum = 50; % Number populations saved per run to track process (spaced evenly)
Problems = {'DTLZ6'};
maxFEs = [20000];
M = 3;
D = 12;
N = 100;
% Perform each problem
for p = 1:length(Problems)
    ProblemName = Problems{p};
    maxFE = maxFEs(p);
    ProblemFunc = str2func(ProblemName);
    common_evals = linspace(0, maxFE, SaveNum); % Adjust to max of maxFE
    IGDProcessAll = zeros(length(Algorithms), length(common_evals)); % To store average IGD process for each alg
    AllGDs = zeros(NumRuns, length(Algorithms));
    AllIGDs = zeros(NumRuns, length(Algorithms));
    AllHVs = zeros(NumRuns, length(Algorithms));
    meanGDs = zeros(1, length(Algorithms));
    stdGDs = zeros(1, length(Algorithms));
    meanIGDs = zeros(1, length(Algorithms));
    stdIGDs = zeros(1, length(Algorithms));
    meanHVs = zeros(1, length(Algorithms));
    stdHVs = zeros(1, length(Algorithms));
    
	% Perform each algorithm
    for a = 1:length(Algorithms)
        %rng('shuffle'); % Set random seed based on current time
        for run = 1:NumRuns
            platemo('algorithm', Algorithms{a}, ...
                    'problem', ProblemFunc, ...
                    'M', M, ...
                    'N', N, ...
                    'maxFE', maxFE, ...
                    'run', run, ...
                    'save', SaveNum); % Save SaveNum populations to file
        end
    end
    
    pro = ProblemFunc('M', M, 'D', D);PF = pro.GetOptimum(2000);

    for a = 1:length(Algorithms)
        alg = Algorithms{a};
        algName = AlgNames{a};
        GDs = zeros(NumRuns,1);
        IGDs = zeros(NumRuns,1);
        HVs = zeros(NumRuns,1);
        % Now load and compute metrics        
        for run = 1:NumRuns
            % Load the saved data (correct path)
            loadFile = fullfile('Data', algName, [algName '_' ProblemName '_M' num2str(M) '_D' num2str(D) '_' num2str(run) '.mat']);
            data = load(loadFile);
            savedPop = data.result; % cell of saved populations
            % Adjust SaveNum to actual number of saved populations
            actualSaveNum = size(savedPop,1);
            % Last population for final metrics
            finalPop = savedPop{end,2};
            % Calculate metrics using PlatEMO
            gd = GD(finalPop, PF);
            igd = IGD(finalPop, PF);
            hv = HV(finalPop, PF);
            GDs(run) = gd;
            IGDs(run) = igd;
            HVs(run) = hv;
        end
        AllGDs(:,a) = GDs;
        AllIGDs(:,a) = IGDs;
        AllHVs(:,a) = HVs;
        % Calculate mean và std
        meanGD = mean(GDs);
        stdGD = std(GDs);
        meanIGD = mean(IGDs);
        stdIGD = std(IGDs);
        meanHV = mean(HVs);
        stdHV = std(HVs);
        meanGDs(a) = meanGD
        stdGDs(a) = stdGD
        meanIGDs(a) = meanIGD
        stdIGDs(a) = stdIGD
        meanHVs(a) = meanHV
        stdHVs(a) = stdHV
    end 
end