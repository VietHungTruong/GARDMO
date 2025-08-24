
classdef GARDMO < ALGORITHM
% <2025> <multi> <real/constrained>
% Genetic Algorithm with Adaptive Random Differential Mutation Operator

    methods
        function main(Algorithm, Problem)
            %% Parameter settings
            ArchiveSize = 500;  % Archive size
            APRDMOSize = 25;    % Population resize
            AllFeasible = false;  % Flag all_feasible
            PopSize = Problem.N;
            %% Generate initial population
            Dec = lhsdesign(PopSize, Problem.D) .* repmat(Problem.upper - Problem.lower, PopSize, 1) + repmat(Problem.lower, PopSize, 1);
            Population = Problem.Evaluation(Dec);
            Archive = Population(1:0);  % Initial archive empty as in Python

            %% Optimization
            while Algorithm.NotTerminated(Archive)
                % Tính fronts cho population
                Fronts = NDSort(Population.objs, Population.cons, inf);
                Rank0 = Population(Fronts == 1);  % Non-dominated front
                RankWorst = Population(Fronts == max(Fronts));  % Last front

                % Generate offspring
                OffDec = zeros(length(Population), Problem.D);
                for i = 1:length(Population)
                    Current = Population(i);

                    % Selection parents
                    if ~isempty(Rank0)
                        Parent1Idx = randi(length(Rank0));
                        Parent1 = Rank0(Parent1Idx);
                    else
                        Parent1 = Population(randi(length(Population)));
                    end                        

                    if ~isempty(RankWorst)
                        Parent2Idx = randi(length(RankWorst));
                        Parent2 = RankWorst(Parent2Idx);
                    else
                        Parent2 = Population(randi(length(Population)));
                    end

                    Parents = [Parent1, Parent2];

                    % Crossover with adaptive probability
                    Progress = Problem.FE / Problem.maxFE;
                    NSGA2Prob = 0.5 + 0.3 * Progress;
                    if rand < NSGA2Prob
                        % SBX crossover on parents
                        ChildDec = SBXCrossover(Parents(1).dec, Parents(2).dec, Problem.lower, Problem.upper, 1.0, 20);
                    else
                        % DE-like crossover on current
                        aaa = randi(Problem.D);
                        ChildDec = Current.dec;
                        for j = 1:Problem.D
                            if rand < rand || j == aaa
                                ChildDec(j) = ChildDec(j) + rand * (Parents(1).dec(j) - Parents(2).dec(j));
                            end
                        end
                        ChildDec = max(min(ChildDec, Problem.upper), Problem.lower);
                    end

                    % Mutation (Polynomial mutation)
                    ChildDec = PolynomialMutation(ChildDec, Problem.lower, Problem.upper, 1/Problem.D, 20);

                    OffDec(i, :) = ChildDec;
                end

                Offspring = Problem.Evaluation(OffDec);

                % Combine all: population + offspring + archive
                AllSols = [Population, Offspring, Archive];

                
                % Adjust population size
                FeasibleCount = sum(sum(max(0, AllSols.cons), 2) == 0);
                FeasibleRatio = FeasibleCount / length(AllSols);
                if ~AllFeasible && FeasibleRatio > 0.8
                    AllFeasible = true;
                end
                if AllFeasible
                    if Problem.FE > Problem.N*10
                        PopSize = APRDMOSize;
                    end
                else
                    PopSize = Problem.N;
                end

                % Non-dominated sorting and crowding for new population
                [Population, FrontNo, CrowdDis] = EnvironmentalSelection(AllSols, PopSize);

                % For archive: select feasible solutions
                CV = sum(max(0, AllSols.cons), 2);
                FeasibleSols = AllSols(CV == 0);
                if ~isempty(FeasibleSols)
                    [Archive, ~, ~] = EnvironmentalSelection(FeasibleSols, ArchiveSize);
                else
                    Archive = AllSols(1:min(length(AllSols), ArchiveSize));
                end
            end

            % Final filter: non-dominated front of archive
            Fronts = NDSort(Archive.objs, Archive.cons, inf);
            Archive = Archive(Fronts == 1);
        end
    end
end

function ChildDec = SBXCrossover(Parent1, Parent2, Lower, Upper, ProC, DisC)
    D = length(Parent1);
    Mu = rand(1, D);
    Beta = zeros(1, D);
    Beta(Mu <= 0.5) = (2 * Mu(Mu <= 0.5)) .^ (1 / (DisC + 1));
    Beta(Mu > 0.5) = (2 - 2 * Mu(Mu > 0.5)) .^ (-1 / (DisC + 1));
    Beta = Beta .* ((-1) .^ randi([0 1], 1, D));
    Beta(rand(1, D) < 0.5) = 1;
    Beta(rand(1, D) > ProC) = 1;
    ChildDec = (Parent1 + Parent2) / 2 + Beta .* (Parent1 - Parent2) / 2;
    ChildDec = max(min(ChildDec, Upper), Lower);
end

function ChildDec = PolynomialMutation(ChildDec, Lower, Upper, ProM, DisM)
    D = length(ChildDec);
    Site = rand(1, D) < ProM;
    Mu = rand(1, D);
    Temp = Site & (Mu <= 0.5);
    Delta = (ChildDec(Temp) - Lower(Temp)) ./ (Upper(Temp) - Lower(Temp));
    ChildDec(Temp) = ChildDec(Temp) + (Upper(Temp) - Lower(Temp)) .* ...
        ((2 * Mu(Temp) + (1 - 2 * Mu(Temp)) .* (1 - Delta) .^ (DisM + 1)) .^ (1 / (DisM + 1)) - 1);
    Temp = Site & (Mu > 0.5);
    Delta = (Upper(Temp) - ChildDec(Temp)) ./ (Upper(Temp) - Lower(Temp));
    ChildDec(Temp) = ChildDec(Temp) + (Upper(Temp) - Lower(Temp)) .* ...
        (1 - (2 * (1 - Mu(Temp)) + 2 * (Mu(Temp) - 0.5) .* (1 - Delta) .^ (DisM + 1)) .^ (1 / (DisM + 1)));
    ChildDec = max(min(ChildDec, Upper), Lower);
end