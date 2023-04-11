%% ResultAnalysis Class
%
% A result analysis object is responsible for storing the analysis results.
%
classdef ResultAnalysis < handle
    %% Public properties
    properties (SetAccess = public, GetAccess = public)
        % Result options
        name    = [];  % vector of curves names
        node    = [];  % vector of curves nodes
        dof     = [];  % vector of curves d.o.f's
        
        % Result storage
        steps = 0;     % number of performed steps
        niter = 0;     % number of total iterations
        lbd   = [];    % vector of load ratios of all steps
        F     = [];    % vector of nodal forces and reactions
        U     = [];    % matrix of nodal displacement vectors of all steps/modes
    end
    
    %% Constructor method
    methods
        %------------------------------------------------------------------
        function res = ResultAnalysis(dof)
            res.dof = dof;
        end
    end

    %% Public methods
    methods
        function plotCurves(this)
            figure
            hold on
            box on, grid on, axis on
            plot(this.U, this.lbd, 'o-k');
            %plot([0.0,1.0,4.5],[0.0,0.5,0.0],'-b')
            %legend('Numerical','Analytical','Interpreter','latex')
            xlabel('Displacement (mm)','Interpreter','latex')
            ylabel('Load factor','Interpreter','latex')
            xaxisproperties= get(gca, 'XAxis');
            xaxisproperties.TickLabelInterpreter = 'latex'; 
            yaxisproperties= get(gca, 'YAxis');
            yaxisproperties.TickLabelInterpreter = 'latex';   
            set(gca,'FontSize',14);
        end

    end
end