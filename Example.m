%% ================ FLUID FLOW EMBEDDED FINITE ELEMENT ====================
%
% Author: Danilo Cavalcanti
%
% Date: April 27th, 2023.
%
%% ========================================================================
%
% Initialize workspace
initWorkspace; 
%
%% ========================== MODEL CREATION ==============================
% --- Mesh of continuum elements ------------------------------------------

% Nodes' coordinates (mm)
NODE = [0.0   0.0;
        2.0   0.0;
        2.0   2.0;
        0.0   2.0];

% Type of elements
type = 'ISOQ4';

% Element's connectivity
ELEM = [1 2 3 4];

% Thickness (mm)
t = 1.0;

% --- Mesh of the fracture elements ---------------------------------------

% Coordinates of the nodes that define the discontinuities (mm)
NODE_D = [0.0  1.5;
          2.0  1.5];

% Fractures definition (by segments)
FRACT = [1 2];

% --- Material properties of the domain -----------------------------------

% Define the material model of the continuum: 'saturated'
matModel = 'saturated';

% Material parameters
K   = 1.1574e-05;   % Hydraulic permeability (m/s)
mu  = 1.0e-6;       % Fluid dynamic viscosity (kPa*s)
gw  = 9.81;         % Specific weight of water (kPa/m)
mat = [K  gw  mu];  % Material parameters vector

% --- Material properties of the fracture ---------------------------------

% Define the traction constitutive law: 'interfaceFlow'

tractionLaw = 'interfaceFlow';  

% Values of the material constitutive model parameters
ct   = 1.0e-5;           % Leakoff at the top (m/(kPa*s))
cb   = 1.0e-5;           % Leakoff at the bottom (m/(kPa*s))
w    = 0.000;            % Initial aperture

% Assemble the vector with the material properties
matfract = [w, mu, ct, cb];

% --- Analysis model ------------------------------------------------------

% Type of analysis: 'Hydro'
anm = 'Hydro';

% --- Boundary conditions --------------------------------------------------

% Define supports
SUPP = zeros(size(NODE,1),1);
% SUPP([1 2 3 4],:) = [1; 1; 1; 1];
SUPP([3 4],:) = [1; 1];

% Define prescribe pressures (kPa)
PRESCDISPL = zeros(size(NODE,1),1);
% PRESCDISPL([1 2 3 4]) = [1;1;0;0];
PRESCDISPL([3 4]) = [0;0];

% Define the load conditions
LOAD = zeros(size(NODE,1),1);
LOAD([1 2]) = [0.05; 0.08];

% --- Order of the integration rule for the domain ------------------------

% Using Gauss quadrature
intOrder = 2;

%% ===================== EFEM FORMULATION SETUP ===========================

% Apply a sub-division of the domain to perform the numerical integration
subDivInt = true;

% Order of the interpolation of the jump displacement field
jumpOrder = 1;

% Level of the enrichment dof ('local' or 'global')
lvlEnrVar = 'global';

% Static condensation
staticCondensation = false;

%% ========================= PRE-PROCESSING ===============================

% Compute the matrix to identify to which element the fracture belongs
IDenr = 1; 

%% ========================= INITIALIZATION ===============================

% Create the model object
mdl = Model(NODE, ELEM, NODE_D, FRACT, t, matModel, mat, tractionLaw, ...
            matfract, anm, type, SUPP, LOAD, ...
            PRESCDISPL, intOrder, subDivInt,...
            jumpOrder, lvlEnrVar, staticCondensation, IDenr);

% Perform the basic pre-computations associated to the model (dof
% definition, etc.)
mdl.preComputations();

% Plot the mesh with the supports
mdl.plotMeshWithBC();

% Create the result object for the analysis
ndPlot  = 3;
dofPlot = 1; % 1 for X and 2 for Y
result  = ResultAnalysis(mdl.ID(ndPlot,dofPlot));

%% ========================== RUN ANALYSIS ================================

% Solve the structural analysis problem
anl = Anl_Linear(result);
anl.process(mdl);

%% ========================= CHECK THE RESULTS ============================

% Print the results in the command window
mdl.printResults();
mdl.plotDeformedMesh();
