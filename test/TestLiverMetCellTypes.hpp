#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "NagaiHondaForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "WildTypeCellMutationState.hpp"
#include "CellAgesWriter.hpp"
#include "CellIdWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellVolumesWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "SmartPointers.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "DeltaNotchSrnModel.hpp"
#include "DeltaNotchTrackingModifier.hpp"

//new celltypes
#include "TCellProliferativeType.hpp"
#include "BackgroundCellProliferativeType.hpp"
#include "SimpleLiverMetCellCycleModel.hpp"

class TestLiverMetCellTypes: public AbstractCellBasedTestSuite
{
public:

    //Set up some global variables
    double minTcellCycleDurationTime = 5;
    double maxTcellCycleDurationTime = 12;

    void TestNodeBasedMonolayerWithNewTypes()
    {
        EXIT_IF_PARALLEL;

        HoneycombMeshGenerator generator(20, 20);
        boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TCellProliferativeType, p_T_type);
        MAKE_PTR(BackgroundCellProliferativeType, p_Background_type);

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            //give each cell a cell cycle model
            SimpleLiverMetCellCycleModel* p_cc_model = new SimpleLiverMetCellCycleModel();
            p_cc_model->SetMinCellCycleDurationForTCell(minTcellCycleDurationTime); //edit some cell cycle properties
            p_cc_model->SetMaxCellCycleDurationForTCell(maxTcellCycleDurationTime);
            p_cc_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state,p_cc_model));

            //give each cell a proliferative type (randomly assign)
            if (RandomNumberGenerator::Instance()->ranf()<0.3)
            {
                p_cell->SetCellProliferativeType(p_T_type);
                double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
                p_cell->SetBirthTime(birth_time);

            }
            else
            {
                p_cell->SetCellProliferativeType(p_Background_type);
                double birth_time = 0; // this doesn't matter because these cells don't proliferate
                p_cell->SetBirthTime(birth_time);
            }
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        //add writer to visualise some cell properties
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestNodeBasedMonolayerWithNewTypes");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(20.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        simulator.Solve();
    }
    
};