#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "RandomNumberGenerator.hpp"
#include "SmartPointers.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"

#include "OffLatticeSimulation.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestChaste : public AbstractCellBasedTestSuite
{
public:
    void TestSpheroidTutorial()
    {
        EXIT_IF_PARALLEL;

        HoneycombMeshGenerator generator(10, 10, 0);
        boost::shared_ptr<MutableMesh<2,2> > p_mesh = generator.GetMesh();

        std::vector<CellPtr> cells;

        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);

        for (unsigned i=0; i<p_mesh->GetNumNodes(); i++)
        {
            SimpleOxygenBasedCellCycleModel* p_model = new SimpleOxygenBasedCellCycleModel;
            p_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state, p_model));
            p_cell->SetCellProliferativeType(p_stem_type);

            p_model->SetStemCellG1Duration(8.0);
            p_model->SetTransitCellG1Duration(8.0);

            double birth_time = - RandomNumberGenerator::Instance()->ranf() *
                                 (  p_model->GetStemCellG1Duration()
                                  + p_model->GetSG2MDuration() );
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        MeshBasedCellPopulation<2> cell_population(*p_mesh, cells);

        MAKE_PTR_ARGS(CellwiseSourceEllipticPde<2>, p_pde, (cell_population, -0.03));

        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (1.0));
        bool is_neumann_bc = false;

        MAKE_PTR_ARGS(EllipticGrowingDomainPdeModifier<2>, p_pde_modifier, (p_pde, p_bc, is_neumann_bc));
        p_pde_modifier->SetDependentVariableName("oxygen");

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.AddSimulationModifier(p_pde_modifier);

        simulator.SetOutputDirectory("SpheroidTutorial");
        simulator.SetEndTime(2.0);
        simulator.SetSamplingTimestepMultiple(10);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_linear_force);
        p_linear_force->SetCutOffLength(3);
        simulator.AddForce(p_linear_force);

        simulator.Solve();
    }
};