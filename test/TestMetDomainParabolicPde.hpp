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
#include "PlaneBoundaryCondition.hpp"

//new celltypes
#include "TCellProliferativeType.hpp"
#include "MetCellProliferativeType.hpp"
#include "NeutrophilProliferativeType.hpp"
#include "FibroblastProliferativeType.hpp"
#include "BackgroundCellProliferativeType.hpp"
#include "SimpleLiverMetCellCycleModel.hpp"

//oxygen pde
#include "TumourCellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "EllipticGrowingDomainPdeModifier.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "CXCL8Pde.hpp"

//new force models
#include "RandomMotionForce.hpp"

//new cell killer
#include "TypeDependentKiller.hpp"


class TestMetDomainParabolic: public AbstractCellBasedTestSuite
{
public:

     //Set up some global variables - (not accurate cell cycle periods)
    double minTcellCycleDurationTime = 10;
    double maxTcellCycleDurationTime = 36;

    double minMetcellCycleDurationTime = 8;
    double maxMetcellCycleDurationTime = 12;

    double minNeutrophilCycleDurationTime = 10;
    double maxNeutrophilCycleDurationTime = 36;

    double minFibroblastCycleDurationTime = 10;
    double maxFibroblastCycleDurationTime = 36;

    double randPertScale=0.01;

    double initial_tumour_radius=1.5;
    int DomainSize = 27;


    void TestDomainParabolic()
    {
        HoneycombMeshGenerator generator(DomainSize, DomainSize);
        boost::shared_ptr<MutableMesh<2,2> > p_generating_mesh = generator.GetMesh();
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(*p_generating_mesh, 1);

        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TCellProliferativeType, p_T_type);
        MAKE_PTR(MetCellProliferativeType, p_met_type);
        MAKE_PTR(NeutrophilProliferativeType, p_neu_type);
        MAKE_PTR(FibroblastProliferativeType, p_fib_type);
        MAKE_PTR(BackgroundCellProliferativeType, p_Background_type);

        c_vector<double,2> centre = zero_vector<double>(2);
        centre(0)=DomainSize/2;
        centre(1)=DomainSize/2;

        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            //give each cell a cell cycle model
            SimpleLiverMetCellCycleModel* p_cc_model = new SimpleLiverMetCellCycleModel();
            p_cc_model->SetMinCellCycleDurationForTCell(minTcellCycleDurationTime); //edit some cell cycle properties
            p_cc_model->SetMaxCellCycleDurationForTCell(maxTcellCycleDurationTime);

            p_cc_model->SetMinCellCycleDurationForMetCell(minMetcellCycleDurationTime); //edit some cell cycle properties
            p_cc_model->SetMaxCellCycleDurationForMetCell(maxMetcellCycleDurationTime);

            p_cc_model->SetMinCellCycleDurationForNeutrophil(minNeutrophilCycleDurationTime); //edit some cell cycle properties
            p_cc_model->SetMaxCellCycleDurationForNeutrophil(maxNeutrophilCycleDurationTime);

            p_cc_model->SetMinCellCycleDurationForFibroblast(minFibroblastCycleDurationTime); //edit some cell cycle properties
            p_cc_model->SetMaxCellCycleDurationForFibroblast(maxFibroblastCycleDurationTime);

            p_cc_model->SetDimension(2);
            CellPtr p_cell(new Cell(p_state,p_cc_model));

            c_vector<double, 2> this_node_location = mesh.GetNode(i)->rGetLocation();
            double dist_to_centre = norm_2(this_node_location-centre);
            double cxcl8_initial_value = 0;
            if (dist_to_centre<initial_tumour_radius)
            {
                p_cell->SetCellProliferativeType(p_met_type);
                double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
                p_cell->SetBirthTime(birth_time);
                cxcl8_initial_value = 1;
            }
            else if (dist_to_centre>4*initial_tumour_radius && dist_to_centre<=6*initial_tumour_radius)
            {   
                double random_number=RandomNumberGenerator::Instance()->ranf();
                if (random_number<0.2)
                { //include some immune cells separated from the initial tumour mass
                    p_cell->SetCellProliferativeType(p_neu_type);
                    double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
                    p_cell->SetBirthTime(birth_time);
                }
                else if (random_number>=0.2 && random_number<0.4)
                {
                    p_cell->SetCellProliferativeType(p_T_type);
                    double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
                    p_cell->SetBirthTime(birth_time); 
                }
                else
                {
                    p_cell->SetCellProliferativeType(p_Background_type);
                    double birth_time = 0;
                    p_cell->SetBirthTime(birth_time);   
                }
            }
            else
            {
                p_cell->SetCellProliferativeType(p_Background_type);
                double birth_time = 0; // this doesn't matter because these cells don't proliferate
                p_cell->SetBirthTime(birth_time);
            }
            p_cell->GetCellData()->SetItem("CXCL8", cxcl8_initial_value);
            cells.push_back(p_cell);
        }

        NodeBasedCellPopulation<2> cell_population(mesh, cells);

        //add writer to visualise some cell properties
        cell_population.AddCellPopulationCountWriter<CellProliferativeTypesCountWriter>();
        cell_population.AddCellWriter<CellIdWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();

        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestMetDomainParabolicPde");
        simulator.SetSamplingTimestepMultiple(10);
        simulator.SetEndTime(10.0);

        MAKE_PTR(GeneralisedLinearSpringForce<2>, p_force);
        p_force->SetCutOffLength(1.5);
        simulator.AddForce(p_force);

        //add bounding box to simulation - no upper bound
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);
        point(0) = DomainSize;
        normal(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);
        point(0) = 0.0;
        point(1) = 0.0;
        normal(0) = 0.0;
        normal(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);
        point(1) = DomainSize;
        normal(1) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);

        //TODO: add an force 
        //MAKE_PTR(RandomMotionForce<2>, p_random_force);
        //p_random_force->SetMovementParameter(randPertScale);
        //simulator.AddForce(p_random_force);

        //add a pde - MODIFY THIS PDE 
        MAKE_PTR_ARGS(CXCL8Pde<2>, p_pde, (cell_population,
                                        /*DuDt*/ 1.0,
                                        /*Diffusion*/ 1,
                                        /*Decay*/ 0.3,
                                        /*Source*/ 7.5,
                                        /*Consumption*/ 1));
        MAKE_PTR_ARGS(ConstBoundaryCondition<2>, p_bc, (0.0));
        bool is_neumann_bc = true;
        ChastePoint<2> lower(0, 0, 0);
        ChastePoint<2> upper(DomainSize, DomainSize, 0);
        boost::shared_ptr<ChasteCuboid<2>> p_domain = boost::make_shared<ChasteCuboid<2> >(lower, upper);
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<2>, p_pde_modifier,
                  (p_pde, p_bc, is_neumann_bc, p_domain, 1.0));
        p_pde_modifier->SetDependentVariableName("CXCL8");
        p_pde_modifier->SetBcsOnBoxBoundary(true);
        simulator.AddSimulationModifier(p_pde_modifier);

        //add a cell killer - MODIFY THIS KILLER
        MAKE_PTR_ARGS(TypeDependentKiller<2>, p_killer,(&cell_population));
        simulator.AddCellKiller(p_killer); 

        simulator.Solve();
    }
};