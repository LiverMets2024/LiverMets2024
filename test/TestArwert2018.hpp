//
// Created by bull on 20/10/20.
//

#ifndef JB_CHASTE_TESTARWERT2018_HPP
#define JB_CHASTE_TESTARWERT2018_HPP

#include <cxxtest/TestSuite.h>
/* Most Chaste code uses PETSc to solve linear algebra problems.  This involves starting PETSc at the beginning of a test-suite
 * and closing it at the end.  (If you never run code in parallel then it is safe to replace PetscSetupAndFinalize.hpp with FakePetscSetup.hpp)
 */

#include "CheckpointArchiveTypes.hpp"
#include "CellBasedSimulationArchiver.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CellsGenerator.hpp"
#include "TransitCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "UniformCellCycleModel.hpp"
#include "HoneycombMeshGenerator.hpp"
#include "GeneralisedLinearSpringForce.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
/* The next header file defines the class for storing the spatial information of cells. */
#include "NodesOnlyMesh.hpp"
/* The next header file defines a node-based {{{CellPopulation}}} class.*/
#include "NodeBasedCellPopulation.hpp"

#include "SimpleRadiationKiller.hpp"

// Write VTU files
#include "NodeLocationWriter.hpp"
#include "NodeVelocityWriter.hpp"
#include "BoundaryNodeWriter.hpp"
#include "CellCycleClockWriter.hpp"

//#include "SimpleOxygenBasedCellCycleModelFreezeWhenHypoxic.hpp"
//#include "UniformCellCycleModelWithQuiescence.hpp"
#include "NoCellCycleModel.hpp"
#include "MacrophagePhenotypeSwitchingCellCycle.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
//#include "MacrophageCellProliferativeType.hpp"
//#include "ChemotacticForceCSF1.hpp"
#include "CellwiseSourceEllipticPde.hpp"
#include "ConstBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
#include "Debug.hpp"

#include "CellwiseSourceParabolicPdeWithDecay.hpp"
#include "CellwiseSourceParabolicPde.hpp"
#include "ParabolicBoxDomainPdeModifier.hpp"
#include "AveragedSourceParabolicPde.hpp"
#include "UniformSourceParabolicPde.hpp"
#include "AveragedSourceParabolicPde_edited.hpp"
#include "AveragedSourceParabolicPde_edited_phenotypeDependent.hpp"

#include "AveragedSourceEllipticPde.hpp"
#include "UniformSourceEllipticPde.hpp"
#include "EllipticBoxDomainPdeModifier.hpp"
#include "EllipticBoxDomainPdeModifierVariableTimestep.hpp"
#include "EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs.hpp"
#include "EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs.hpp"
#include "ApoptoticCellKiller.hpp"
#include "BoundaryConditionsContainer.hpp"

//#include "DiffusionForceChooseD.hpp"
//#include "BiasedBrownianMotion.hpp"

#include "CellMutationStatesCountWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"

#include "CellVolumesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellMutationStatesWriter.hpp"

#include "GeneralisedLinearSpringForce.hpp"//DifferentialAdhesionForApoptosisAndMacrophages.hpp"
#include "GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages.hpp"
#include "ExternalPressureForceOnConcaveHull.hpp"

#include "CellPropertyCollection.hpp"
#include "CellLabel.hpp"

#include "BoundaryCellWriter.hpp"

#include "PlaneBoundaryCondition.hpp"

#include "SpheroidComparisonCellCycle.hpp"
#include "ContactInhibitionCellCycleModel.hpp"

#include "PhenotypeBasedMacrophagePhagocytosisModifier.hpp"
#include "ChemotacticForce_SpecifyNutrientAndCellType.hpp"
#include "ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType.hpp"

#include "Arwert2018_MacrophageCellCycle.hpp"
#include "GreenspanAndVolumeCellCycle.hpp"
#include "Arwert2018_AddMacrophagesModifier.hpp"
#include "Arwert2018_AddMacrophagesFromPointVesselModifier.hpp"
#include "Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier.hpp"
#include "EgfBasedMacrophageTumourCellBindingModifier.hpp"

#include "DiffusionForceChooseD.hpp"

#include "ApplyRadiationProtocolModifier.hpp"
#include "VolumeTrackingModifier.hpp"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/lexical_cast.hpp>



class TestArwert2018 : public CxxTest::TestSuite
{
public:

    void dontTestStromaOnly() throw(Exception)
    {
        /*
         * Replicate the system described in Arwert et al (2018): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5946803/pdf/main.pdf
         *
         * Macrophages "extravasate" with phenotype 0 due to CSF1
         * They're recruited up the CSF1 gradient into the tumour
         * Tumour expression of TGF-b causes their phenotype to change towards 1
         * (Interpreted as higher expression of CXCR4)
         * Macrophages with phenotype close to 1 are more sensitive to CXCL12 expressed by perivascular fibroblasts
         * Hence macrophages move back towards the vessel while also expressing EGF (hence attracting tumour cells)
         *
         * Change is unidirectional - i.e., phenotype can increase towards 1, but not away from it
         */

        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        int visualisationOutputFrequencyPerHour = 2;
        const int DIM = 2;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<DIM>*> nodes;
        // Add some nodes

        MARK;

        // Domain size
        double height = 15;
        double width = 15;

        // Initial cell fill size
        double initialCellsHeight = 15;
        double initialStromalCellsWidth = 15;
//        double initialTumourCellsWidth = 5;
        // Stroma Nodes
        unsigned nodeNum=0;
        unsigned numStromaNodes = 0;
        for (double x=0; x<=initialStromalCellsWidth; x++)
        {
            for (double y=0; y<=initialCellsHeight; y++) {
                if (DIM == 2) {
                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                    nodeNum++;
                    numStromaNodes++;
                }
            }

        }
//        // Tumour Nodes
//        unsigned numTumourNodes = 0;
//        for (double x=initialStromalCellsWidth+1; x<=initialStromalCellsWidth+initialTumourCellsWidth; x++)
//        {
//            for (double y=0; y<=initialCellsHeight; y++) {
//                if (DIM == 2) {
//                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
//                    nodeNum++;
//                    numTumourNodes++;
//                }
//            }
//
//        }


        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

//        double avgCellCycleDurationTumour = 16;
        double avgCellCycleDurationStroma = 32;
//        unsigned cancerLabelColour = 3;
//        unsigned macrophageLabelColour = 7;
        MAKE_PTR_ARGS(CellLabel, p_stroma_label, (1));
//        MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));
        for (unsigned i=0; i<numStromaNodes; i++)
        {
            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationStroma*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationStroma*1.25);
            p_model->SetHypoxicConcentration(0.2);
            p_model->SetNecroticConcentration(0.05);
            p_model->SetQuiescenceVolumeProportion(0.75);

//            ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
//            p_model->SetDimension(DIM);
//            p_model->SetQuiescentVolumeFraction(0.7);
//            p_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
//            p_cell->GetCellData()->SetItem("csf1", 0);
//            p_cell->GetCellData()->SetItem("egf", 0);
//            p_cell->GetCellData()->SetItem("cxcl12", 0);
//            p_cell->GetCellData()->SetItem("tgf", 0);
//            p_cell->GetCellData()->SetItem("phenotype", -2);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_stroma_label);
            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationStroma*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);

        }
//        for (unsigned i=numStromaNodes; i<numStromaNodes+numTumourNodes; i++)
//        {
//
//            SpheroidComparisonCellCycle* p_model = new SpheroidComparisonCellCycle;
//            p_model->SetDimension(DIM);
//            p_model->SetMinCellCycleDuration(avgCellCycleDurationTumour*0.75);
//            p_model->SetMaxCellCycleDuration(avgCellCycleDurationTumour*1.25);
//            p_model->SetHypoxicConcentration(0.15);
//            p_model->SetNecroticConcentration(0.1);
//
//            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
//            p_cell->SetCellProliferativeType(p_stem_type);
//            p_cell->GetCellData()->SetItem("oxygen", 1);
//            p_cell->GetCellData()->SetItem("csf1", 0);
//            p_cell->GetCellData()->SetItem("egf", 0);
//            p_cell->GetCellData()->SetItem("cxcl12", 0);
//            p_cell->GetCellData()->SetItem("tgf", 0);
//            p_cell->GetCellData()->SetItem("phenotype", -1);
//            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
//            p_cell->AddCellProperty(p_cancer_label);
////            MARK;
//            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
//                               (avgCellCycleDurationTumour*0.75);
//            p_cell->SetBirthTime(birthTime);
//            cells.push_back(p_cell);
//        }
//
//        // Add macrophages
//        double numMacrophagesToAdd = 10;
//        double normalize;
//        std::vector<double> randomCoords = {0,0};
//
//        MAKE_PTR_ARGS(CellLabel, p_macrophage_label, (macrophageLabelColour));
//        for (int j = 0; j < numMacrophagesToAdd; j++) {
//
//            // Make Macrophage
//            Arwert2018_MacrophageCellCycle *p_model = new Arwert2018_MacrophageCellCycle;
//            p_model->SetDimension(DIM);
//            p_model->SetPhagocytosisCooldownPeriod(24.0);
//            CellPtr pNewCell(new Cell(p_state, p_model));
//            pNewCell->SetCellProliferativeType(p_diff_type);
//            pNewCell->GetCellData()->SetItem("oxygen", 1);
//            pNewCell->GetCellData()->SetItem("csf1", 0);
//            pNewCell->GetCellData()->SetItem("egf", 0);
//            pNewCell->GetCellData()->SetItem("cxcl12", 0);
//            pNewCell->GetCellData()->SetItem("tgf", 0);
//            pNewCell->GetCellData()->SetItem("phenotype", 0);
//            pNewCell->AddCellProperty(p_macrophage_label);
//
//            Node<DIM> *p_new_node;
//            if (DIM == 2) {
//                p_new_node = new Node<DIM>(nodeNum, false, 0.5 , 2*j + 5);
//                nodes.push_back(p_new_node);
//                nodeNum++;
//            }
//            cells.push_back(pNewCell);
//        }


        NodesOnlyMesh<DIM> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);


        // Make cell population (2D)
        NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

        // Write summary statistic files
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();
        cell_population.AddCellWriter<CellCycleClockWriter>();


        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<DIM> lower(0, 0, 0);
        ChastePoint<DIM> upper(width,height, 0);
        MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));


        //
        // OXYGEN
        //
        double consumptionRate = -0.03;//was -0.03
        double diffusionCoefficient = 1;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_oxygen_pde, (cell_population, consumptionRate,diffusionCoefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_oxygen_bc, (1));
        bool is_neumann_bc = false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        int updateIntervalForPdeInTimesteps = 120/2;
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs<DIM>, p_oxygen_pde_modifier, (p_oxygen_pde, p_oxygen_bc, is_neumann_bc, p_cuboid, 1.0));
        p_oxygen_pde_modifier->SetTimestepInterval(1);
        p_oxygen_pde_modifier->SetDependentVariableName("oxygen");
        p_oxygen_pde_modifier->SetBcsOnBoxBoundary(false); //was false

//
//        //
//        // CXCL12
//        //
//        double cxcl12_decay_rate = -0.02;
//        MAKE_PTR_ARGS(UniformSourceEllipticPde<DIM>, p_cxcl12_pde, (cxcl12_decay_rate));
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_cxcl12_bc, (1));
//        bool is_neumann_bc_cxcl12= false; // Dirichlet BCs
//        // Create a PDE modifier and set the name of the dependent variable in the PDE
//        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs<DIM>, p_cxcl12_pde_modifier, (p_cxcl12_pde, p_cxcl12_bc, is_neumann_bc_cxcl12, p_cuboid, 1.0));
//        p_cxcl12_pde_modifier->SetTimestepInterval(1);
//        p_cxcl12_pde_modifier->SetDependentVariableName("cxcl12");
//        p_cxcl12_pde_modifier->SetBcsOnBoxBoundary(false); //was false
//        p_cxcl12_pde_modifier->SetOutputGradient(true);
//
//        //
//        // CSF-1
//        //
//        double dudtCoefficient = 1.0;
//        double csfDiffusionCoefficient = 0.4;
//        double sourceCoefficient = 0.2;
//        double decayRate_csf = 0.1;
//        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_csf1_pde, (cell_population, dudtCoefficient,csfDiffusionCoefficient,sourceCoefficient));
//        p_csf1_pde->SetCellTypesAsSource({cancerLabelColour});
//        p_csf1_pde->SetDecayRate(decayRate_csf);
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_csf1_bc, (0));
//        bool is_neumann_bc_csf1 = true;
//        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_csf1_pde_modifier, (p_csf1_pde, p_csf1_bc, is_neumann_bc_csf1, p_cuboid, 1.0));
//        p_csf1_pde_modifier->SetDependentVariableName("csf1");
//        p_csf1_pde_modifier->SetBcsOnBoxBoundary(true);
//        p_csf1_pde_modifier->SetOutputGradient(true);
//
//        //
//        // TGF
//        //
//        double dudtCoefficient_tgf = 1.0;
//        double diffusionCoefficient_tgf = 0.1;
//        double sourceCoefficient_tgf = 0.1;
//        double decayRate_tgf = 0.1;
//        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_tgf_pde, (cell_population, dudtCoefficient_tgf,diffusionCoefficient_tgf,sourceCoefficient_tgf));
//        p_tgf_pde->SetCellTypesAsSource({cancerLabelColour});
//        p_tgf_pde->SetDecayRate(decayRate_tgf);
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_tgf_bc, (0));
//        bool is_neumann_bc_tgf = true;
//        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_tgf_pde_modifier, (p_tgf_pde, p_tgf_bc, is_neumann_bc_tgf, p_cuboid, 1.0));
//        p_tgf_pde_modifier->SetDependentVariableName("tgf");
//        p_tgf_pde_modifier->SetBcsOnBoxBoundary(true);
//        p_tgf_pde_modifier->SetOutputGradient(true);
//
//        //
//        // EGF
//        //
//        double dudtCoefficient_egf = 1.0;
//        double diffusionCoefficient_egf = 0.2;
//        double sourceCoefficient_egf = 0.2;
//        double decayRate_egf = 0.1;
//        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_egf_pde, (cell_population, dudtCoefficient_egf,diffusionCoefficient_egf,sourceCoefficient_egf));
//        p_egf_pde->SetCellTypesAsSource({macrophageLabelColour});
//        p_egf_pde->SetDecayRate(decayRate_egf);
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_egf_bc, (0));
//        bool is_neumann_bc_egf = true;
//        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_egf_pde_modifier, (p_egf_pde, p_egf_bc, is_neumann_bc_egf, p_cuboid, 0.5));
//        p_egf_pde_modifier->SetDependentVariableName("egf");
//        p_egf_pde_modifier->SetBcsOnBoxBoundary(true);
//        p_egf_pde_modifier->SetOutputGradient(true);



        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<DIM> simulator(cell_population);
        simulator.AddSimulationModifier(p_oxygen_pde_modifier);
//        simulator.AddSimulationModifier(p_csf1_pde_modifier);
//        simulator.AddSimulationModifier(p_egf_pde_modifier);
//        simulator.AddSimulationModifier(p_tgf_pde_modifier);
//        simulator.AddSimulationModifier(p_cxcl12_pde_modifier);
        std::stringstream output_directory;
        output_directory << "Arwert2018/I_StromalCellsOnly/";
        simulator.SetOutputDirectory(output_directory.str());
        simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
        simulator.SetEndTime(200);

//        MAKE_PTR(PhenotypeBasedMacrophagePhagocytosisModifier<DIM>, p_phagoytosis_modifier);
//        p_phagoytosis_modifier->SetCellTypesToEat({cancerLabelColour});
//        simulator.AddSimulationModifier(p_phagoytosis_modifier);


        // Add Brownian motion for all cells
//		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(0.01);
//		simulator.AddForce(p_diffusion_force);

        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

//        // Macrophages go up csf1 gradient
//        // todo make this dependent on phenotype close to 0?
//        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCSF1);
//        p_chemotaxis_MacrophageToCSF1->SetNutrient("csf1");
//        p_chemotaxis_MacrophageToCSF1->SetCellTypesForChemotaxis({macrophageLabelColour});
//        p_chemotaxis_MacrophageToCSF1->SetSensitivity(3);
//        simulator.AddForce(p_chemotaxis_MacrophageToCSF1);
//
//        // Tumour cells go up egf gradient
//        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_TumourCellToEGF);
//        p_chemotaxis_TumourCellToEGF->SetNutrient("egf");
//        p_chemotaxis_TumourCellToEGF->SetCellTypesForChemotaxis({cancerLabelColour});
//        p_chemotaxis_TumourCellToEGF->SetSensitivity(5);
//        simulator.AddForce(p_chemotaxis_TumourCellToEGF);
//
//        // Macrophages go up cxcl12 gradient, depending on phenotype (close to 1)
//        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCXCL12);
//        p_chemotaxis_MacrophageToCXCL12->SetNutrient("cxcl12");
//        p_chemotaxis_MacrophageToCXCL12->SetCellTypesForChemotaxis({macrophageLabelColour});
//        p_chemotaxis_MacrophageToCXCL12->SetSensitivity(5);
//        simulator.AddForce(p_chemotaxis_MacrophageToCXCL12);


        MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(5.0);
        p_force->SetCutOffLength(1.5);
//        p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
//        p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        simulator.AddForce(p_force);
        MARK;

        // Killer which removes apoptotic cells
        MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
        simulator.AddCellKiller(p_apoptosis_killer);

        // Boundary conditions: Add a wall along the x axis
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        // Boundary conditions: Add a wall along the bottom
        c_vector<double,2> normal2 = zero_vector<double>(2);
        c_vector<double,2> point2 = zero_vector<double>(2);
        normal2(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        // Boundary conditions: Add a wall along the top
        c_vector<double,2> normal3 = zero_vector<double>(2);
        c_vector<double,2> point3 = zero_vector<double>(2);
        normal3(1) = 1.0;
        point3(1) = height;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point3, normal3));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        // Boundary conditions: Add a wall along the right hand side
        c_vector<double,2> point4 = zero_vector<double>(2);
        c_vector<double,2> normal4 = zero_vector<double>(2);
        point4(0) = width;
        normal4(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point4, normal4));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);


        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
//		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

    }


    void dontTestStromaAndTumour() throw(Exception)
    {
        /*
         * Replicate the system described in Arwert et al (2018): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5946803/pdf/main.pdf
         *
         * Macrophages "extravasate" with phenotype 0 due to CSF1
         * They're recruited up the CSF1 gradient into the tumour
         * Tumour expression of TGF-b causes their phenotype to change towards 1
         * (Interpreted as higher expression of CXCR4)
         * Macrophages with phenotype close to 1 are more sensitive to CXCL12 expressed by perivascular fibroblasts
         * Hence macrophages move back towards the vessel while also expressing EGF (hence attracting tumour cells)
         *
         * Change is unidirectional - i.e., phenotype can increase towards 1, but not away from it
         */

        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        int visualisationOutputFrequencyPerHour = 2;
        const int DIM = 2;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<DIM>*> nodes;
        // Add some nodes

        MARK;

        // Domain size
        double height = 30;
        double width = 15;

        // Initial cell fill size
        double initialCellsHeight = 30;
        double initialStromalCellsWidth = 15;
//        double initialTumourCellsWidth = 5;
        // Stroma Nodes
        unsigned nodeNum=0;
        unsigned numStromaNodes = 0;
        // For staggering rows
        unsigned rownum = 1;
        double offset = 0;
        for (double x=0; x<=initialStromalCellsWidth; x=x+0.75)//x++)
        {
            offset = 0;
            if (rownum % 2)
            {
                offset = 0.75/2;
            }
            for (double y=0+offset; y<=initialCellsHeight; y=y+0.75)//y++)
            {
                if (DIM == 2) {
                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                    nodeNum++;
                    numStromaNodes++;
                }
            }
            rownum++;
        }
        // Tumour Nodes
        unsigned numTumourNodes = 0;
        // Seed one tumour cell
        nodes.push_back(new Node<DIM>(nodeNum, false, initialStromalCellsWidth-2.5, height*0.5));
        nodeNum++;
        numTumourNodes++;
//
//        for (double x=initialStromalCellsWidth+1; x<=initialStromalCellsWidth+initialTumourCellsWidth; x++)
//        {
//            for (double y=0; y<=initialCellsHeight; y++) {
//                if (DIM == 2) {
//                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
//                    nodeNum++;
//                    numTumourNodes++;
//                }
//            }
//
//        }


        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        double avgCellCycleDurationTumour = 24;
        double avgCellCycleDurationStroma = 32;
        unsigned cancerLabelColour = 3;
//        unsigned macrophageLabelColour = 7;
        MAKE_PTR_ARGS(CellLabel, p_stroma_label, (1));
        MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));
        for (unsigned i=0; i<numStromaNodes; i++)
        {
            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationStroma*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationStroma*1.25);
            p_model->SetHypoxicConcentration(0.2);
            p_model->SetNecroticConcentration(0.05);
            p_model->SetQuiescenceVolumeProportion(0.75);

//            ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
//            p_model->SetDimension(DIM);
//            p_model->SetQuiescentVolumeFraction(0.7);
//            p_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
//            p_cell->GetCellData()->SetItem("csf1", 0);
//            p_cell->GetCellData()->SetItem("egf", 0);
//            p_cell->GetCellData()->SetItem("cxcl12", 0);
//            p_cell->GetCellData()->SetItem("tgf", 0);
//            p_cell->GetCellData()->SetItem("phenotype", -2);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_stroma_label);
            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationStroma*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);

        }

        for (unsigned i=numStromaNodes; i<numStromaNodes+numTumourNodes; i++)
        {

            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationTumour*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationTumour*1.25);
            p_model->SetHypoxicConcentration(0.05);
            p_model->SetNecroticConcentration(0.05);
            // Tumour cells have much higher tolerance to contact inhibition
            p_model->SetQuiescenceVolumeProportion(0.6);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
//            p_cell->GetCellData()->SetItem("csf1", 0);
//            p_cell->GetCellData()->SetItem("egf", 0);
//            p_cell->GetCellData()->SetItem("cxcl12", 0);
//            p_cell->GetCellData()->SetItem("tgf", 0);
//            p_cell->GetCellData()->SetItem("phenotype", -1);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_cancer_label);
//            MARK;
            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationTumour*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);
        }

//        // Add macrophages
//        double numMacrophagesToAdd = 10;
//        double normalize;
//        std::vector<double> randomCoords = {0,0};
//
//        MAKE_PTR_ARGS(CellLabel, p_macrophage_label, (macrophageLabelColour));
//        for (int j = 0; j < numMacrophagesToAdd; j++) {
//
//            // Make Macrophage
//            Arwert2018_MacrophageCellCycle *p_model = new Arwert2018_MacrophageCellCycle;
//            p_model->SetDimension(DIM);
//            p_model->SetPhagocytosisCooldownPeriod(24.0);
//            CellPtr pNewCell(new Cell(p_state, p_model));
//            pNewCell->SetCellProliferativeType(p_diff_type);
//            pNewCell->GetCellData()->SetItem("oxygen", 1);
//            pNewCell->GetCellData()->SetItem("csf1", 0);
//            pNewCell->GetCellData()->SetItem("egf", 0);
//            pNewCell->GetCellData()->SetItem("cxcl12", 0);
//            pNewCell->GetCellData()->SetItem("tgf", 0);
//            pNewCell->GetCellData()->SetItem("phenotype", 0);
//            pNewCell->AddCellProperty(p_macrophage_label);
//
//            Node<DIM> *p_new_node;
//            if (DIM == 2) {
//                p_new_node = new Node<DIM>(nodeNum, false, 0.5 , 2*j + 5);
//                nodes.push_back(p_new_node);
//                nodeNum++;
//            }
//            cells.push_back(pNewCell);
//        }


        NodesOnlyMesh<DIM> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);


        // Make cell population (2D)
        NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

        // Write summary statistic files
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();
        cell_population.AddCellWriter<CellCycleClockWriter>();


        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<DIM> lower(0, 0, 0);
        ChastePoint<DIM> upper(width,height, 0);
        MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));


        //
        // OXYGEN
        //
        double consumptionRate = -0.03;//was -0.03
        double diffusionCoefficient = 1;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_oxygen_pde, (cell_population, consumptionRate,diffusionCoefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_oxygen_bc, (1));
        bool is_neumann_bc = false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        int updateIntervalForPdeInTimesteps = 120/2;
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs<DIM>, p_oxygen_pde_modifier, (p_oxygen_pde, p_oxygen_bc, is_neumann_bc, p_cuboid, 1.0));
        p_oxygen_pde_modifier->SetTimestepInterval(1);
        p_oxygen_pde_modifier->SetDependentVariableName("oxygen");
        p_oxygen_pde_modifier->SetBcsOnBoxBoundary(false); //was false

//
//        //
//        // CXCL12
//        //
//        double cxcl12_decay_rate = -0.02;
//        MAKE_PTR_ARGS(UniformSourceEllipticPde<DIM>, p_cxcl12_pde, (cxcl12_decay_rate));
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_cxcl12_bc, (1));
//        bool is_neumann_bc_cxcl12= false; // Dirichlet BCs
//        // Create a PDE modifier and set the name of the dependent variable in the PDE
//        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs<DIM>, p_cxcl12_pde_modifier, (p_cxcl12_pde, p_cxcl12_bc, is_neumann_bc_cxcl12, p_cuboid, 1.0));
//        p_cxcl12_pde_modifier->SetTimestepInterval(1);
//        p_cxcl12_pde_modifier->SetDependentVariableName("cxcl12");
//        p_cxcl12_pde_modifier->SetBcsOnBoxBoundary(false); //was false
//        p_cxcl12_pde_modifier->SetOutputGradient(true);
//
//        //
//        // CSF-1
//        //
//        double dudtCoefficient = 1.0;
//        double csfDiffusionCoefficient = 0.4;
//        double sourceCoefficient = 0.2;
//        double decayRate_csf = 0.1;
//        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_csf1_pde, (cell_population, dudtCoefficient,csfDiffusionCoefficient,sourceCoefficient));
//        p_csf1_pde->SetCellTypesAsSource({cancerLabelColour});
//        p_csf1_pde->SetDecayRate(decayRate_csf);
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_csf1_bc, (0));
//        bool is_neumann_bc_csf1 = true;
//        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_csf1_pde_modifier, (p_csf1_pde, p_csf1_bc, is_neumann_bc_csf1, p_cuboid, 1.0));
//        p_csf1_pde_modifier->SetDependentVariableName("csf1");
//        p_csf1_pde_modifier->SetBcsOnBoxBoundary(true);
//        p_csf1_pde_modifier->SetOutputGradient(true);
//
//        //
//        // TGF
//        //
//        double dudtCoefficient_tgf = 1.0;
//        double diffusionCoefficient_tgf = 0.1;
//        double sourceCoefficient_tgf = 0.1;
//        double decayRate_tgf = 0.1;
//        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_tgf_pde, (cell_population, dudtCoefficient_tgf,diffusionCoefficient_tgf,sourceCoefficient_tgf));
//        p_tgf_pde->SetCellTypesAsSource({cancerLabelColour});
//        p_tgf_pde->SetDecayRate(decayRate_tgf);
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_tgf_bc, (0));
//        bool is_neumann_bc_tgf = true;
//        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_tgf_pde_modifier, (p_tgf_pde, p_tgf_bc, is_neumann_bc_tgf, p_cuboid, 1.0));
//        p_tgf_pde_modifier->SetDependentVariableName("tgf");
//        p_tgf_pde_modifier->SetBcsOnBoxBoundary(true);
//        p_tgf_pde_modifier->SetOutputGradient(true);
//
//        //
//        // EGF
//        //
//        double dudtCoefficient_egf = 1.0;
//        double diffusionCoefficient_egf = 0.2;
//        double sourceCoefficient_egf = 0.2;
//        double decayRate_egf = 0.1;
//        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_egf_pde, (cell_population, dudtCoefficient_egf,diffusionCoefficient_egf,sourceCoefficient_egf));
//        p_egf_pde->SetCellTypesAsSource({macrophageLabelColour});
//        p_egf_pde->SetDecayRate(decayRate_egf);
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_egf_bc, (0));
//        bool is_neumann_bc_egf = true;
//        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_egf_pde_modifier, (p_egf_pde, p_egf_bc, is_neumann_bc_egf, p_cuboid, 0.5));
//        p_egf_pde_modifier->SetDependentVariableName("egf");
//        p_egf_pde_modifier->SetBcsOnBoxBoundary(true);
//        p_egf_pde_modifier->SetOutputGradient(true);



        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<DIM> simulator(cell_population);
        simulator.AddSimulationModifier(p_oxygen_pde_modifier);
//        simulator.AddSimulationModifier(p_csf1_pde_modifier);
//        simulator.AddSimulationModifier(p_egf_pde_modifier);
//        simulator.AddSimulationModifier(p_tgf_pde_modifier);
//        simulator.AddSimulationModifier(p_cxcl12_pde_modifier);
        std::stringstream output_directory;
        output_directory << "Arwert2018/II_StromaAndTumour_3/";
        simulator.SetOutputDirectory(output_directory.str());
        simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
        simulator.SetEndTime(1000);

//        MAKE_PTR(PhenotypeBasedMacrophagePhagocytosisModifier<DIM>, p_phagoytosis_modifier);
//        p_phagoytosis_modifier->SetCellTypesToEat({cancerLabelColour});
//        simulator.AddSimulationModifier(p_phagoytosis_modifier);


        // Add Brownian motion for all cells
//		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(0.01);
//		simulator.AddForce(p_diffusion_force);

        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

//        // Macrophages go up csf1 gradient
//        // todo make this dependent on phenotype close to 0?
//        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCSF1);
//        p_chemotaxis_MacrophageToCSF1->SetNutrient("csf1");
//        p_chemotaxis_MacrophageToCSF1->SetCellTypesForChemotaxis({macrophageLabelColour});
//        p_chemotaxis_MacrophageToCSF1->SetSensitivity(3);
//        simulator.AddForce(p_chemotaxis_MacrophageToCSF1);
//
//        // Tumour cells go up egf gradient
//        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_TumourCellToEGF);
//        p_chemotaxis_TumourCellToEGF->SetNutrient("egf");
//        p_chemotaxis_TumourCellToEGF->SetCellTypesForChemotaxis({cancerLabelColour});
//        p_chemotaxis_TumourCellToEGF->SetSensitivity(5);
//        simulator.AddForce(p_chemotaxis_TumourCellToEGF);
//
//        // Macrophages go up cxcl12 gradient, depending on phenotype (close to 1)
//        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCXCL12);
//        p_chemotaxis_MacrophageToCXCL12->SetNutrient("cxcl12");
//        p_chemotaxis_MacrophageToCXCL12->SetCellTypesForChemotaxis({macrophageLabelColour});
//        p_chemotaxis_MacrophageToCXCL12->SetSensitivity(5);
//        simulator.AddForce(p_chemotaxis_MacrophageToCXCL12);


        MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(5.0);
        p_force->SetCutOffLength(1.5);
//        p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
//        p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        simulator.AddForce(p_force);
        MARK;

        // Killer which removes apoptotic cells
        MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
        simulator.AddCellKiller(p_apoptosis_killer);

        // Boundary conditions: Add a wall along the x axis
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        // Boundary conditions: Add a wall along the bottom
        c_vector<double,2> normal2 = zero_vector<double>(2);
        c_vector<double,2> point2 = zero_vector<double>(2);
        normal2(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        // Boundary conditions: Add a wall along the top
        c_vector<double,2> normal3 = zero_vector<double>(2);
        c_vector<double,2> point3 = zero_vector<double>(2);
        normal3(1) = 1.0;
        point3(1) = height;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point3, normal3));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        // Boundary conditions: Add a wall along the right hand side
        c_vector<double,2> point4 = zero_vector<double>(2);
        c_vector<double,2> normal4 = zero_vector<double>(2);
        point4(0) = width;
        normal4(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point4, normal4));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);


        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
//		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

    }

    void dontTestSaveSpeedyMacrophages() throw(Exception)
    {
        /*
         * Replicate the system described in Arwert et al (2018): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5946803/pdf/main.pdf
         *
         * Macrophages "extravasate" with phenotype 0 due to CSF1
         * They're recruited up the CSF1 gradient into the tumour
         * Tumour expression of TGF-b causes their phenotype to change towards 1
         * (Interpreted as higher expression of CXCR4)
         * Macrophages with phenotype close to 1 are more sensitive to CXCL12 expressed by perivascular fibroblasts
         * Hence macrophages move back towards the vessel while also expressing EGF (hence attracting tumour cells)
         *
         * Change is unidirectional - i.e., phenotype can increase towards 1, but not away from it
         */

        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        int visualisationOutputFrequencyPerHour = 2;
        const int DIM = 2;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<DIM>*> nodes;
        // Add some nodes

        // Domain size
        double height = 30;
        double width = 15;

        // Initial cell fill size
        double initialCellsHeight = 30;
        double initialStromalCellsWidth = 15;
//        double initialTumourCellsWidth = 5;
        // Stroma Nodes
        unsigned nodeNum=0;
        unsigned numStromaNodes = 0;
        // For staggering rows
        unsigned rownum = 1;
        double offset = 0;
        for (double x=0; x<=initialStromalCellsWidth; x=x+0.75)//x++)
        {
            offset = 0;
            if (rownum % 2)
            {
                offset = 0.75/2;
            }
            for (double y=0+offset; y<=initialCellsHeight; y=y+0.75)//y++)
            {
                if (DIM == 2) {
                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                    nodeNum++;
                    numStromaNodes++;
                }
            }
            rownum++;
        }
        // Tumour Nodes
        unsigned numTumourNodes = 0;
        // Seed one tumour cell
        nodes.push_back(new Node<DIM>(nodeNum, false, initialStromalCellsWidth-2.5, height*0.5));
        nodeNum++;
        numTumourNodes++;
//
//        for (double x=initialStromalCellsWidth+1; x<=initialStromalCellsWidth+initialTumourCellsWidth; x++)
//        {
//            for (double y=0; y<=initialCellsHeight; y++) {
//                if (DIM == 2) {
//                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
//                    nodeNum++;
//                    numTumourNodes++;
//                }
//            }
//
//        }


        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        double avgCellCycleDurationTumour = 24;
        double avgCellCycleDurationStroma = 32;
        unsigned cancerLabelColour = 3;
        unsigned macrophageLabelColour = 7;
        MAKE_PTR_ARGS(CellLabel, p_stroma_label, (1));
        MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));
        for (unsigned i=0; i<numStromaNodes; i++)
        {
            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationStroma*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationStroma*1.25);
            p_model->SetHypoxicConcentration(0.2);
            p_model->SetNecroticConcentration(0.05);
            p_model->SetQuiescenceVolumeProportion(0.75);

//            ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
//            p_model->SetDimension(DIM);
//            p_model->SetQuiescentVolumeFraction(0.7);
//            p_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -2);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_stroma_label);
            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationStroma*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);

        }

        for (unsigned i=numStromaNodes; i<numStromaNodes+numTumourNodes; i++)
        {

            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationTumour*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationTumour*1.25);
            p_model->SetHypoxicConcentration(0.05);
            p_model->SetNecroticConcentration(0.05);
            // Tumour cells have much higher tolerance to contact inhibition
            p_model->SetQuiescenceVolumeProportion(0.6);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -1);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_cancer_label);
//            MARK;
            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationTumour*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);
        }

        // Add macrophages
        double numMacrophagesToAdd = 0;
        double normalize;
        std::vector<double> randomCoords = {0,0};

        MAKE_PTR_ARGS(CellLabel, p_macrophage_label, (macrophageLabelColour));
        for (int j = 0; j < numMacrophagesToAdd; j++) {

            // Make Macrophage
            Arwert2018_MacrophageCellCycle *p_model = new Arwert2018_MacrophageCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetPhagocytosisCooldownPeriod(24.0);
            CellPtr pNewCell(new Cell(p_state, p_model));
            pNewCell->SetCellProliferativeType(p_diff_type);
            pNewCell->GetCellData()->SetItem("oxygen", 1);
            pNewCell->GetCellData()->SetItem("csf1", 0);
            pNewCell->GetCellData()->SetItem("egf", 0);
            pNewCell->GetCellData()->SetItem("cxcl12", 0);
            pNewCell->GetCellData()->SetItem("tgf", 0);
            pNewCell->GetCellData()->SetItem("phenotype", 0);
            pNewCell->AddCellProperty(p_macrophage_label);

            Node<DIM> *p_new_node;
            if (DIM == 2) {
                p_new_node = new Node<DIM>(nodeNum, false, 0.5 , 2*j + 5);
                nodes.push_back(p_new_node);
                nodeNum++;
            }
            cells.push_back(pNewCell);
        }


        NodesOnlyMesh<DIM> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);


        // Make cell population (2D)
        NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

        // Write summary statistic files
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();
        cell_population.AddCellWriter<CellCycleClockWriter>();


        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<DIM> lower(0, 0, 0);
        ChastePoint<DIM> upper(width,height, 0);
        MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));


        //
        // OXYGEN
        //
        double consumptionRate = -0.03;//was -0.03
        double diffusionCoefficient = 1;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_oxygen_pde, (cell_population, consumptionRate,diffusionCoefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_oxygen_bc, (1));
        bool is_neumann_bc = false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        int updateIntervalForPdeInTimesteps = 120/2;
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs<DIM>, p_oxygen_pde_modifier, (p_oxygen_pde, p_oxygen_bc, is_neumann_bc, p_cuboid, 1.0));
        p_oxygen_pde_modifier->SetTimestepInterval(1);
        p_oxygen_pde_modifier->SetDependentVariableName("oxygen");
        p_oxygen_pde_modifier->SetBcsOnBoxBoundary(false); //was false


        //
        // CXCL12
        //
        double cxcl12_decay_rate = -0.02;
        MAKE_PTR_ARGS(UniformSourceEllipticPde<DIM>, p_cxcl12_pde, (cxcl12_decay_rate));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_cxcl12_bc, (1));
        bool is_neumann_bc_cxcl12= false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs<DIM>, p_cxcl12_pde_modifier, (p_cxcl12_pde, p_cxcl12_bc, is_neumann_bc_cxcl12, p_cuboid, 1.0));
        p_cxcl12_pde_modifier->SetTimestepInterval(1);
        p_cxcl12_pde_modifier->SetDependentVariableName("cxcl12");
        p_cxcl12_pde_modifier->SetBcsOnBoxBoundary(false); //was false
        p_cxcl12_pde_modifier->SetOutputGradient(true);

        //
        // CSF-1
        //
        double dudtCoefficient = 1.0;
        double csfDiffusionCoefficient = 1.0;
        double sourceCoefficient = 0.3;
        double decayRate_csf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_csf1_pde, (cell_population, dudtCoefficient,csfDiffusionCoefficient,sourceCoefficient));
        p_csf1_pde->SetCellTypesAsSource({cancerLabelColour});
        p_csf1_pde->SetDecayRate(decayRate_csf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_csf1_bc, (0));
        bool is_neumann_bc_csf1 = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_csf1_pde_modifier, (p_csf1_pde, p_csf1_bc, is_neumann_bc_csf1, p_cuboid, 1.0));
        p_csf1_pde_modifier->SetDependentVariableName("csf1");
        p_csf1_pde_modifier->SetBcsOnBoxBoundary(true);
        p_csf1_pde_modifier->SetOutputGradient(true);

        //
        // TGF
        //
        double dudtCoefficient_tgf = 1.0;
        double diffusionCoefficient_tgf = 0.1;
        double sourceCoefficient_tgf = 0.1;
        double decayRate_tgf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_tgf_pde, (cell_population, dudtCoefficient_tgf,diffusionCoefficient_tgf,sourceCoefficient_tgf));
        p_tgf_pde->SetCellTypesAsSource({cancerLabelColour});
        p_tgf_pde->SetDecayRate(decayRate_tgf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_tgf_bc, (0));
        bool is_neumann_bc_tgf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_tgf_pde_modifier, (p_tgf_pde, p_tgf_bc, is_neumann_bc_tgf, p_cuboid, 1.0));
        p_tgf_pde_modifier->SetDependentVariableName("tgf");
        p_tgf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_tgf_pde_modifier->SetOutputGradient(true);

        //
        // EGF
        //
        double dudtCoefficient_egf = 1.0;
        double diffusionCoefficient_egf = 0.2;
        double sourceCoefficient_egf = 1.0;
        double decayRate_egf = 1.0;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_egf_pde, (cell_population, dudtCoefficient_egf,diffusionCoefficient_egf,sourceCoefficient_egf));
        p_egf_pde->SetCellTypesAsSource({macrophageLabelColour});
        p_egf_pde->SetDecayRate(decayRate_egf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_egf_bc, (0));
        bool is_neumann_bc_egf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_egf_pde_modifier, (p_egf_pde, p_egf_bc, is_neumann_bc_egf, p_cuboid, 0.5));
        p_egf_pde_modifier->SetDependentVariableName("egf");
        p_egf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_egf_pde_modifier->SetOutputGradient(true);



        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<DIM> simulator(cell_population);
        simulator.AddSimulationModifier(p_oxygen_pde_modifier);
        simulator.AddSimulationModifier(p_csf1_pde_modifier);
        simulator.AddSimulationModifier(p_egf_pde_modifier);
        simulator.AddSimulationModifier(p_tgf_pde_modifier);
        simulator.AddSimulationModifier(p_cxcl12_pde_modifier);
        std::stringstream output_directory;
        output_directory << "Arwert2018/III_SpeedyMacrophages/";
        simulator.SetOutputDirectory(output_directory.str());
        simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
        simulator.SetEndTime(300);

        MAKE_PTR(PhenotypeBasedMacrophagePhagocytosisModifier<DIM>, p_phagoytosis_modifier);
        p_phagoytosis_modifier->SetCellTypesToEat({cancerLabelColour});
        simulator.AddSimulationModifier(p_phagoytosis_modifier);


////         Add Brownian motion for all cells
//		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(0.01);
//		simulator.AddForce(p_diffusion_force);

        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Macrophages go up csf1 gradient
        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCSF1);
        p_chemotaxis_MacrophageToCSF1->SetNutrient("csf1");
        p_chemotaxis_MacrophageToCSF1->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCSF1->SetSensitivity(10);
        simulator.AddForce(p_chemotaxis_MacrophageToCSF1);

        // Tumour cells go up egf gradient
        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_TumourCellToEGF);
        p_chemotaxis_TumourCellToEGF->SetNutrient("egf");
        p_chemotaxis_TumourCellToEGF->SetCellTypesForChemotaxis({cancerLabelColour});
        p_chemotaxis_TumourCellToEGF->SetSensitivity(10);
        simulator.AddForce(p_chemotaxis_TumourCellToEGF);

        // Macrophages go up cxcl12 gradient, depending on phenotype (close to 1)
        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCXCL12);
        p_chemotaxis_MacrophageToCXCL12->SetNutrient("cxcl12");
        p_chemotaxis_MacrophageToCXCL12->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCXCL12->SetSensitivity(10);
        simulator.AddForce(p_chemotaxis_MacrophageToCXCL12);


        MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(5.0);
        p_force->SetCutOffLength(1.5);
//        p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
//        p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        simulator.AddForce(p_force);
        MARK;

        // Killer which removes apoptotic cells
        MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
        simulator.AddCellKiller(p_apoptosis_killer);

        // Boundary conditions: Add a wall along the x axis
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        // Boundary conditions: Add a wall along the bottom
        c_vector<double,2> normal2 = zero_vector<double>(2);
        c_vector<double,2> point2 = zero_vector<double>(2);
        normal2(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        // Boundary conditions: Add a wall along the top
        c_vector<double,2> normal3 = zero_vector<double>(2);
        c_vector<double,2> point3 = zero_vector<double>(2);
        normal3(1) = 1.0;
        point3(1) = height;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point3, normal3));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        // Boundary conditions: Add a wall along the right hand side
        c_vector<double,2> point4 = zero_vector<double>(2);
        c_vector<double,2> normal4 = zero_vector<double>(2);
        point4(0) = width;
        normal4(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point4, normal4));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);


        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

    }

    void dontTestLoadSpeedyMacrophages() throw(Exception)
    {
        // Load in the saved simulation
        const int DIM = 2;
        double height = 30;
        double width = 15;

        double startTime = 300;
        double newDuration = 48;
        double visualisationOutputFrequencyPerHour = 60;

        MARK;
        OffLatticeSimulation<DIM>* p_simulator = CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Load("Arwert2018/III_SpeedyMacrophages", startTime);
        MARK;
        p_simulator->SetEndTime(startTime + newDuration);
        MARK;
        p_simulator->SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
        MARK;

        // We remake the cuboid for the PDE FE Mesh
        ChastePoint<DIM> lower(0,0,0);
        ChastePoint<DIM> upper(width,height,0);
        MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));
        MARK;
        // Now loop over simulation modifiers
        for (typename std::vector<boost::shared_ptr<AbstractCellBasedSimulationModifier<DIM> > >::iterator iter = p_simulator->GetSimulationModifiers()->begin();  iter != p_simulator->GetSimulationModifiers()->end(); ++iter)
        {
            MARK;
            PRINT_VARIABLE((*iter)->GetIdentifier());
            // Attempt to convert modifier to EllipticBoxDomainPdeModifier. If modifier is of a different type, this should return NULL
            boost::shared_ptr<EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs<DIM> > p_newModifier = boost::dynamic_pointer_cast<EllipticBoxDomainPdeModifier_VariableTimestep_MixedBCs<DIM> >(*iter);

            // If not NULL, we have found the correct modifier
            if(p_newModifier)
            {
                // Now regenerate the FeMesh.
                p_newModifier->GenerateFeMesh(p_cuboid,1.0);
                p_newModifier->SetTimestepInterval(1);
            }
            else{
                // OK, now do it again in case it's parabolic

                // Attempt to convert modifier to ParabolicBoxDomainPdeModifier. If modifier is of a different type, this should return NULL
                boost::shared_ptr<ParabolicBoxDomainPdeModifier<DIM> > p_newModifier = boost::dynamic_pointer_cast<ParabolicBoxDomainPdeModifier<DIM> >(*iter);
                if(p_newModifier)
                {
                    // Now regenerate the FeMesh.
                    p_newModifier->GenerateFeMesh(p_cuboid,1.0);
                }
            }

        }
        MARK;

//         Add Brownian motion for cancer cells
        unsigned stromaLabelColour = 1;
        unsigned cancerLabelColour = 3;
        unsigned macrophageLabelColour = 7;
//		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(0.01);
//        p_diffusion_force->SetCellTypesForDiffusion({stromaLabelColour,cancerLabelColour,6});
//		p_simulator->AddForce(p_diffusion_force);

        MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_forceMac);
        p_diffusion_forceMac->SetDiffusionScalingConstant(1);
        p_diffusion_forceMac->SetCellTypesForDiffusion({macrophageLabelColour});
        p_simulator->AddForce(p_diffusion_forceMac);

        p_simulator->rGetCellPopulation().AddCellWriter<CellCycleClockWriter>();

        // Add macrophages
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        double numMacrophagesToAdd = 10;
        double normalize;
        std::vector<double> randomCoords = {0,0};
        MAKE_PTR_ARGS(CellLabel, p_macrophage_label, (macrophageLabelColour));
        for (int j = 0; j < numMacrophagesToAdd; j++) {

            // Make Macrophage
            Arwert2018_MacrophageCellCycle *p_model = new Arwert2018_MacrophageCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetPhagocytosisCooldownPeriod(1.0);
            CellPtr pNewCell(new Cell(p_state, p_model));
            pNewCell->SetCellProliferativeType(p_diff_type);
            pNewCell->GetCellData()->SetItem("oxygen", 1);
            pNewCell->GetCellData()->SetItem("csf1", 0);
            pNewCell->GetCellData()->SetItem("egf", 0);
            pNewCell->GetCellData()->SetItem("cxcl12", 0);
            pNewCell->GetCellData()->SetItem("tgf", 0);
            pNewCell->GetCellData()->SetItem("phenotype", 0);
            pNewCell->AddCellProperty(p_macrophage_label);

            Node<DIM> *p_new_node;
            if (DIM == 2) {
                // 100000 is a placeholder that gets overridden when the node is added to the mesh
                p_new_node = new Node<DIM>(100000, false, 0.5 , 2*j + 5);
                p_new_node->SetRadius(0.5);

                auto nbcp = dynamic_cast<NodeBasedCellPopulation<DIM>* >(&p_simulator->rGetCellPopulation());
                unsigned new_node_index = nbcp->rGetMesh().AddNode(p_new_node);

                // Update cells vector
                p_simulator->rGetCellPopulation().rGetCells().push_back(pNewCell);

                // Update mappings between cells and location indices
                p_simulator->rGetCellPopulation().SetCellUsingLocationIndex(new_node_index, pNewCell);
            }
        }

        p_simulator->Solve();
        MARK;
        CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(p_simulator);
        MARK;
        delete p_simulator;

    }

    void dontTestMacrophageInfiltrationIntoSpheroid() throw(Exception)
    {
        /*
         * Replicate the system described in Arwert et al (2018): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5946803/pdf/main.pdf
         *
         * Macrophages "extravasate" with phenotype 0 due to CSF1
         * They're recruited up the CSF1 gradient into the tumour
         * Tumour expression of TGF-b causes their phenotype to change towards 1
         * (Interpreted as higher expression of CXCR4)
         * Macrophages with phenotype close to 1 are more sensitive to CXCL12 expressed by perivascular fibroblasts
         * Hence macrophages move back towards the vessel while also expressing EGF (hence attracting tumour cells)
         *
         * Change is unidirectional - i.e., phenotype can increase towards 1, but not away from it
         */

        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        int visualisationOutputFrequencyPerHour = 2;
        const int DIM = 2;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<DIM>*> nodes;
        // Add some nodes
        double cubeDomainDistanceToBoundary = 30;

        // Tumour Nodes
        unsigned nodeNum=0;
        double initialRadius = 5.0;
        for (double x=-initialRadius; x<initialRadius+1; x++)
        {
            for (double y=-initialRadius; y<initialRadius+1; y++) {
                    if (pow(x, 2) + pow(y, 2) < pow(initialRadius, 2)) {
                        nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                        nodeNum++;
                    }
            }

        }



        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        double avgCellCycleDurationTumour = 24;
        unsigned cancerLabelColour = 3;
        unsigned macrophageLabelColour = 7;
        unsigned complexLabelColour = 8;
        MAKE_PTR_ARGS(CellLabel, p_macrophage_label, (macrophageLabelColour));
        MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));

        for (unsigned i=0; i<nodeNum; i++)
        {

            SpheroidComparisonCellCycle* p_model = new SpheroidComparisonCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationTumour*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationTumour*1.25);
            p_model->SetHypoxicConcentration(0.2);
            p_model->SetNecroticConcentration(0.1);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -1);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_cancer_label);

            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationTumour*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);
        }


        NodesOnlyMesh<DIM> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);


        // Make cell population (2D)
        NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

        // Write summary statistic files
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();
        cell_population.AddCellWriter<CellCycleClockWriter>();


        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        c_vector<double,DIM> centroid = cell_population.GetCentroidOfCellPopulation();
        ChastePoint<DIM> lower(-cubeDomainDistanceToBoundary, -cubeDomainDistanceToBoundary,-cubeDomainDistanceToBoundary);
        ChastePoint<DIM> upper(cubeDomainDistanceToBoundary, cubeDomainDistanceToBoundary,cubeDomainDistanceToBoundary);
        MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));



        //
        // OXYGEN
        //
        double consumptionRate = -0.03;//was -0.03
        double diffusionCoefficient = 1;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_oxygen_pde, (cell_population, consumptionRate,diffusionCoefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_oxygen_bc, (1));
        bool is_neumann_bc = false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
//        int updateIntervalForPdeInTimesteps = 120/2;
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_oxygen_pde_modifier, (p_oxygen_pde, p_oxygen_bc, is_neumann_bc, p_cuboid, 1.0));
        p_oxygen_pde_modifier->SetTimestepInterval(60);
        p_oxygen_pde_modifier->SetDependentVariableName("oxygen");
        p_oxygen_pde_modifier->SetBcsOnBoxBoundary(false); //was false


        //
        // CXCL12
        //
        double cxcl12_decay_rate = -0.02;
        MAKE_PTR_ARGS(UniformSourceEllipticPde<DIM>, p_cxcl12_pde, (cxcl12_decay_rate));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_cxcl12_bc, (1));
        bool is_neumann_bc_cxcl12= false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifierVariableTimestep<DIM>, p_cxcl12_pde_modifier, (p_cxcl12_pde, p_cxcl12_bc, is_neumann_bc_cxcl12, p_cuboid, 1.0));
        p_cxcl12_pde_modifier->SetTimestepInterval(60);
        p_cxcl12_pde_modifier->SetDependentVariableName("cxcl12");
        p_cxcl12_pde_modifier->SetBcsOnBoxBoundary(false); //was false
        p_cxcl12_pde_modifier->SetOutputGradient(true);


        //
        // CSF-1
        //
        double dudtCoefficient = 1.0;
        double csfDiffusionCoefficient = 0.4;
        double sourceCoefficient = 0.2;
        double decayRate_csf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_csf1_pde, (cell_population, dudtCoefficient,csfDiffusionCoefficient,sourceCoefficient));
        p_csf1_pde->SetCellTypesAsSource({cancerLabelColour});
        p_csf1_pde->SetDecayRate(decayRate_csf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_csf1_bc, (0));
        bool is_neumann_bc_csf1 = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_csf1_pde_modifier, (p_csf1_pde, p_csf1_bc, is_neumann_bc_csf1, p_cuboid, 1.0));
        p_csf1_pde_modifier->SetDependentVariableName("csf1");
        p_csf1_pde_modifier->SetBcsOnBoxBoundary(true);
        p_csf1_pde_modifier->SetOutputGradient(true);

        //
        // TGF
        //
        double dudtCoefficient_tgf = 1.0;
        double diffusionCoefficient_tgf = 0.1;
        double sourceCoefficient_tgf = 0.1;
        double decayRate_tgf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_tgf_pde, (cell_population, dudtCoefficient_tgf,diffusionCoefficient_tgf,sourceCoefficient_tgf));
        p_tgf_pde->SetCellTypesAsSource({cancerLabelColour});
        p_tgf_pde->SetDecayRate(decayRate_tgf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_tgf_bc, (0));
        bool is_neumann_bc_tgf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_tgf_pde_modifier, (p_tgf_pde, p_tgf_bc, is_neumann_bc_tgf, p_cuboid, 1.0));
        p_tgf_pde_modifier->SetDependentVariableName("tgf");
        p_tgf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_tgf_pde_modifier->SetOutputGradient(true);

        //
        // EGF
        //
//        double dudtCoefficient_egf = 1.0;
//        double diffusionCoefficient_egf = 0.2;
//        double sourceCoefficient_egf = 0.2;
//        double decayRate_egf = 0.1;
//        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_egf_pde, (cell_population, dudtCoefficient_egf,diffusionCoefficient_egf,sourceCoefficient_egf));
//        p_egf_pde->SetCellTypesAsSource({macrophageLabelColour});
//        p_egf_pde->SetDecayRate(decayRate_egf);
//        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_egf_bc, (0));
//        bool is_neumann_bc_egf = true;
//        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_egf_pde_modifier, (p_egf_pde, p_egf_bc, is_neumann_bc_egf, p_cuboid, 1.0));
//        p_egf_pde_modifier->SetDependentVariableName("egf");
//        p_egf_pde_modifier->SetBcsOnBoxBoundary(true);
//        p_egf_pde_modifier->SetOutputGradient(true);





        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<DIM> simulator(cell_population);
        std::stringstream output_directory;
        output_directory << "Arwert2018/IV_Spheroid_reordered/";
        simulator.SetOutputDirectory(output_directory.str());
        simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
        double growthTime = 150;
        double infiltrationDuration = 150;
        simulator.SetEndTime(growthTime+infiltrationDuration);

        simulator.AddSimulationModifier(p_oxygen_pde_modifier);
        simulator.AddSimulationModifier(p_cxcl12_pde_modifier);
        simulator.AddSimulationModifier(p_csf1_pde_modifier);
        simulator.AddSimulationModifier(p_tgf_pde_modifier);
//        simulator.AddSimulationModifier(p_egf_pde_modifier);



        MAKE_PTR(Arwert2018_AddMacrophagesModifier<DIM>, p_addMacs_modifier);
        p_addMacs_modifier->SetTimeToAddMacrophages(growthTime);
        p_addMacs_modifier->SetNumberOfMacrophagesToAdd(100);
        p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(infiltrationDuration);
        simulator.AddSimulationModifier(p_addMacs_modifier);

        MAKE_PTR(EgfBasedMacrophageTumourCellBindingModifier<DIM>, p_binding_modifier);
        p_binding_modifier->SetMacrophageLabel(macrophageLabelColour);
        p_binding_modifier->SetComplexLabel(complexLabelColour);
        p_binding_modifier->SetCellTypesToBindWith({cancerLabelColour});
//        simulator.AddSimulationModifier(p_binding_modifier);

        MAKE_PTR(PhenotypeBasedMacrophagePhagocytosisModifier<DIM>, p_phagoytosis_modifier);
        p_phagoytosis_modifier->SetCellTypesToEat({cancerLabelColour});
//        simulator.AddSimulationModifier(p_phagoytosis_modifier);


        // Add Brownian motion for all cells
//		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(0.01);
//		simulator.AddForce(p_diffusion_force);

        // Add Brownian motion for macrophages
        MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_forceMac);
        p_diffusion_forceMac->SetDiffusionScalingConstant(1);
        p_diffusion_forceMac->SetCellTypesForDiffusion({macrophageLabelColour});
        simulator.AddForce(p_diffusion_forceMac);

        // Add Brownian motion for macrophage/cancer cell complexes
        MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_forceComplex);
        p_diffusion_forceMac->SetDiffusionScalingConstant(0.5);
        p_diffusion_forceMac->SetCellTypesForDiffusion({complexLabelColour});
        simulator.AddForce(p_diffusion_forceComplex);


        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Macrophages go up csf1 gradient
        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCSF1);
        p_chemotaxis_MacrophageToCSF1->SetNutrient("csf1");
        p_chemotaxis_MacrophageToCSF1->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCSF1->SetSensitivity(2.5);
        simulator.AddForce(p_chemotaxis_MacrophageToCSF1);

//        // Tumour cells go up egf gradient
//        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_TumourCellToEGF);
//        p_chemotaxis_TumourCellToEGF->SetNutrient("egf");
//        p_chemotaxis_TumourCellToEGF->SetCellTypesForChemotaxis({cancerLabelColour});
//        p_chemotaxis_TumourCellToEGF->SetSensitivity(10);
//        simulator.AddForce(p_chemotaxis_TumourCellToEGF);

        // Macrophages go up cxcl12 gradient, depending on phenotype (close to 1)
        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCXCL12);
        p_chemotaxis_MacrophageToCXCL12->SetNutrient("cxcl12");
        p_chemotaxis_MacrophageToCXCL12->SetCellTypesForChemotaxis({macrophageLabelColour,complexLabelColour});
        p_chemotaxis_MacrophageToCXCL12->SetSensitivity(10);
        simulator.AddForce(p_chemotaxis_MacrophageToCXCL12);


        MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(15.0);
        p_force->SetCutOffLength(1.5);
//        p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
//        p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        simulator.AddForce(p_force);
        MARK;

        // Killer which removes apoptotic cells
        MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
        simulator.AddCellKiller(p_apoptosis_killer);

        MAKE_PTR(ExternalPressureForceOnConcaveHull<DIM>, p_pressure);
        p_pressure->SetPressure(5);
        simulator.AddForce(p_pressure);



        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);



        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }

    }

    void dontTestInVivoGeometryMacrophageInfiltration() throw(Exception)
    {
        /*
         * Replicate the system described in Arwert et al (2018): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5946803/pdf/main.pdf
         *
         * Macrophages "extravasate" with phenotype 0 due to CSF1
         * They're recruited up the CSF1 gradient into the tumour
         * Tumour expression of TGF-b causes their phenotype to change towards 1
         * (Interpreted as higher expression of CXCR4)
         * Macrophages with phenotype close to 1 are more sensitive to CXCL12 expressed by perivascular fibroblasts
         * Hence macrophages move back towards the vessel while also expressing EGF (hence attracting tumour cells)
         *
         * Change is unidirectional - i.e., phenotype can increase towards 1, but not away from it
         */

        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        int visualisationOutputFrequencyPerHour = 2;
        const int DIM = 2;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<DIM>*> nodes;
        // Add some nodes

        // Domain size
        double domainHeight = 25;
        double domainWidth = 25;

        // Fill the domain with stromal cells

        // Stroma Nodes
        unsigned nodeNum=0;
        unsigned numStromaNodes = 0;
        // For staggering rows
        unsigned rownum = 1;
        double offset;
        for (double x=0; x<=domainWidth; x=x+0.75)//x++)
        {
            offset = 0;
            if (rownum % 2)
            {
                offset = 0.75/2;
            }
            for (double y=0+offset; y<=domainHeight; y=y+0.75)//y++)
            {
                if (DIM == 2) {
                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                    nodeNum++;
                    numStromaNodes++;
                }
            }
            rownum++;
        }
        // Tumour Nodes
        unsigned numTumourNodes = 0;
        // Seed a cluster of tumour cells
        // Make sure it's not right at the edge
        double x = domainWidth*0.75;//0.8*domainWidth*RandomNumberGenerator::Instance()->ranf() + 0.1*domainWidth;
        double y = domainHeight*0.74;//0.8*domainHeight*RandomNumberGenerator::Instance()->ranf() + 0.1*domainHeight;
        nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x-0.25, y));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x, y-0.25));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x-0.25, y-0.25));
        nodeNum++;
        numTumourNodes = numTumourNodes+4;


        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        double avgCellCycleDurationTumour = 24;
        double avgCellCycleDurationStroma = 32;

        unsigned cancerLabelColour = 3;
        unsigned macrophageLabelColour = 7;
        MAKE_PTR_ARGS(CellLabel, p_stroma_label, (1));
        MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));
        for (unsigned i=0; i<numStromaNodes; i++)
        {
            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationStroma*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationStroma*1.25);
            p_model->SetHypoxicConcentration(0.2);
            p_model->SetNecroticConcentration(0.05);
            p_model->SetQuiescenceVolumeProportion(0.75);

//            ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
//            p_model->SetDimension(DIM);
//            p_model->SetQuiescentVolumeFraction(0.7);
//            p_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -2);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_stroma_label);
            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationStroma*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);

        }

        for (unsigned i=numStromaNodes; i<numStromaNodes+numTumourNodes; i++)
        {

            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationTumour*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationTumour*1.25);
            p_model->SetHypoxicConcentration(0.05);
            p_model->SetNecroticConcentration(0.05);
            // Tumour cells have much higher tolerance to contact inhibition
            p_model->SetQuiescenceVolumeProportion(0.6);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -1);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_cancer_label);

            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationTumour*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);
        }



        NodesOnlyMesh<DIM> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);


        // Make cell population (2D)
        NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

        // Write summary statistic files
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();
        cell_population.AddCellWriter<CellCycleClockWriter>();


        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<DIM> lower(0, 0, 0);
        ChastePoint<DIM> upper(domainWidth,domainHeight, 0);
        MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));


        //
        // OXYGEN
        //

        // We have Neumann BCs along the boundaries of the simulation domain
        // Dirichlet BCs are set at blood vessel point sources
        // Choose BVs to be random mesh points, rather than picking random points and then meshing
        std::vector<ChastePoint<DIM> > vesselLocations;
        // Small x
        vesselLocations.push_back(ChastePoint<DIM> (6,10));
        vesselLocations.push_back(ChastePoint<DIM> (8,12));
        vesselLocations.push_back(ChastePoint<DIM> (5,16));
        vesselLocations.push_back(ChastePoint<DIM> (7,20));
        vesselLocations.push_back(ChastePoint<DIM> (5,23));
//        vesselLocations.push_back(ChastePoint<DIM> (10,24));
//        vesselLocations.push_back(ChastePoint<DIM> (9,26));
//        vesselLocations.push_back(ChastePoint<DIM> (5,31));
//        vesselLocations.push_back(ChastePoint<DIM> (7,35));

        // Small y
        vesselLocations.push_back(ChastePoint<DIM> (7,4));
        vesselLocations.push_back(ChastePoint<DIM> (12,6));
        vesselLocations.push_back(ChastePoint<DIM> (16,7));
        vesselLocations.push_back(ChastePoint<DIM> (20,3));
        vesselLocations.push_back(ChastePoint<DIM> (23,6));
//        vesselLocations.push_back(ChastePoint<DIM> (28,5));
//        vesselLocations.push_back(ChastePoint<DIM> (32,8));
//        vesselLocations.push_back(ChastePoint<DIM> (36,9));

        double consumptionRate = -0.03;//was -0.03
        double diffusionCoefficient = 1;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_oxygen_pde, (cell_population, consumptionRate,diffusionCoefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_oxygen_bc, (1));
        bool is_neumann_bc = false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        int updateIntervalForPdeInTimesteps = 120/2;
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>, p_oxygen_pde_modifier, (p_oxygen_pde, p_oxygen_bc, is_neumann_bc, p_cuboid, 1.0));
        p_oxygen_pde_modifier->SetTimestepInterval(1);
        p_oxygen_pde_modifier->SetVesselLocations(vesselLocations);
        p_oxygen_pde_modifier->SetDependentVariableName("oxygen");
        p_oxygen_pde_modifier->SetBcsOnBoxBoundary(false); //was false


        //
        // CXCL12
        //
        double cxcl12_decay_rate = -0.02;
        MAKE_PTR_ARGS(UniformSourceEllipticPde<DIM>, p_cxcl12_pde, (cxcl12_decay_rate));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_cxcl12_bc, (1));
        bool is_neumann_bc_cxcl12= false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>, p_cxcl12_pde_modifier, (p_cxcl12_pde, p_cxcl12_bc, is_neumann_bc_cxcl12, p_cuboid, 1.0));
        p_cxcl12_pde_modifier->SetTimestepInterval(1);
        p_cxcl12_pde_modifier->SetDependentVariableName("cxcl12");
        p_cxcl12_pde_modifier->SetBcsOnBoxBoundary(false); //was false
        p_cxcl12_pde_modifier->SetOutputGradient(true);
        p_cxcl12_pde_modifier->SetVesselLocations(vesselLocations);
        //
        // CSF-1
        //
        double dudtCoefficient = 1.0;
        double csfDiffusionCoefficient = 1.0;
        double sourceCoefficient = 0.25;
        double decayRate_csf = 0.02;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_csf1_pde, (cell_population, dudtCoefficient,csfDiffusionCoefficient,sourceCoefficient));
        p_csf1_pde->SetCellTypesAsSource({cancerLabelColour});
        p_csf1_pde->SetDecayRate(decayRate_csf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_csf1_bc, (0));
        bool is_neumann_bc_csf1 = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_csf1_pde_modifier, (p_csf1_pde, p_csf1_bc, is_neumann_bc_csf1, p_cuboid, 1.0));
        p_csf1_pde_modifier->SetDependentVariableName("csf1");
        p_csf1_pde_modifier->SetBcsOnBoxBoundary(true);
        p_csf1_pde_modifier->SetOutputGradient(true);

        //
        // TGF
        //
        double dudtCoefficient_tgf = 1.0;
        double diffusionCoefficient_tgf = 0.1;
        double sourceCoefficient_tgf = 0.1;
        double decayRate_tgf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_tgf_pde, (cell_population, dudtCoefficient_tgf,diffusionCoefficient_tgf,sourceCoefficient_tgf));
        p_tgf_pde->SetCellTypesAsSource({cancerLabelColour});
        p_tgf_pde->SetDecayRate(decayRate_tgf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_tgf_bc, (0));
        bool is_neumann_bc_tgf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_tgf_pde_modifier, (p_tgf_pde, p_tgf_bc, is_neumann_bc_tgf, p_cuboid, 1.0));
        p_tgf_pde_modifier->SetDependentVariableName("tgf");
        p_tgf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_tgf_pde_modifier->SetOutputGradient(true);

        //
        // EGF
        //
        double dudtCoefficient_egf = 1.0;
        double diffusionCoefficient_egf = 0.2;
        double sourceCoefficient_egf = 0.2;
        double decayRate_egf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited_phenotypeDependent<DIM>, p_egf_pde, (cell_population, dudtCoefficient_egf,diffusionCoefficient_egf,sourceCoefficient_egf));
        p_egf_pde->SetCellTypesAsSource({macrophageLabelColour});
        p_egf_pde->SetDecayRate(decayRate_egf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_egf_bc, (0));
        bool is_neumann_bc_egf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_egf_pde_modifier, (p_egf_pde, p_egf_bc, is_neumann_bc_egf, p_cuboid, 0.5));
        p_egf_pde_modifier->SetDependentVariableName("egf");
        p_egf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_egf_pde_modifier->SetOutputGradient(true);



        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<DIM> simulator(cell_population);
        simulator.AddSimulationModifier(p_oxygen_pde_modifier);
        simulator.AddSimulationModifier(p_csf1_pde_modifier);
        simulator.AddSimulationModifier(p_egf_pde_modifier);
        simulator.AddSimulationModifier(p_tgf_pde_modifier);
        simulator.AddSimulationModifier(p_cxcl12_pde_modifier);
        std::stringstream output_directory;
        output_directory << "Arwert2018/V_PointVessels_v/";
        simulator.SetOutputDirectory(output_directory.str());
        simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
        double burnInPeriod = 240;
        double infiltrationDuration = 24*28;
        double endTime = burnInPeriod+infiltrationDuration;
        simulator.SetEndTime(endTime);


        MAKE_PTR(Arwert2018_AddMacrophagesFromPointVesselModifier<DIM>, p_addMacs_modifier);
        p_addMacs_modifier->SetTgfbThresholdForMacrophagePhenotypeSwitch(1);
        p_addMacs_modifier->SetTimeToAddMacrophages(burnInPeriod);
        p_addMacs_modifier->SetNumberOfMacrophagesToAdd(50);
        p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(infiltrationDuration);
        p_addMacs_modifier->SetVesselLocations(vesselLocations);
        simulator.AddSimulationModifier(p_addMacs_modifier);

//        MAKE_PTR(PhenotypeBasedMacrophagePhagocytosisModifier<DIM>, p_phagoytosis_modifier);
//        p_phagoytosis_modifier->SetCellTypesToEat({cancerLabelColour});
//        simulator.AddSimulationModifier(p_phagoytosis_modifier);


        // Add Brownian motion for all cells
//		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(0.01);
//		simulator.AddForce(p_diffusion_force);

        // Add Brownian motion for macrophages
        MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_forceMac);
        p_diffusion_forceMac->SetDiffusionScalingConstant(0.01);
        p_diffusion_forceMac->SetCellTypesForDiffusion({macrophageLabelColour});
        simulator.AddForce(p_diffusion_forceMac);

        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Macrophages go up csf1 gradient
        // todo make this dependent on phenotype close to 0?
        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCSF1);
        p_chemotaxis_MacrophageToCSF1->SetNutrient("csf1");
        p_chemotaxis_MacrophageToCSF1->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCSF1->SetSensitivity(1);
        simulator.AddForce(p_chemotaxis_MacrophageToCSF1);

        // Tumour cells go up egf gradient
        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_TumourCellToEGF);
        p_chemotaxis_TumourCellToEGF->SetNutrient("egf");
        p_chemotaxis_TumourCellToEGF->SetCellTypesForChemotaxis({cancerLabelColour});
        p_chemotaxis_TumourCellToEGF->SetSensitivity(2.5);
        simulator.AddForce(p_chemotaxis_TumourCellToEGF);

        // Macrophages go up cxcl12 gradient, depending on phenotype (close to 1)
        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCXCL12);
        p_chemotaxis_MacrophageToCXCL12->SetNutrient("cxcl12");
        p_chemotaxis_MacrophageToCXCL12->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCXCL12->SetSensitivity(1);
        simulator.AddForce(p_chemotaxis_MacrophageToCXCL12);


        MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(5.0);
        p_force->SetCutOffLength(1.5);
//        p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
//        p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        simulator.AddForce(p_force);
        MARK;

        // Killer which removes apoptotic cells
        MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
        simulator.AddCellKiller(p_apoptosis_killer);

        // Boundary conditions: Add a wall along the x axis
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        // Boundary conditions: Add a wall along the bottom
        c_vector<double,2> normal2 = zero_vector<double>(2);
        c_vector<double,2> point2 = zero_vector<double>(2);
        normal2(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        // Boundary conditions: Add a wall along the top
        c_vector<double,2> normal3 = zero_vector<double>(2);
        c_vector<double,2> point3 = zero_vector<double>(2);
        normal3(1) = 1.0;
        point3(1) = domainHeight;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point3, normal3));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        // Boundary conditions: Add a wall along the right hand side
        c_vector<double,2> point4 = zero_vector<double>(2);
        c_vector<double,2> normal4 = zero_vector<double>(2);
        point4(0) = domainWidth;
        normal4(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point4, normal4));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);



        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
//		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void dontTestAlmostFinalSimulation() throw(Exception)
    {
        /*
         * Replicate the system described in Arwert et al (2018): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5946803/pdf/main.pdf
         *
         * Macrophages "extravasate" with phenotype 0 due to CSF1
         * They're recruited up the CSF1 gradient into the tumour
         * Tumour expression of TGF-b causes their phenotype to change towards 1
         * (Interpreted as higher expression of CXCR4)
         * Macrophages with phenotype close to 1 are more sensitive to CXCL12 expressed by perivascular fibroblasts
         * Hence macrophages move back towards the vessel while also expressing EGF (hence attracting tumour cells)
         *
         * Change is unidirectional - i.e., phenotype can increase towards 1, but not away from it
         */

        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        int visualisationOutputFrequencyPerHour = 2;
        const int DIM = 2;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<DIM>*> nodes;
        // Add some nodes

        // Domain size
        double domainHeight = 40;
        double domainWidth = 40;

        // Fill the domain with stromal cells

        // Stroma Nodes
        unsigned nodeNum=0;
        unsigned numStromaNodes = 0;
        // For staggering rows
        unsigned rownum = 1;
        double offset;
        for (double x=0; x<=domainWidth; x=x+0.75)//x++)
        {
            offset = 0;
            if (rownum % 2)
            {
                offset = 0.75/2;
            }
            for (double y=0+offset; y<=domainHeight; y=y+0.75)//y++)
            {
                if (DIM == 2) {
                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                    nodeNum++;
                    numStromaNodes++;
                }
            }
            rownum++;
        }
        // Tumour Nodes
        unsigned numTumourNodes = 0;
        // Seed a cluster of tumour cells in the centre
        // Make sure it's not right at the edge
        double x = domainWidth*0.5;//0.8*domainWidth*RandomNumberGenerator::Instance()->ranf() + 0.1*domainWidth;
        double y = domainHeight*0.5;//0.8*domainHeight*RandomNumberGenerator::Instance()->ranf() + 0.1*domainHeight;
        nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x-0.25, y));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x, y-0.25));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x-0.25, y-0.25));
        nodeNum++;
        numTumourNodes = numTumourNodes+4;


        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        double avgCellCycleDurationTumour = 24;
        double avgCellCycleDurationStroma = 32;

        unsigned cancerLabelColour = 3;
        unsigned macrophageLabelColour = 7;
        MAKE_PTR_ARGS(CellLabel, p_stroma_label, (1));
        MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));
        for (unsigned i=0; i<numStromaNodes; i++)
        {
            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationStroma*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationStroma*1.25);
            p_model->SetHypoxicConcentration(0.2);
            p_model->SetNecroticConcentration(0.05);
            p_model->SetQuiescenceVolumeProportion(0.75);

//            ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
//            p_model->SetDimension(DIM);
//            p_model->SetQuiescentVolumeFraction(0.7);
//            p_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -2);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_stroma_label);
            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationStroma*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);

        }

        for (unsigned i=numStromaNodes; i<numStromaNodes+numTumourNodes; i++)
        {

            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationTumour*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationTumour*1.25);
            p_model->SetHypoxicConcentration(0.05);
            p_model->SetNecroticConcentration(0.05);
            // Tumour cells have much higher tolerance to contact inhibition
            p_model->SetQuiescenceVolumeProportion(0.6);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -1);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_cancer_label);

            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationTumour*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);
        }



        NodesOnlyMesh<DIM> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);


        // Make cell population (2D)
        NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

        // Write summary statistic files
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();
        cell_population.AddCellWriter<CellCycleClockWriter>();


        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<DIM> lower(0, 0, 0);
        ChastePoint<DIM> upper(domainWidth,domainHeight, 0);
        MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));


        //
        // OXYGEN
        //

        // We have Neumann BCs along the boundaries of the simulation domain
        // Dirichlet BCs are set at blood vessel point sources
        // Choose BVs to be random mesh points, rather than picking random points and then meshing
        std::vector<ChastePoint<DIM> > vesselLocations;
        // Under current settings, O2 diffuses approx 7 cell diameters before quiescence, approx 15 before necrosis
        // So ensure that centre of the domain is max of 15 cell diameters from a BV
        vesselLocations.push_back(ChastePoint<DIM> (37,14));
        vesselLocations.push_back(ChastePoint<DIM> (35,34));
        vesselLocations.push_back(ChastePoint<DIM> (33,6));
        vesselLocations.push_back(ChastePoint<DIM> (30,19));
        vesselLocations.push_back(ChastePoint<DIM> (26,32));
        vesselLocations.push_back(ChastePoint<DIM> (20,8));
        vesselLocations.push_back(ChastePoint<DIM> (19,33));
        vesselLocations.push_back(ChastePoint<DIM> (13,37));
        vesselLocations.push_back(ChastePoint<DIM> (10,7));
        vesselLocations.push_back(ChastePoint<DIM> (8,21));
        vesselLocations.push_back(ChastePoint<DIM> (5,31));
        vesselLocations.push_back(ChastePoint<DIM> (4,14));

        double consumptionRate = -0.03;//was -0.03
        double diffusionCoefficient = 1;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_oxygen_pde, (cell_population, consumptionRate,diffusionCoefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_oxygen_bc, (1));
        bool is_neumann_bc = false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        int updateIntervalForPdeInTimesteps = 120/2;
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>, p_oxygen_pde_modifier, (p_oxygen_pde, p_oxygen_bc, is_neumann_bc, p_cuboid, 1.0));
        p_oxygen_pde_modifier->SetTimestepInterval(1);
        p_oxygen_pde_modifier->SetVesselLocations(vesselLocations);
        p_oxygen_pde_modifier->SetDependentVariableName("oxygen");
        p_oxygen_pde_modifier->SetBcsOnBoxBoundary(false); //was false


        //
        // CXCL12
        //
        double cxcl12_decay_rate = -0.02;
        MAKE_PTR_ARGS(UniformSourceEllipticPde<DIM>, p_cxcl12_pde, (cxcl12_decay_rate));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_cxcl12_bc, (1));
        bool is_neumann_bc_cxcl12= false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>, p_cxcl12_pde_modifier, (p_cxcl12_pde, p_cxcl12_bc, is_neumann_bc_cxcl12, p_cuboid, 1.0));
        p_cxcl12_pde_modifier->SetTimestepInterval(1);
        p_cxcl12_pde_modifier->SetDependentVariableName("cxcl12");
        p_cxcl12_pde_modifier->SetBcsOnBoxBoundary(false); //was false
        p_cxcl12_pde_modifier->SetOutputGradient(true);
        p_cxcl12_pde_modifier->SetVesselLocations(vesselLocations);
        //
        // CSF-1
        //
        double dudtCoefficient = 1.0;
        double csfDiffusionCoefficient = 1.0;
        double sourceCoefficient = 0.25;
        double decayRate_csf = 0.02;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_csf1_pde, (cell_population, dudtCoefficient,csfDiffusionCoefficient,sourceCoefficient));
        p_csf1_pde->SetCellTypesAsSource({cancerLabelColour});
        p_csf1_pde->SetDecayRate(decayRate_csf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_csf1_bc, (0));
        bool is_neumann_bc_csf1 = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_csf1_pde_modifier, (p_csf1_pde, p_csf1_bc, is_neumann_bc_csf1, p_cuboid, 1.0));
        p_csf1_pde_modifier->SetDependentVariableName("csf1");
        p_csf1_pde_modifier->SetBcsOnBoxBoundary(true);
        p_csf1_pde_modifier->SetOutputGradient(true);

        //
        // TGF
        //
        double dudtCoefficient_tgf = 1.0;
        double diffusionCoefficient_tgf = 0.1;
        double sourceCoefficient_tgf = 0.1;
        double decayRate_tgf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_tgf_pde, (cell_population, dudtCoefficient_tgf,diffusionCoefficient_tgf,sourceCoefficient_tgf));
        p_tgf_pde->SetCellTypesAsSource({cancerLabelColour});
        p_tgf_pde->SetDecayRate(decayRate_tgf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_tgf_bc, (0));
        bool is_neumann_bc_tgf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_tgf_pde_modifier, (p_tgf_pde, p_tgf_bc, is_neumann_bc_tgf, p_cuboid, 1.0));
        p_tgf_pde_modifier->SetDependentVariableName("tgf");
        p_tgf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_tgf_pde_modifier->SetOutputGradient(true);

        //
        // EGF
        //
        double dudtCoefficient_egf = 1.0;
        double diffusionCoefficient_egf = 0.2;
        double sourceCoefficient_egf = 0.2;
        double decayRate_egf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited_phenotypeDependent<DIM>, p_egf_pde, (cell_population, dudtCoefficient_egf,diffusionCoefficient_egf,sourceCoefficient_egf));
        p_egf_pde->SetCellTypesAsSource({macrophageLabelColour});
        p_egf_pde->SetDecayRate(decayRate_egf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_egf_bc, (0));
        bool is_neumann_bc_egf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_egf_pde_modifier, (p_egf_pde, p_egf_bc, is_neumann_bc_egf, p_cuboid, 0.5));
        p_egf_pde_modifier->SetDependentVariableName("egf");
        p_egf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_egf_pde_modifier->SetOutputGradient(true);



        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<DIM> simulator(cell_population);
        simulator.AddSimulationModifier(p_oxygen_pde_modifier);
        simulator.AddSimulationModifier(p_csf1_pde_modifier);
        simulator.AddSimulationModifier(p_egf_pde_modifier);
        simulator.AddSimulationModifier(p_tgf_pde_modifier);
        simulator.AddSimulationModifier(p_cxcl12_pde_modifier);
        std::stringstream output_directory;
        output_directory << "Arwert2018/VI_PointVessels_i/";
        simulator.SetOutputDirectory(output_directory.str());
        simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
        double burnInPeriod = 24*7;
        double infiltrationDuration = 24*28;
        double endTime = burnInPeriod+infiltrationDuration;
        simulator.SetEndTime(endTime);


        MAKE_PTR(Arwert2018_AddMacrophagesFromPointVesselModifier<DIM>, p_addMacs_modifier);
        p_addMacs_modifier->SetTgfbThresholdForMacrophagePhenotypeSwitch(1);
        p_addMacs_modifier->SetTimeToAddMacrophages(burnInPeriod);
        p_addMacs_modifier->SetNumberOfMacrophagesToAdd(50);
        p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(infiltrationDuration);
        p_addMacs_modifier->SetVesselLocations(vesselLocations);
        simulator.AddSimulationModifier(p_addMacs_modifier);

//        MAKE_PTR(PhenotypeBasedMacrophagePhagocytosisModifier<DIM>, p_phagoytosis_modifier);
//        p_phagoytosis_modifier->SetCellTypesToEat({cancerLabelColour});
//        simulator.AddSimulationModifier(p_phagoytosis_modifier);


        // Add Brownian motion for all cells
//		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(0.01);
//		simulator.AddForce(p_diffusion_force);

        // Add Brownian motion for macrophages
        MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_forceMac);
        p_diffusion_forceMac->SetDiffusionScalingConstant(0.01);
        p_diffusion_forceMac->SetCellTypesForDiffusion({macrophageLabelColour});
        simulator.AddForce(p_diffusion_forceMac);

        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Macrophages go up csf1 gradient
        // todo make this dependent on phenotype close to 0?
        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCSF1);
        p_chemotaxis_MacrophageToCSF1->SetNutrient("csf1");
        p_chemotaxis_MacrophageToCSF1->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCSF1->SetSensitivity(1);
        simulator.AddForce(p_chemotaxis_MacrophageToCSF1);

        // Tumour cells go up egf gradient
        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_TumourCellToEGF);
        p_chemotaxis_TumourCellToEGF->SetNutrient("egf");
        p_chemotaxis_TumourCellToEGF->SetCellTypesForChemotaxis({cancerLabelColour});
        p_chemotaxis_TumourCellToEGF->SetSensitivity(2.5);
        simulator.AddForce(p_chemotaxis_TumourCellToEGF);

        // Macrophages go up cxcl12 gradient, depending on phenotype (close to 1)
        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCXCL12);
        p_chemotaxis_MacrophageToCXCL12->SetNutrient("cxcl12");
        p_chemotaxis_MacrophageToCXCL12->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCXCL12->SetSensitivity(1);
        simulator.AddForce(p_chemotaxis_MacrophageToCXCL12);


        MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(5.0);
        p_force->SetCutOffLength(1.5);
//        p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
//        p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        simulator.AddForce(p_force);
        MARK;

        // Killer which removes apoptotic cells
        MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
        simulator.AddCellKiller(p_apoptosis_killer);

        // Boundary conditions: Add a wall along the x axis
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        // Boundary conditions: Add a wall along the bottom
        c_vector<double,2> normal2 = zero_vector<double>(2);
        c_vector<double,2> point2 = zero_vector<double>(2);
        normal2(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        // Boundary conditions: Add a wall along the top
        c_vector<double,2> normal3 = zero_vector<double>(2);
        c_vector<double,2> point3 = zero_vector<double>(2);
        normal3(1) = 1.0;
        point3(1) = domainHeight;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point3, normal3));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        // Boundary conditions: Add a wall along the right hand side
        c_vector<double,2> point4 = zero_vector<double>(2);
        c_vector<double,2> normal4 = zero_vector<double>(2);
        point4(0) = domainWidth;
        normal4(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point4, normal4));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);



        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
//		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }

    void TestFinalSimulation() throw(Exception)
    {
        /*
         * Replicate the system described in Arwert et al (2018): https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5946803/pdf/main.pdf
         *
         * Macrophages "extravasate" with phenotype 0 due to CSF1
         * They're recruited up the CSF1 gradient into the tumour
         * Tumour expression of TGF-b causes their phenotype to change towards 1
         * (Interpreted as higher expression of CXCR4)
         * Macrophages with phenotype close to 1 are more sensitive to CXCL12 expressed by perivascular fibroblasts
         * Hence macrophages move back towards the vessel while also expressing EGF (hence attracting tumour cells)
         *
         * Change is unidirectional - i.e., phenotype can increase towards 1, but not away from it
         */

        /** The next line is needed because we cannot currently run node based simulations in parallel. */
        EXIT_IF_PARALLEL;
        SimulationTime::Instance()->SetStartTime(0.0);
        RandomNumberGenerator::Instance()->Reseed(time(NULL));
        CellPropertyRegistry::Instance()->Clear();
        CellId::ResetMaxCellId();

        int visualisationOutputFrequencyPerHour = 2;
        const int DIM = 2;

        // Generate Mesh:
        // Make Vector
        std::vector<Node<DIM>*> nodes;
        // Add some nodes

        // Domain size
        double domainHeight = 40;
        double domainWidth = 40;

        // Fill the domain with stromal cells

        // Stroma Nodes
        unsigned nodeNum=0;
        unsigned numStromaNodes = 0;
        // For staggering rows
        unsigned rownum = 1;
        double offset;
        for (double x=0; x<=domainWidth; x=x+0.75)//x++)
        {
            offset = 0;
            if (rownum % 2)
            {
                offset = 0.75/2;
            }
            for (double y=0+offset; y<=domainHeight; y=y+0.75)//y++)
            {
                if (DIM == 2) {
                    nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
                    nodeNum++;
                    numStromaNodes++;
                }
            }
            rownum++;
        }
        // Tumour Nodes
        unsigned numTumourNodes = 0;
        // Seed a cluster of tumour cells in the centre
        // Make sure it's not right at the edge
        double x = domainWidth*0.5;//0.8*domainWidth*RandomNumberGenerator::Instance()->ranf() + 0.1*domainWidth;
        double y = domainHeight*0.5;//0.8*domainHeight*RandomNumberGenerator::Instance()->ranf() + 0.1*domainHeight;
        nodes.push_back(new Node<DIM>(nodeNum, false, x, y));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x-0.25, y));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x, y-0.25));
        nodeNum++;
        nodes.push_back(new Node<DIM>(nodeNum, false, x-0.25, y-0.25));
        nodeNum++;
        numTumourNodes = numTumourNodes+4;


        // Make cell pointers
        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(StemCellProliferativeType, p_stem_type);
        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);

        double avgCellCycleDurationTumour = 24;
        double avgCellCycleDurationStroma = 32;

        unsigned cancerLabelColour = 3;
        unsigned macrophageLabelColour = 7;
        MAKE_PTR_ARGS(CellLabel, p_stroma_label, (1));
        MAKE_PTR_ARGS(CellLabel, p_cancer_label, (cancerLabelColour));
        for (unsigned i=0; i<numStromaNodes; i++)
        {
            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationStroma*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationStroma*1.25);
            p_model->SetHypoxicConcentration(0.2);
            p_model->SetNecroticConcentration(0.05);
            p_model->SetQuiescenceVolumeProportion(0.75);

//            ContactInhibitionCellCycleModel* p_model = new ContactInhibitionCellCycleModel;
//            p_model->SetDimension(DIM);
//            p_model->SetQuiescentVolumeFraction(0.7);
//            p_model->SetEquilibriumVolume(1.0);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -2);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_stroma_label);
            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationStroma*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);

        }

        for (unsigned i=numStromaNodes; i<numStromaNodes+numTumourNodes; i++)
        {

            GreenspanAndVolumeCellCycle* p_model = new GreenspanAndVolumeCellCycle;
            p_model->SetDimension(DIM);
            p_model->SetMinCellCycleDuration(avgCellCycleDurationTumour*0.75);
            p_model->SetMaxCellCycleDuration(avgCellCycleDurationTumour*1.25);
            p_model->SetHypoxicConcentration(0.05);
            p_model->SetNecroticConcentration(0.05);
            // Tumour cells have much higher tolerance to contact inhibition
            p_model->SetQuiescenceVolumeProportion(0.6);

            CellPtr p_cell(new Cell(p_state, p_model)); // Default cell radius is 0.5
            p_cell->SetCellProliferativeType(p_stem_type);
            p_cell->GetCellData()->SetItem("oxygen", 1);
            p_cell->GetCellData()->SetItem("csf1", 0);
            p_cell->GetCellData()->SetItem("egf", 0);
            p_cell->GetCellData()->SetItem("cxcl12", 0);
            p_cell->GetCellData()->SetItem("tgf", 0);
            p_cell->GetCellData()->SetItem("phenotype", -1);
            p_cell->SetApoptosisTime(48); // Apoptosis time in hours - how long before a cell is removed?
            p_cell->AddCellProperty(p_cancer_label);

            double birthTime = - RandomNumberGenerator::Instance()->ranf() *
                               (avgCellCycleDurationTumour*0.75);
            p_cell->SetBirthTime(birthTime);
            cells.push_back(p_cell);
        }



        NodesOnlyMesh<DIM> mesh;
        // Cut off length: 1.5 cell radii
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);


        // Make cell population (2D)
        NodeBasedCellPopulation<DIM> cell_population(mesh, cells);
        cell_population.SetAbsoluteMovementThreshold(DBL_MAX);//Set big movement threshold

        // Write summary statistic files
        cell_population.AddPopulationWriter<NodeLocationWriter>();
        cell_population.AddPopulationWriter<BoundaryNodeWriter>();
        cell_population.AddCellWriter<CellProliferativeTypesWriter>();
        cell_population.AddCellWriter<BoundaryCellWriter>();
        cell_population.AddCellWriter<CellCycleClockWriter>();


        // Create a ChasteCuboid on which to base the finite element mesh used to solve the PDE
        ChastePoint<DIM> lower(0, 0, 0);
        ChastePoint<DIM> upper(domainWidth,domainHeight, 0);
        MAKE_PTR_ARGS(ChasteCuboid<DIM>, p_cuboid, (lower, upper));


        //
        // OXYGEN
        //

        // We have Neumann BCs along the boundaries of the simulation domain
        // Dirichlet BCs are set at blood vessel point sources
        // Choose BVs to be random mesh points, rather than picking random points and then meshing
        std::vector<ChastePoint<DIM> > vesselLocations;
        // Under current settings, O2 diffuses approx 7 cell diameters before quiescence, approx 15 before necrosis
        // So ensure that centre of the domain is max of 15 cell diameters from a BV
        vesselLocations.push_back(ChastePoint<DIM> (37,14));
        vesselLocations.push_back(ChastePoint<DIM> (35,34));
        vesselLocations.push_back(ChastePoint<DIM> (33,6));
        vesselLocations.push_back(ChastePoint<DIM> (30,19));
        vesselLocations.push_back(ChastePoint<DIM> (26,32));
        vesselLocations.push_back(ChastePoint<DIM> (20,8));
        vesselLocations.push_back(ChastePoint<DIM> (19,33));
        vesselLocations.push_back(ChastePoint<DIM> (13,37));
        vesselLocations.push_back(ChastePoint<DIM> (10,7));
        vesselLocations.push_back(ChastePoint<DIM> (8,21));
        vesselLocations.push_back(ChastePoint<DIM> (5,31));
        vesselLocations.push_back(ChastePoint<DIM> (4,14));

        double consumptionRate = -0.03;//was -0.03
        double diffusionCoefficient = 1;
        MAKE_PTR_ARGS(AveragedSourceEllipticPde<DIM>, p_oxygen_pde, (cell_population, consumptionRate,diffusionCoefficient));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_oxygen_bc, (1));
        bool is_neumann_bc = false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        int updateIntervalForPdeInTimesteps = 120/2;
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>, p_oxygen_pde_modifier, (p_oxygen_pde, p_oxygen_bc, is_neumann_bc, p_cuboid, 1.0));
        p_oxygen_pde_modifier->SetTimestepInterval(1);
        p_oxygen_pde_modifier->SetVesselLocations(vesselLocations);
        p_oxygen_pde_modifier->SetDependentVariableName("oxygen");
        p_oxygen_pde_modifier->SetBcsOnBoxBoundary(false); //was false


        //
        // CXCL12
        //
        double cxcl12_decay_rate = -0.02;
        MAKE_PTR_ARGS(UniformSourceEllipticPde<DIM>, p_cxcl12_pde, (cxcl12_decay_rate));
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_cxcl12_bc, (1));
        bool is_neumann_bc_cxcl12= false; // Dirichlet BCs
        // Create a PDE modifier and set the name of the dependent variable in the PDE
        MAKE_PTR_ARGS(EllipticBoxDomainPdeModifier_VariableTimestep_PointVesselBCs<DIM>, p_cxcl12_pde_modifier, (p_cxcl12_pde, p_cxcl12_bc, is_neumann_bc_cxcl12, p_cuboid, 1.0));
        p_cxcl12_pde_modifier->SetTimestepInterval(1);
        p_cxcl12_pde_modifier->SetDependentVariableName("cxcl12");
        p_cxcl12_pde_modifier->SetBcsOnBoxBoundary(false); //was false
        p_cxcl12_pde_modifier->SetOutputGradient(true);
        p_cxcl12_pde_modifier->SetVesselLocations(vesselLocations);
        //
        // CSF-1
        //
        double dudtCoefficient = 1.0;
        double csfDiffusionCoefficient = 1.0;
        double sourceCoefficient = 0.25;
        double decayRate_csf = 0.02;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_csf1_pde, (cell_population, dudtCoefficient,csfDiffusionCoefficient,sourceCoefficient));
        p_csf1_pde->SetCellTypesAsSource({cancerLabelColour});
        p_csf1_pde->SetDecayRate(decayRate_csf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_csf1_bc, (0));
        bool is_neumann_bc_csf1 = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_csf1_pde_modifier, (p_csf1_pde, p_csf1_bc, is_neumann_bc_csf1, p_cuboid, 1.0));
        p_csf1_pde_modifier->SetDependentVariableName("csf1");
        p_csf1_pde_modifier->SetBcsOnBoxBoundary(true);
        p_csf1_pde_modifier->SetOutputGradient(true);

        //
        // TGF
        //
        double dudtCoefficient_tgf = 1.0;
        double diffusionCoefficient_tgf = 0.1;
        double sourceCoefficient_tgf = 0.1;
        double decayRate_tgf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited<DIM>, p_tgf_pde, (cell_population, dudtCoefficient_tgf,diffusionCoefficient_tgf,sourceCoefficient_tgf));
        p_tgf_pde->SetCellTypesAsSource({cancerLabelColour});
        p_tgf_pde->SetDecayRate(decayRate_tgf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_tgf_bc, (0));
        bool is_neumann_bc_tgf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_tgf_pde_modifier, (p_tgf_pde, p_tgf_bc, is_neumann_bc_tgf, p_cuboid, 1.0));
        p_tgf_pde_modifier->SetDependentVariableName("tgf");
        p_tgf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_tgf_pde_modifier->SetOutputGradient(true);

        //
        // EGF
        //
        double dudtCoefficient_egf = 1.0;
        double diffusionCoefficient_egf = 0.2;
        double sourceCoefficient_egf = 0.2;
        double decayRate_egf = 0.1;
        MAKE_PTR_ARGS(AveragedSourceParabolicPde_edited_phenotypeDependent<DIM>, p_egf_pde, (cell_population, dudtCoefficient_egf,diffusionCoefficient_egf,sourceCoefficient_egf));
        p_egf_pde->SetCellTypesAsSource({macrophageLabelColour});
        p_egf_pde->SetDecayRate(decayRate_egf);
        MAKE_PTR_ARGS(ConstBoundaryCondition<DIM>, p_egf_bc, (0));
        bool is_neumann_bc_egf = true;
        MAKE_PTR_ARGS(ParabolicBoxDomainPdeModifier<DIM>, p_egf_pde_modifier, (p_egf_pde, p_egf_bc, is_neumann_bc_egf, p_cuboid, 0.5));
        p_egf_pde_modifier->SetDependentVariableName("egf");
        p_egf_pde_modifier->SetBcsOnBoxBoundary(true);
        p_egf_pde_modifier->SetOutputGradient(true);



        /* We then pass in the cell population into an {{{OffLatticeSimulation}}},
         * (this time with dimension 3) and set the output directory, output multiple and end time. */
        OffLatticeSimulation<DIM> simulator(cell_population);
        simulator.AddSimulationModifier(p_oxygen_pde_modifier);
        simulator.AddSimulationModifier(p_csf1_pde_modifier);
        simulator.AddSimulationModifier(p_egf_pde_modifier);
        simulator.AddSimulationModifier(p_tgf_pde_modifier);
        simulator.AddSimulationModifier(p_cxcl12_pde_modifier);
        std::stringstream output_directory;
        output_directory << "Arwert2018/VII_Extravasation/";
        simulator.SetOutputDirectory(output_directory.str());
        simulator.SetSamplingTimestepMultiple(120/visualisationOutputFrequencyPerHour);
        double burnInPeriod = 24*7;
        double infiltrationDuration = 24*28;
        double endTime = burnInPeriod+infiltrationDuration;
        simulator.SetEndTime(endTime);


        MAKE_PTR(Arwert2018_AddMacrophagesFromPointVesselInResponseToCSFModifier<DIM>, p_addMacs_modifier);
        p_addMacs_modifier->SetTgfbThresholdForMacrophagePhenotypeSwitch(1);
        p_addMacs_modifier->SetTimeToAddMacrophages(burnInPeriod);
//        p_addMacs_modifier->SetNumberOfMacrophagesToAdd(50);
        p_addMacs_modifier->SetHalfMaximalExtravasationCsf1Conc(0.5);
        p_addMacs_modifier->SetMaximalProbOfExtravasationPerHour(0.05);
        p_addMacs_modifier->SetNumberOfHoursToRunSimulationAfterAddingMacrophages(infiltrationDuration);
        p_addMacs_modifier->SetVesselLocations(vesselLocations);
        simulator.AddSimulationModifier(p_addMacs_modifier);

//        MAKE_PTR(PhenotypeBasedMacrophagePhagocytosisModifier<DIM>, p_phagoytosis_modifier);
//        p_phagoytosis_modifier->SetCellTypesToEat({cancerLabelColour});
//        simulator.AddSimulationModifier(p_phagoytosis_modifier);


        // Add Brownian motion for all cells
//		MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_force);
//		p_diffusion_force->SetDiffusionScalingConstant(0.01);
//		simulator.AddForce(p_diffusion_force);

        // Add Brownian motion for macrophages
        MAKE_PTR(DiffusionForceChooseD<DIM>, p_diffusion_forceMac);
        p_diffusion_forceMac->SetDiffusionScalingConstant(0.01);
        p_diffusion_forceMac->SetCellTypesForDiffusion({macrophageLabelColour});
        simulator.AddForce(p_diffusion_forceMac);

        MAKE_PTR(VolumeTrackingModifier<DIM>, p_modifier);
        simulator.AddSimulationModifier(p_modifier);

        // Macrophages go up csf1 gradient
        // todo make this dependent on phenotype close to 0?
        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCSF1);
        p_chemotaxis_MacrophageToCSF1->SetNutrient("csf1");
        p_chemotaxis_MacrophageToCSF1->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCSF1->SetSensitivity(1);
        simulator.AddForce(p_chemotaxis_MacrophageToCSF1);

        // Tumour cells go up egf gradient
        MAKE_PTR(ChemotacticForce_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_TumourCellToEGF);
        p_chemotaxis_TumourCellToEGF->SetNutrient("egf");
        p_chemotaxis_TumourCellToEGF->SetCellTypesForChemotaxis({cancerLabelColour});
        p_chemotaxis_TumourCellToEGF->SetSensitivity(2.5);
        simulator.AddForce(p_chemotaxis_TumourCellToEGF);

        // Macrophages go up cxcl12 gradient, depending on phenotype (close to 1)
        MAKE_PTR(ChemotacticForce_PhenotypeDependent_SpecifyNutrientAndCellType<DIM>, p_chemotaxis_MacrophageToCXCL12);
        p_chemotaxis_MacrophageToCXCL12->SetNutrient("cxcl12");
        p_chemotaxis_MacrophageToCXCL12->SetCellTypesForChemotaxis({macrophageLabelColour});
        p_chemotaxis_MacrophageToCXCL12->SetSensitivity(1);
        simulator.AddForce(p_chemotaxis_MacrophageToCXCL12);


        MAKE_PTR(GeneralisedLinearSpringForceDifferentialAdhesionForApoptosisAndMacrophages<DIM>, p_force);
        p_force->SetMeinekeSpringStiffness(5.0);
        p_force->SetCutOffLength(1.5);
//        p_force->SetHeterotypicLabelledSpringConstantMultiplier(1.0);
//        p_force->SetHomotypicLabelledSpringConstantMultiplier(1.0);
        simulator.AddForce(p_force);
        MARK;

        // Killer which removes apoptotic cells
        MAKE_PTR_ARGS(ApoptoticCellKiller<DIM>, p_apoptosis_killer, (&cell_population));
        simulator.AddCellKiller(p_apoptosis_killer);

        // Boundary conditions: Add a wall along the x axis
        c_vector<double,2> point = zero_vector<double>(2);
        c_vector<double,2> normal = zero_vector<double>(2);
        normal(0) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc1, (&cell_population, point, normal));
        simulator.AddCellPopulationBoundaryCondition(p_bc1);

        // Boundary conditions: Add a wall along the bottom
        c_vector<double,2> normal2 = zero_vector<double>(2);
        c_vector<double,2> point2 = zero_vector<double>(2);
        normal2(1) = -1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc2, (&cell_population, point2, normal2));
        simulator.AddCellPopulationBoundaryCondition(p_bc2);

        // Boundary conditions: Add a wall along the top
        c_vector<double,2> normal3 = zero_vector<double>(2);
        c_vector<double,2> point3 = zero_vector<double>(2);
        normal3(1) = 1.0;
        point3(1) = domainHeight;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc3, (&cell_population, point3, normal3));
        simulator.AddCellPopulationBoundaryCondition(p_bc3);

        // Boundary conditions: Add a wall along the right hand side
        c_vector<double,2> point4 = zero_vector<double>(2);
        c_vector<double,2> normal4 = zero_vector<double>(2);
        point4(0) = domainWidth;
        normal4(0) = 1.0;
        MAKE_PTR_ARGS(PlaneBoundaryCondition<2>, p_bc4, (&cell_population, point4, normal4));
        simulator.AddCellPopulationBoundaryCondition(p_bc4);



        /* To run the simulation, we call {{{Solve()}}}. */
        simulator.Solve();
//		CellBasedSimulationArchiver<DIM, OffLatticeSimulation<DIM> >::Save(&simulator);

        /* To avoid memory leaks, we conclude by deleting any pointers that we created in the test.*/
        for (unsigned i=0; i<nodes.size(); i++)
        {
            delete nodes[i];
        }
    }
};

#endif //JB_CHASTE_TESTARWERT2018_HPP
