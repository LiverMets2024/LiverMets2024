#include "TypeDependentForce.hpp"
#include "NeutrophilProliferativeType.hpp"

	template<unsigned DIM>
	TypeDependentForce<DIM>::TypeDependentForce()
	    : AbstractForce<DIM>(),
	          mMovementParameter(0.01)
	{
	}
	
	template<unsigned DIM>
	TypeDependentForce<DIM>::~TypeDependentForce()
	{
	}
	
template<unsigned DIM>
	void TypeDependentForce<DIM>::SetMovementParameter(double movementParameter)
	{
	    assert(movementParameter > 0.0);
	    mMovementParameter = movementParameter;
	}
	
	template<unsigned DIM>
	double TypeDependentForce<DIM>::GetMovementParameter()
	{
	    return mMovementParameter;
	}
	
	template<unsigned DIM>
	void TypeDependentForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
	{
	    double dt = SimulationTime::Instance()->GetTimeStep();
	
	    // Iterate over the nodes
	    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
	         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
	         ++node_iter)
	    {
	                c_vector<double, DIM> force_contribution;
	        for (unsigned i=0; i<DIM; i++)
	        {
	            /*
42	             * The force on this cell is scaled with the timestep such that when it is
43	             * used in the discretised equation of motion for the cell, we obtain the
44	             * correct formula
45	             *
46	             * x_new = x_old + sqrt(2*D*dt)*W
47	             *
48	             * where W is a standard normal random variable.
49	             */

				// edit here for new force modifier
	            double xi = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
	
	            force_contribution[i] = (sqrt(2.0*mMovementParameter*dt)/dt)*xi;
	        }
	        node_iter->AddAppliedForceContribution(force_contribution);
	    }
	}
	
	template<unsigned DIM>
	void TypeDependentForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
	{
	    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";
	
	    // Call direct parent class
	    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
	}
	
	/////////////////////////////////////////////////////////////////////////////
	// Explicit instantiation
	/////////////////////////////////////////////////////////////////////////////
	
	template class TypeDependentForce<1>;
	template class TypeDependentForce<2>;
	template class TypeDependentForce<3>;
	
	// Serialization for Boost >= 1.36
	#include "SerializationExportWrapperForCpp.hpp"
	EXPORT_TEMPLATE_CLASS_SAME_DIMS(TypeDependentForce)