#include "TypeDependentKiller.hpp"
#include "CellLabel.hpp"
#include "ApoptoticCellProperty.hpp"
#include "Cell.hpp"
#include "AbstractCellPopulation.hpp"
#include "SmartPointers.hpp"
#include "CellLabel.hpp"

#include "TCellProliferativeType.hpp"
#include "MetCellProliferativeType.hpp"
#include "NeutrophilProliferativeType.hpp"
#include "FibroblastProliferativeType.hpp"
#include "BackgroundCellProliferativeType.hpp"

template <unsigned DIM>
TypeDependentKiller<DIM>::TypeDependentKiller(AbstractCellPopulation<DIM> *pCellPopulation)
    : AbstractCellKiller<DIM>(pCellPopulation)
{
}

template <unsigned SPACE_DIM>
void TypeDependentKiller<SPACE_DIM>::CheckAndLabelSingleCellForApoptosis(CellPtr pCell)
{
}

template <unsigned DIM>
void TypeDependentKiller<DIM>::CheckAndLabelCellsForApoptosisOrDeath()
{   
    const double age_tolerance = 30;//Just some arbitary death age

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->mpCellPopulation->Begin();
         cell_iter != this->mpCellPopulation->End();
         ++cell_iter)
    {
            // only iterate over the live cells
        if (!(cell_iter->template HasCellProperty<BackgroundCellProliferativeType>()))
        {   
            //add some interesting conditions to cause cell death 
            if (cell_iter->GetAge()>age_tolerance && !cell_iter->HasApoptosisBegun())
            {
                        cell_iter->StartApoptosis();
                        this->mpCellPopulation->Update();
            }
        }
    }
}

template <unsigned DIM>
void TypeDependentKiller<DIM>::OutputCellKillerParameters(out_stream &rParamsFile)
{

    // Call method on direct parent class
    AbstractCellKiller<DIM>::OutputCellKillerParameters(rParamsFile);
}

// Explicit instantiation
template class TypeDependentKiller<1>;
template class TypeDependentKiller<2>;
template class TypeDependentKiller<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TypeDependentKiller)
