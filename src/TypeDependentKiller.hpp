#ifndef TYPEDEPENDENTKILLER_HPP_
#define TYPEDEPENDENTKILLER_HPP_

#include "AbstractCellKiller.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Simple cell killer than induces apoptosis if a cell is surrounded by labelled cells.
 * A cell is denoted 'dead' using a label
  */
template<unsigned DIM>
class TypeDependentKiller : public AbstractCellKiller<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractCellKiller<DIM> >(*this);
       
    }

public:

    /**
     * Default constructor.
     *
     * @param pCellPopulation pointer to the cell population
     */
    TypeDependentKiller(AbstractCellPopulation<DIM>* pCellPopulation);

    void CheckAndLabelSingleCellForApoptosis(CellPtr pCell);

    /**
     * Loop over cells and start apoptosis randomly, based on the user-set
     * probability.
     */
      void CheckAndLabelCellsForApoptosisOrDeath();

    /**
     * Overridden OutputCellKillerParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputCellKillerParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(TypeDependentKiller)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a TargetedCellKiller.
 */
template<class Archive, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const TypeDependentKiller<DIM> * t, const unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractCellPopulation<DIM>* const p_cell_population = t->GetCellPopulation();
    ar << p_cell_population;
    
}

/**
 * De-serialize constructor parameters and initialise a TargetedCellKiller.
 */
template<class Archive, unsigned DIM>
inline void load_construct_data(
    Archive & ar, TypeDependentKiller<DIM> * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractCellPopulation<DIM>* p_cell_population;
    ar >> p_cell_population;
    

    // Invoke inplace constructor to initialise instance
    ::new(t)TypeDependentKiller<DIM>(p_cell_population);
}
}
} // namespace ...

#endif /*TYPEDEPENDENTKILLER_HPP_*/
