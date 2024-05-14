#ifndef SIMPLELIVERMETCELLCYCLEMODEL_HPP_
#define SIMPLELIVERMETCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A stochastic cell-cycle model where cells divide with a stochastic cell cycle duration
 * with the length of the cell cycle drawn from a uniform distribution
 * on [mMinCellCycleDuration, mMaxCellCycleDuration].
 *
 * If the cell is differentiated, then the cell cycle duration is set to be infinite,
 * so that the cell will never divide.
 */
class SimpleLiverMetCellCycleModel : public AbstractSimpleCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    /**
     * The minimum cell cycle duration. Used to define the uniform distribution.
     * Defaults to 12 hours.
     */
    double mMinCellCycleDurationTCell;

    /**
     * The maximum cell cycle duration. Used to define the uniform distribution.
     * Defaults to 14 hours.
     */
    double mMaxCellCycleDurationTCell;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;
        archive & mMinCellCycleDurationTCell;
        archive & mMaxCellCycleDurationTCell;
    }

protected:

    /**
     * Protected copy-constructor for use by CreateCellCycleModel().
     *
     * The only way for external code to create a copy of a cell cycle model
     * is by calling that method, to ensure that a model of the correct subclass is created.
     * This copy-constructor helps subclasses to ensure that all member variables are correctly copied when this happens.
     *
     * This method is called by child classes to set member variables for a daughter cell upon cell division.
     * Note that the parent cell cycle model will have had ResetForDivision() called just before CreateCellCycleModel() is called,
     * so performing an exact copy of the parent is suitable behaviour. Any daughter-cell-specific initialisation
     * can be done in InitialiseDaughterCell().
     *
     * @param rModel the cell cycle model to copy.
     */
    SimpleLiverMetCellCycleModel(const SimpleLiverMetCellCycleModel& rModel);

public:

    /**
     * Constructor - just a default, mBirthTime is set in the AbstractCellCycleModel class.
     */
    SimpleLiverMetCellCycleModel();

    /**
     * Overridden SetCellCycleDuration() method to add stochastic cell cycle times
     */
    void SetCellCycleDuration();

    /**
     * Overridden builder method to create new copies of
     * this cell-cycle model.
     *
     * @return new cell-cycle model
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @return mMinCellCycleDuration
     */
    double GetMinCellCycleDurationForTCell();

    /**
     * Set mMinCellCycleDuration.
     *
     * @param minCellCycleDurationTCell
     */
    void SetMinCellCycleDurationForTCell(double minCellCycleDurationTCell);

    /**
     * @return mMaxCellCycleDuration
     */
    double GetMaxCellCycleDurationForTCell();

    /**
     * Set mMaxCellCycleDuration.
     *
     * @param maxCellCycleDurationTCell
     */
    void SetMaxCellCycleDurationForTCell(double maxCellCycleDurationTCell);

    /**
     * Overridden GetAverageTransitCellCycleTime() method.
     *
     * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
     */
    double GetAverageCellCycleTimeForTCell();

        /**
     * Overridden GetAverageTransitCellCycleTime() method.
     *
     * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
     */
    double GetAverageTransitCellCycleTime();

    /**
     * Overridden GetAverageStemCellCycleTime() method.
     *
     * @return the average of mMinCellCycleDuration and mMaxCellCycleDuration
     */
    double GetAverageStemCellCycleTime();

 
    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(SimpleLiverMetCellCycleModel)

#endif /*SIMPLELIVERMETCELLCYCLEMODEL_HPP_*/
