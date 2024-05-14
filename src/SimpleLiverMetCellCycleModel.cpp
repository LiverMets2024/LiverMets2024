#include "SimpleLiverMetCellCycleModel.hpp"
#include "TCellProliferativeType.hpp"
#include "BackgroundCellProliferativeType.hpp"


SimpleLiverMetCellCycleModel::SimpleLiverMetCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mMinCellCycleDurationTCell(12.0), // Hours
      mMaxCellCycleDurationTCell(14.0)  // Hours
{
}

SimpleLiverMetCellCycleModel::SimpleLiverMetCellCycleModel(const SimpleLiverMetCellCycleModel& rModel)
   : AbstractSimpleCellCycleModel(rModel),
     mMinCellCycleDurationTCell(rModel.mMinCellCycleDurationTCell),
     mMaxCellCycleDurationTCell(rModel.mMaxCellCycleDurationTCell)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variable mCellCycleDuration is initialized in the
     * AbstractSimpleCellCycleModel constructor.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     *
     * Note that mCellCycleDuration is (re)set as soon as
     * InitialiseDaughterCell() is called on the new cell-cycle model.
     */
}

AbstractCellCycleModel* SimpleLiverMetCellCycleModel::CreateCellCycleModel()
{
    return new SimpleLiverMetCellCycleModel(*this);
}

void SimpleLiverMetCellCycleModel::SetCellCycleDuration()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();
    //add modify for other cell types
    if (mpCell->GetCellProliferativeType()->IsType<TCellProliferativeType>())
    {
        mCellCycleDuration = mMinCellCycleDurationTCell + (mMaxCellCycleDurationTCell - mMinCellCycleDurationTCell) * p_gen->ranf(); 
    }
    else
    {
        mCellCycleDuration = DBL_MAX; // U[MinCCD,MaxCCD]
    }
}

double SimpleLiverMetCellCycleModel::GetMinCellCycleDurationForTCell()
{
    return mMinCellCycleDurationTCell;
}

void SimpleLiverMetCellCycleModel::SetMinCellCycleDurationForTCell(double minCellCycleDurationTCell)
{
    mMinCellCycleDurationTCell = minCellCycleDurationTCell;
}

double SimpleLiverMetCellCycleModel::GetMaxCellCycleDurationForTCell()
{
    return mMaxCellCycleDurationTCell;
}

void SimpleLiverMetCellCycleModel::SetMaxCellCycleDurationForTCell(double maxCellCycleDurationTCell)
{
    mMaxCellCycleDurationTCell = maxCellCycleDurationTCell;
}

double SimpleLiverMetCellCycleModel::GetAverageCellCycleTimeForTCell()
{
    return 0.5*(mMinCellCycleDurationTCell + mMaxCellCycleDurationTCell);
}

double SimpleLiverMetCellCycleModel::GetAverageTransitCellCycleTime()
{
    return 0.5*(mMinCellCycleDurationTCell + mMaxCellCycleDurationTCell);
}

double SimpleLiverMetCellCycleModel::GetAverageStemCellCycleTime()
{
    return 0.5*(mMinCellCycleDurationTCell + mMaxCellCycleDurationTCell);
}


void SimpleLiverMetCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MinCellCycleDurationTCell>" << mMinCellCycleDurationTCell << "</MinCellCycleDurationTCell>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDurationTCell>" << mMaxCellCycleDurationTCell << "</MaxCellCycleDurationTCell>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleLiverMetCellCycleModel)
