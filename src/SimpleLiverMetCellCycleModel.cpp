#include "SimpleLiverMetCellCycleModel.hpp"
#include "TCellProliferativeType.hpp"
#include "BackgroundCellProliferativeType.hpp"
#include "MetCellProliferativeType.hpp"
#include "NeutrophilProliferativeType.hpp"
#include "FibroblastProliferativeType.hpp"

SimpleLiverMetCellCycleModel::SimpleLiverMetCellCycleModel()
    : AbstractSimpleCellCycleModel(),
      mMinCellCycleDurationTCell(12.0), // Hours
      mMaxCellCycleDurationTCell(14.0),  // Hours
      mMinCellCycleDurationMetCell(12.0), // Hours
      mMaxCellCycleDurationMetCell(14.0),  // Hours
      mMinCellCycleDurationNeutrophil(12.0), // Hours
      mMaxCellCycleDurationNeutrophil(14.0),  // Hours
      mMinCellCycleDurationFibroblast(12.0), // Hours
      mMaxCellCycleDurationFibroblast(14.0)  // Hours

{
}

SimpleLiverMetCellCycleModel::SimpleLiverMetCellCycleModel(const SimpleLiverMetCellCycleModel& rModel)
   : AbstractSimpleCellCycleModel(rModel),
     mMinCellCycleDurationTCell(rModel.mMinCellCycleDurationTCell),
     mMaxCellCycleDurationTCell(rModel.mMaxCellCycleDurationTCell),
     mMinCellCycleDurationMetCell(rModel.mMinCellCycleDurationMetCell),
     mMaxCellCycleDurationMetCell(rModel.mMaxCellCycleDurationMetCell),
     mMinCellCycleDurationNeutrophil(rModel.mMinCellCycleDurationNeutrophil),
     mMaxCellCycleDurationNeutrophil(rModel.mMaxCellCycleDurationNeutrophil),
     mMinCellCycleDurationFibroblast(rModel.mMinCellCycleDurationFibroblast),
     mMaxCellCycleDurationFibroblast(rModel.mMaxCellCycleDurationFibroblast)
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
    else if (mpCell->GetCellProliferativeType()->IsType<MetCellProliferativeType>())
    {
        mCellCycleDuration = mMinCellCycleDurationMetCell + (mMaxCellCycleDurationMetCell - mMinCellCycleDurationMetCell) * p_gen->ranf(); 
    }
    else if (mpCell->GetCellProliferativeType()->IsType<NeutrophilProliferativeType>())
    {
        mCellCycleDuration = mMinCellCycleDurationNeutrophil + (mMaxCellCycleDurationNeutrophil - mMinCellCycleDurationNeutrophil) * p_gen->ranf(); 
    }
    else if (mpCell->GetCellProliferativeType()->IsType<FibroblastProliferativeType>())
    {
        mCellCycleDuration = mMinCellCycleDurationFibroblast + (mMaxCellCycleDurationFibroblast - mMinCellCycleDurationFibroblast) * p_gen->ranf(); 
    }
    else
    {
        mCellCycleDuration = DBL_MAX; // U[MinCCD,MaxCCD]
    }
}


//T-cell cell cycle helper functions
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

//Met-cell cell cycle helper functions
double SimpleLiverMetCellCycleModel::GetMinCellCycleDurationForMetCell()
{
    return mMinCellCycleDurationMetCell;
}

void SimpleLiverMetCellCycleModel::SetMinCellCycleDurationForMetCell(double minCellCycleDurationMetCell)
{
    mMinCellCycleDurationMetCell = minCellCycleDurationMetCell;
}

double SimpleLiverMetCellCycleModel::GetMaxCellCycleDurationForMetCell()
{
    return mMaxCellCycleDurationMetCell;
}

void SimpleLiverMetCellCycleModel::SetMaxCellCycleDurationForMetCell(double maxCellCycleDurationMetCell)
{
    mMaxCellCycleDurationMetCell= maxCellCycleDurationMetCell;
}

double SimpleLiverMetCellCycleModel::GetAverageCellCycleTimeForMetCell()
{
    return 0.5*(mMinCellCycleDurationMetCell + mMaxCellCycleDurationMetCell);
}


//Neutrophil cell cycle helper functions
double SimpleLiverMetCellCycleModel::GetMinCellCycleDurationForNeutrophil()
{
    return mMinCellCycleDurationNeutrophil;
}

void SimpleLiverMetCellCycleModel::SetMinCellCycleDurationForNeutrophil(double minCellCycleDurationNeutrophil)
{
    mMinCellCycleDurationNeutrophil = minCellCycleDurationNeutrophil;
}

double SimpleLiverMetCellCycleModel::GetMaxCellCycleDurationForNeutrophil()
{
    return mMaxCellCycleDurationNeutrophil;
}

void SimpleLiverMetCellCycleModel::SetMaxCellCycleDurationForNeutrophil(double maxCellCycleDurationNeutrophil)
{
    mMaxCellCycleDurationNeutrophil= maxCellCycleDurationNeutrophil;
}

double SimpleLiverMetCellCycleModel::GetAverageCellCycleTimeForNeutrophil()
{
    return 0.5*(mMinCellCycleDurationNeutrophil + mMaxCellCycleDurationNeutrophil);
}


//Fibroblast cell cycle helper functions
double SimpleLiverMetCellCycleModel::GetMinCellCycleDurationForFibroblast()
{
    return mMinCellCycleDurationFibroblast;
}

void SimpleLiverMetCellCycleModel::SetMinCellCycleDurationForFibroblast(double minCellCycleDurationFibroblast)
{
    mMinCellCycleDurationFibroblast = minCellCycleDurationFibroblast;
}

double SimpleLiverMetCellCycleModel::GetMaxCellCycleDurationForFibroblast()
{
    return mMaxCellCycleDurationFibroblast;
}

void SimpleLiverMetCellCycleModel::SetMaxCellCycleDurationForFibroblast(double maxCellCycleDurationFibroblast)
{
    mMaxCellCycleDurationFibroblast= maxCellCycleDurationFibroblast;
}

double SimpleLiverMetCellCycleModel::GetAverageCellCycleTimeForFibroblast()
{
    return 0.5*(mMinCellCycleDurationFibroblast + mMaxCellCycleDurationFibroblast);
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
    *rParamsFile << "\t\t\t<MinCellCycleDurationMetCell" << mMinCellCycleDurationMetCell << "</MinCellCycleDurationMetCell\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDurationMetCell>" << mMaxCellCycleDurationMetCell << "</MaxCellCycleDurationMetCell\n";
    *rParamsFile << "\t\t\t<MinCellCycleDurationNeutrophil>" << mMinCellCycleDurationNeutrophil << "</MinCellCycleDurationNeutrophil>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDurationNeutrophil>" << mMaxCellCycleDurationNeutrophil << "</MaxCellCycleDurationNeutrophil>\n";
    *rParamsFile << "\t\t\t<MinCellCycleDurationFibroblast>" << mMinCellCycleDurationFibroblast << "</MinCellCycleDurationFibroblast>\n";
    *rParamsFile << "\t\t\t<MaxCellCycleDurationFibroblast>" << mMaxCellCycleDurationFibroblast << "</MaxCellCycleDurationFibroblast>\n";

    // Call method on direct parent class
    AbstractSimpleCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(SimpleLiverMetCellCycleModel)
