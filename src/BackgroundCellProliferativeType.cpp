#include "BackgroundCellProliferativeType.hpp"

BackgroundCellProliferativeType::BackgroundCellProliferativeType()
    : AbstractCellProliferativeType(11)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(BackgroundCellProliferativeType)