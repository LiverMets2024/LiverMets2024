#include "MetCellProliferativeType.hpp"

MetCellProliferativeType::MetCellProliferativeType()
    : AbstractCellProliferativeType(12)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(MetCellProliferativeType)