#include "TCellProliferativeType.hpp"

TCellProliferativeType::TCellProliferativeType()
    : AbstractCellProliferativeType(10)
{}

#include "SerializationExportWrapperForCpp.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(TCellProliferativeType)