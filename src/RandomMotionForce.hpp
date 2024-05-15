#ifndef RANDOMMOTIONFORCE_HPP_
#define RANDOMMOTIONFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A force class to model random cell movement.
 */
template<unsigned DIM>
class RandomMotionForce : public AbstractForce<DIM>
{
private :

    /**
     * Random Movement Parameter.
     */
    double mMovementParameter;
        /**
     * Archiving.
     */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mMovementParameter;
    }

public :

    /**
45	     * Constructor.
46	     */
    RandomMotionForce();

    /**
50	     * Destructor.
        */
    ~RandomMotionForce();

    /**
     * Set the diffusion constant for the cells.
     *
     * @param movementParameter the movement parameter to use
     */
    void SetMovementParameter(double movementParameter);

    /**
     * Get the random motion coefficient.
     *
     * @return mMovementParameter
     */
    double GetMovementParameter();

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the tissue
     *
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(RandomMotionForce)

#endif /*RANDOMMOTIONFORCE_HPP_*/