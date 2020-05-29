#ifndef RadTseTrajGadget_H
#define RadTseTrajGadget_H
#pragma once

#include "gadgetron_radtse_export.h"
#include "Gadget.h"
#include "GadgetMRIHeaders.h"
#include "hoNDArray.h"

#include <ismrmrd/ismrmrd.h>
#include <complex>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <ismrmrd/xml.h>
#include <boost/optional.hpp>

namespace Gadgetron {

    class EXPORTGADGETS_RADTSE RadTseTrajGadget :
            public Gadget2<ISMRMRD::AcquisitionHeader, hoNDArray<std::complex<float> > > {

    public:
        GADGET_DECLARE(RadTseTrajGadget);

        RadTseTrajGadget();

        virtual ~RadTseTrajGadget();

    protected:

        virtual int process_config(ACE_Message_Block *mb);

        virtual int process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
                            GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2);

    private:
        long radial_reordering;
        long num_radial_spokes;
        long num_samples_per_spoke;
        long num_echo_train_length;
        long num_shots;
        hoNDArray<float> traj_and_dcf;
    };
}
#endif //SpiralToGenericGadget_H
