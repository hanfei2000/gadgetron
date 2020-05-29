#pragma once

#include "gadgetron_radtse_export.h"
#include "GenericReconBase.h"

namespace Gadgetron {

    class EXPORTGADGETS_RADTSE RadTseEchoShareGadget : public GenericReconDataBase
    {
    public:
        GADGET_DECLARE(RadTseEchoShareGadget);
        typedef GenericReconDataBase BaseClass;

        RadTseEchoShareGadget();
        ~RadTseEchoShareGadget();

        /// ------------------------------------------------------------------------------------
        /// parameters to control the reconstruction
        /// ------------------------------------------------------------------------------------
        GADGET_PROPERTY(direct_view_sharing,            int,    "direct view sharing output",                       0);
		GADGET_PROPERTY(tiered_view_sharing,            bool,   "Oversampling used in NFFT",                        1);
        GADGET_PROPERTY(tiered_view_sharing_us_factor,  double, "under-sampling factor for Tiered View Sharing",    1.0 );
		
    protected:
    
        // default interface function
        virtual int process_config(ACE_Message_Block* mb);
        virtual int process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1);

        virtual int directViewSharing(IsmrmrdReconData m, int num_shared_views);
        virtual int tieredViewSharing(IsmrmrdReconData m, double us);
    
    private:
        int num_etl, num_shots, num_radial_views, num_samples_per_spoke;
        double esp_us;
    };
}
