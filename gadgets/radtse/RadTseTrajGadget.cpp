#include "RadTseTrajGadget.h"
#include "ismrmrd/xml.h"

#include <algorithm>
#include <vector>

namespace Gadgetron {

    RadTseTrajGadget::RadTseTrajGadget()
    : radial_reordering(0)
    , num_radial_spokes(0)
    , num_samples_per_spoke(0)
    , num_echo_train_length(0)
    , num_shots(0) { }

    RadTseTrajGadget::~RadTseTrajGadget() {}

    int RadTseTrajGadget::process_config(ACE_Message_Block *mb) {
        // Start parsing the ISMRMRD XML header

        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(mb->rd_ptr(), h);

        if (h.encoding.size() != 1) {
            GDEBUG("This Gadget only supports one encoding space\n");
            return GADGET_FAIL;
        }   
        if (h.encoding[0].trajectory != ISMRMRD::TrajectoryType::RADIAL ){
            GDEBUG("This Gadget only supports radial sampling\n");
            return GADGET_FAIL;
        }
        ISMRMRD::TrajectoryDescription traj_desc;
        traj_desc = *h.encoding[0].trajectoryDescription;
        num_radial_spokes = traj_desc.userParameterLong[0].value;
        
        num_samples_per_spoke = h.encoding[0].encodedSpace.matrixSize.x;
        num_echo_train_length = h.userParameters->userParameterLong[64].value;
        num_shots = num_radial_spokes / num_echo_train_length;

        GDEBUG(" Sampling Specs: %d spokes, %d samples per spoke, etl: %d, %d shots", num_radial_spokes, num_samples_per_spoke, num_echo_train_length, num_shots);
        
        //other wip parameters ...
        //radial_reordering = h.userParameters->userParameterLong[0].value;
        //GDEBUG("Radial Reordering (wipalFree[0]):             %d\n", radial_reordering);

        return GADGET_OK;
    }

    int RadTseTrajGadget::
    process(GadgetContainerMessage<ISMRMRD::AcquisitionHeader> *m1,
            GadgetContainerMessage<hoNDArray<std::complex<float> > > *m2) {
        // Noise should have been consumed by the noise adjust, but just in case...
        //

        bool is_noise = m1->getObjectPtr()->isFlagSet(ISMRMRD::ISMRMRD_ACQ_IS_NOISE_MEASUREMENT);
        if (is_noise) {
            m1->release();
            return GADGET_OK;
        }
        ISMRMRD::AcquisitionHeader * hdr = m1->getObjectPtr();
        //hanfei...test
        //if(hdr->idx.slice>0){
        //   m1->release();
        //    return GADGET_OK;
        //}

        // Delete previously attached trajectories
        if (m2->cont()) {
            m2->cont()->release();
        }

        //calculate trajectory and dcf
        //
        size_t col = 512;
        size_t view = 128;
        double angle_offset = 0.0;
        
        //calculate_trajectory_and_dcf(col, view, angle_offset);
        //
        
        //temporary 
        std::vector<size_t> traj_dims = {3,col};
        hoNDArray < float > traj_and_dcf_single_line = hoNDArray<float>(traj_dims);

        double cur_theta = (hdr->idx.kspace_encode_step_1) * M_PI / num_radial_spokes;
        float kx_pos, ky_pos;
        for (int col = -num_samples_per_spoke/2; col < num_samples_per_spoke/2; col++){
            int idx = col+num_samples_per_spoke/2;
            //traj
            kx_pos = float(-col*sin(cur_theta)/num_samples_per_spoke);
            ky_pos = float(col*cos(cur_theta)/num_samples_per_spoke); 

            if(abs(kx_pos) >0.5 || abs(ky_pos) >0.5)
                col = 512;

            traj_and_dcf_single_line[idx*3] = kx_pos;          //kx
            traj_and_dcf_single_line[idx*3+1] =  ky_pos;         //ky            
            
            //dcf
            if(col==0)
                traj_and_dcf_single_line[idx*3+2] = float(M_PI / (4*num_radial_spokes));           //dcf
            else
                traj_and_dcf_single_line[idx*3+2] = float(M_PI * abs(col)/num_radial_spokes); 

        }
    
        //finalize
        GadgetContainerMessage<hoNDArray<float> > *cont = new GadgetContainerMessage<hoNDArray<float> >();
        *(cont->getObjectPtr()) = traj_and_dcf_single_line;
        m2->cont(cont);

        //We need to make sure that the trajectory dimensions are attached.
        m1->getObjectPtr()->trajectory_dimensions = 3;

        if (this->next()->putq(m1) < 0) {
            GDEBUG("Failed to put job on queue.\n");
            return GADGET_FAIL;
        }

        return GADGET_OK;
    }

    int 
    RadTseTrajGadget::calculate_trajectory_and_dcf(unsigned int col, unsigned short view, double angle_offset)
    {   
        std::vector<size_t> traj_dims = {3,col*view};
        traj_and_dcf = hoNDArray<float>(traj_dims);

        for (size_t i = 0; i < col*view; i++){
            traj_and_dcf[i*3] = float(i + .0);
            traj_and_dcf[i*3+1] = float(i + 0.1);
            traj_and_dcf[i*3+2] = float(1 + 0.2);
        }

        return 1;
    } 

    GADGET_FACTORY_DECLARE(RadTseTrajGadget)
}
