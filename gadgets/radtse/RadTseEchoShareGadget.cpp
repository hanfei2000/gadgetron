
#include "RadTseEchoShareGadget.h"
#include <iomanip>

#include "hoNDArray_reductions.h"
#include "mri_core_def.h"
#include "hoNDArray_utils.h"
#include "hoNDArray_elemwise.h"
#include "hoNDFFT.h"
#include "mri_core_utility.h"

#include "RadTseTools.h"

namespace Gadgetron 
{
    RadTseEchoShareGadget::RadTseEchoShareGadget() 
    :   BaseClass()
    ,   num_etl(1)
    ,   num_shots(1)
    ,   num_radial_views(1)
    ,   num_samples_per_spoke(1)
    ,   esp_us(1000)
    {
    }

    RadTseEchoShareGadget::~RadTseEchoShareGadget()
    {
    }

    int RadTseEchoShareGadget::process_config(ACE_Message_Block* mb)
    {
        GADGET_CHECK_RETURN(BaseClass::process_config(mb) == GADGET_OK, GADGET_FAIL);

        ISMRMRD::IsmrmrdHeader h;
        ISMRMRD::deserialize(mb->rd_ptr(), h);

        ISMRMRD::TrajectoryDescription traj_desc;
        traj_desc = *h.encoding[0].trajectoryDescription;
        if (traj_desc.identifier != "RadialTSE932") {
            GDEBUG("This Gadget only supports Radial Tse data\n");
            return GADGET_FAIL;
        }   
        
        num_radial_views = traj_desc.userParameterLong[0].value;
        
        num_samples_per_spoke = h.encoding[0].encodedSpace.matrixSize.x;
        num_etl = h.userParameters->userParameterLong[64].value;
        num_shots = num_radial_views / num_etl;
        
        return GADGET_OK;
    }

    int RadTseEchoShareGadget::process(Gadgetron::GadgetContainerMessage< IsmrmrdReconData >* m1)
    {
         auto m_copy = *m1->getObjectPtr();
        // ---------------------
        // send out composite image
        // ---------------------
        if (this->next()->putq(m1) == -1)
        {
            GERROR("GenericReconKSpaceFilteringGadget::process, passing data on to next gadget");
            return GADGET_FAIL;
        }
        
        //-----------------------
        // send out a 4 folds undersampled image
        //-----------------------
        int directVS = direct_view_sharing.value();
        if(directVS){
            if(directViewSharing(m_copy,directVS)!=GADGET_OK)
                return GADGET_FAIL;
        }
        if(tiered_view_sharing.value()){
            double us = tiered_view_sharing_us_factor.value();
            if(tieredViewSharing(m_copy, us)!=GADGET_OK)
                return GADGET_FAIL;
        }
        return GADGET_OK;
    }
    GADGET_FACTORY_DECLARE(RadTseEchoShareGadget)
}

int RadTseEchoShareGadget::directViewSharing(IsmrmrdReconData m, int num_shared_views){

    std::vector<std::vector<int>> view;
    view.resize(num_etl);
    for (int i = 0; i < num_etl; i++)
        view[i].resize(num_shots);
    
    RadTse::newViewSelect(view, num_radial_views, num_etl);
    std::vector<int> selectedViews;
    
    for (int echo=0; echo<=num_etl-num_shared_views; echo++ )
    {
        selectedViews.clear();
        for(int i=echo; i<echo+num_shared_views; i++){
            for(int j=0; j<num_shots; j++)
                selectedViews.push_back(view[i][j]);
        }
        auto m2 = new GadgetContainerMessage<IsmrmrdReconData>();
        *m2->getObjectPtr() = m;
        auto & datasets = m2->getObjectPtr()->rbit_;
        auto & buffer = datasets[0];

        auto & dataOld = buffer.data_.data_;
        auto data_dims = dataOld.dimensions();
        auto dataNew_dims = data_dims;
        dataNew_dims[1] = selectedViews.size();
        hoNDArray < complex_float_t > dataNew = hoNDArray<complex_float_t>(dataNew_dims);
        auto oldDataPtr = dataOld.get_data_ptr();
        auto newDataPtr = dataNew.get_data_ptr();
        int cnt = 0;
        for (size_t v = 0; v <  data_dims[1]; v++){
            auto it = std::find (selectedViews.begin(), selectedViews.end(), v);
            if (it != selectedViews.end()){
                for(size_t c=0; c<data_dims[3]; c++){
                    for(size_t d=0; d<data_dims[0]; d++){
                        newDataPtr[c*data_dims[0]*dataNew_dims[1] + cnt * data_dims[0] + d] = oldDataPtr[c*data_dims[0]*data_dims[1] + v * data_dims[0] + d]; 
                    }
                }
                cnt ++;
            }
        }
        auto traj_dims = buffer.data_.trajectory_->dimensions();
        auto trajNew_dims = traj_dims;
        trajNew_dims[2] = selectedViews.size();
        Core::optional<hoNDArray<float>> trajNew = Core::optional<hoNDArray<float>>(trajNew_dims);
        float* trajOldPtr = buffer.data_.trajectory_->get_data_ptr();
        float* trajNewPtr = trajNew->get_data_ptr();
        cnt = 0;
        for (size_t v = 0; v <  traj_dims[2]; v++){
            auto it = std::find (selectedViews.begin(), selectedViews.end(), v);
            if (it != selectedViews.end()){
                for(size_t d=0; d<traj_dims[1]; d++){
                    trajNewPtr[cnt*traj_dims[1]*3 + d*3 + 0] = trajOldPtr[v*traj_dims[1]*3 + d*3 + 0]; 
                    trajNewPtr[cnt*traj_dims[1]*3 + d*3 + 1] = trajOldPtr[v*traj_dims[1]*3 + d*3 + 1]; 
                    trajNewPtr[cnt*traj_dims[1]*3 + d*3 + 2] = trajOldPtr[v*traj_dims[1]*3 + d*3 + 2]; 
                }
                cnt++;
            }    
        }
        //replace data
        dataNew_dims[0] = dataNew_dims[0]*dataNew_dims[1];
        dataNew_dims[1] = 1;
        dataNew.reshape(dataNew_dims);
        buffer.data_.data_ = dataNew;

        trajNew_dims[1] = trajNew_dims[1]*trajNew_dims[2];
        trajNew_dims[2] = 1;
        trajNew->reshape(trajNew_dims);
        buffer.data_.trajectory_ = trajNew;
        
        GDEBUG_STREAM("Sending out a copy echo images: " <<echo <<"\n");

        if (this->next()->putq(m2) == GADGET_FAIL)
            return GADGET_FAIL;
    }
    return GADGET_OK;
}

int RadTseEchoShareGadget::tieredViewSharing(IsmrmrdReconData m, double us){
    for (int echo=0; echo<num_etl; echo++ )
    {
        hoNDArray < double > mask = hoNDArray<double>({num_samples_per_spoke, num_radial_views});
        hoNDArray < double > dcf = hoNDArray<double>({num_samples_per_spoke, num_radial_views});
        auto maskPtr = mask.get_data_ptr();
        auto dcfPtr = dcf.get_data_ptr();
        for(int i=0; i<mask.get_number_of_elements(); i++){
            maskPtr[i] = 0;
            dcfPtr[i] = 0;
        }

        int numSamplesSelected;
        numSamplesSelected = mask.get_number_of_elements() - RadTse::prepareData(maskPtr, num_samples_per_spoke, num_radial_views, num_etl, echo, dcfPtr, us);
        auto m2 = new GadgetContainerMessage<IsmrmrdReconData>();
        *m2->getObjectPtr() = m;
        auto & datasets = m2->getObjectPtr()->rbit_;
        auto & buffer = datasets[0];

        auto & dataFull = buffer.data_.data_;
        auto data_dims = dataFull.dimensions();
        auto newData_dims = data_dims;
        newData_dims[0] = numSamplesSelected;
        newData_dims[1] = 1;
        
        hoNDArray < complex_float_t > dataNew = hoNDArray<complex_float_t>(newData_dims);
        auto dataFullPtr = dataFull.get_data_ptr();
        auto dataNewPtr = dataNew.get_data_ptr();
        int cnt = 0;
        for (size_t v = 0; v <  data_dims[1]; v++){
            for(size_t d = 0; d<data_dims[0]; d++){
                if(maskPtr[v*num_samples_per_spoke+d]<0.5){
                    for(size_t c=0; c<data_dims[3]; c++){
                        dataNewPtr[cnt+numSamplesSelected*c] = dataFull[c*data_dims[0]*data_dims[1] + v * data_dims[0] + d]; 
                    }
                    cnt++;
                }

            }
        }

        auto traj_dims = buffer.data_.trajectory_->dimensions();
        auto newTraj_dims = traj_dims;
        newTraj_dims[1] = numSamplesSelected;
        newTraj_dims[2] = 1;
        Core::optional<hoNDArray<float>> trajNew = Core::optional<hoNDArray<float>>(newTraj_dims);
        float* trajFullPtr = buffer.data_.trajectory_->get_data_ptr();
        float* trajNewPtr = trajNew->get_data_ptr();
        cnt = 0;
        for (size_t v = 0; v <  traj_dims[2]; v++){
            for(size_t d=0; d < traj_dims[1]; d++){
                if(maskPtr[v*num_samples_per_spoke+d]<0.5){
                    trajNewPtr[cnt*3 + 0] = trajFullPtr[v*traj_dims[1]*3 + d*3 + 0]; 
                    trajNewPtr[cnt*3 + 1] = trajFullPtr[v*traj_dims[1]*3 + d*3 + 1]; 
                    trajNewPtr[cnt*3 + 2] = dcfPtr[v*num_samples_per_spoke+d]; 
                    cnt++;
                }
            }    
        }

        //replace data
        buffer.data_.data_ = dataNew;
        buffer.data_.trajectory_ = trajNew;
        
        GDEBUG_STREAM("Sending out a copy echo images: " <<echo <<"\n");

        if (this->next()->putq(m2) == GADGET_FAIL)
            return GADGET_FAIL;
    }
    return GADGET_OK;
}