#include "WeightsCalculator.h"

#include <set>

#include "common/AcquisitionBuffer.h"
#include "common/grappa_common.h"

#ifdef USE_CUDA
#include "gpu/WeightsCore.h"
#endif

#include "cpu/WeightsCore.h"

#include "SliceAccumulator.h"
#include "Unmixing.h"

#include "Gadget.h"

namespace {
    using namespace Gadgetron;
    using namespace Gadgetron::Core;
    namespace Grappa = Gadgetron::Grappa;

    // A similar function should be available in the std library at some point.
    template <class T, std::size_t N>
    std::array<std::remove_cv_t<T>, N> to_array(T (&a)[N])
    {
        std::array<std::remove_cv_t<T>, N> array{};
        std::copy(std::begin(a), std::end(a), array.begin());
        return array;
    }

    std::vector<Grappa::Slice> take_available_slices(InputChannel<Grappa::Slice> &input) {

        std::vector<Grappa::Slice> slices{};

        slices.emplace_back(input.pop());

        while(auto opt_slice = input.try_pop()) {
            slices.emplace_back(std::move(*opt_slice));
        }

        GDEBUG_STREAM("Read " << slices.size() << " available slice(s).");

        return std::move(slices);
    }

    class SupportMonitor {
    public:
        explicit SupportMonitor(const Context &context) {
            auto e_limits = context.header.encoding[0].encodingLimits;
            auto slices = e_limits.slice ? e_limits.slice->maximum + 1u : 1u;

            regions = std::vector<std::array<uint16_t, 4>>(slices, {0,0,std::numeric_limits<uint16_t>::max(),0});
        }

        void operator()(const Grappa::AnnotatedAcquisition &acquisition) {
            auto old_region = regions[slice_of(acquisition)];
            regions[slice_of(acquisition)] = std::array<uint16_t, 4> {
                0,
                uint16_t(samples_in(acquisition) - 1u),
                std::min(old_region[2], line_of(acquisition)),
                std::max(old_region[3], line_of(acquisition))
            };
        }

        std::array<uint16_t, 4> region_of_support(size_t slice) const {
            return regions[slice];
        }

        void clear(size_t slice) {
            regions[slice] = {0,0,std::numeric_limits<uint16_t>::max(),0};
        }

    private:
        std::vector<std::array<uint16_t, 4>> regions;
    };

    class AccelerationMonitor {
    public:
        AccelerationMonitor(const Context &context) {
            auto e_limits = context.header.encoding[0].encodingLimits;
            auto slices = e_limits.slice ? e_limits.slice->maximum + 1u : 1u;

            previous_line = std::vector<optional<size_t>>(slices, none);
            acceleration = std::vector<optional<size_t>>(slices, none);
        }

        void operator()(const Grappa::AnnotatedAcquisition &acquisition) {

            if(previous_line[slice_of(acquisition)]) {
                if (line_of(acquisition) < previous_line[slice_of(acquisition)].value()) {
                    acceleration[slice_of(acquisition)] = none;
                }
                else {
                    acceleration[slice_of(acquisition)] = line_of(acquisition) - previous_line[slice_of(acquisition)].value();
                }
            }
            previous_line[slice_of(acquisition)] = line_of(acquisition);
        }

        size_t acceleration_factor(size_t slice) const {
            return acceleration[slice].value();
        }

        void clear(size_t slice) {
            previous_line[slice] = acceleration[slice] = none;
        }

    private:
        std::vector<optional<size_t>> previous_line;
        std::vector<optional<size_t>> acceleration;
    };

    class DirectionMonitor {
    public:
        explicit DirectionMonitor(Grappa::AcquisitionBuffer &buffer, SupportMonitor &support, AccelerationMonitor &acceleration)
        : buffer(buffer), support(support), acceleration(acceleration) {
            position = read_dir = phase_dir = slice_dir = {0.0, 0.0, 0.0};
        }

        void operator()(const Grappa::AnnotatedAcquisition &acquisition) {

            auto header = std::get<ISMRMRD::AcquisitionHeader>(acquisition);

            if (position == to_array(header.position) &&
                read_dir == to_array(header.read_dir) &&
                phase_dir == to_array(header.phase_dir) &&
                slice_dir == to_array(header.slice_dir)) {
                return;
            }

            position = to_array(header.position);
            read_dir = to_array(header.read_dir);
            phase_dir = to_array(header.phase_dir);
            slice_dir = to_array(header.slice_dir);

            clear(slice_of(acquisition));
        }

        void clear(size_t slice) {
            buffer.clear(slice);
            support.clear(slice);
            acceleration.clear(slice);
        }


    private:
        Grappa::AcquisitionBuffer &buffer;
        SupportMonitor &support;
        AccelerationMonitor &acceleration;

        std::array<float, 3> position, read_dir, phase_dir, slice_dir;
    };
}

namespace Gadgetron::Grappa {

    template<class WeightsCore>
    Grappa::Weights create_weights(
            uint16_t index,
            const AcquisitionBuffer &buffer,
            uint16_t n_combined_channels,
            uint16_t n_uncombined_channels,
            const SupportMonitor &support_monitor,
            const AccelerationMonitor &acceleration_monitor,
            WeightsCore &core
    ) {
        return Grappa::Weights{
                {
                        index,
                        n_combined_channels,
                        n_uncombined_channels
                },
                core.calculate_weights(
                        buffer.view(index),
                        support_monitor.region_of_support(index),
                        acceleration_monitor.acceleration_factor(index),
                        n_combined_channels,
                        n_uncombined_channels
                )
        };
    }

    template<class WeightsCore>
    WeightsCalculator<WeightsCore>::WeightsCalculator(
            const Context &context,
            const std::unordered_map<std::string, std::string> &props
    ) : ChannelGadget<Grappa::Slice>(props), context(context) {}

    template<class WeightsCore>
    void WeightsCalculator<WeightsCore>::process(InputChannel<Grappa::Slice> &in, OutputChannel &out) {

        std::set<uint16_t> updated_slices{};
        uint16_t n_combined_channels = 0, n_uncombined_channels = 0;

        AcquisitionBuffer buffer{context};
        SupportMonitor support_monitor{context};
        AccelerationMonitor acceleration_monitor{context};

        buffer.add_pre_update_callback(DirectionMonitor{buffer, support_monitor, acceleration_monitor});
        buffer.add_post_update_callback([&](auto &acq) { updated_slices.insert(slice_of(acq)); });
        buffer.add_post_update_callback([&](auto &acq) { acceleration_monitor(acq); });
        buffer.add_post_update_callback([&](auto &acq) { support_monitor(acq); });
        buffer.add_post_update_callback([&](auto &acq) {
            n_combined_channels = combined_channels(acq);
            n_uncombined_channels = uncombined_channels(acq);
        });

        WeightsCore core{
                {coil_map_estimation_ks, coil_map_estimation_power},
                {block_size_samples, block_size_lines, convolution_kernel_threshold}
        };

        while (true) {
            auto slices = take_available_slices(in);
            buffer.add(slices);

            for (auto index : updated_slices) {

                if (!buffer.is_fully_sampled(index)) continue;

                out.push(create_weights(
                        index,
                        buffer,
                        n_combined_channels,
                        n_uncombined_channels,
                        support_monitor,
                        acceleration_monitor,
                        core
                ));
            }
            updated_slices.clear();
        }
    }

    using cpuWeightsCalculator = WeightsCalculator<CPU::WeightsCore>;
    GADGETRON_GADGET_EXPORT(cpuWeightsCalculator);

#ifdef USE_CUDA
    using gpuWeightsCalculator = WeightsCalculator<GPU::WeightsCore>;
    GADGETRON_GADGET_EXPORT(gpuWeightsCalculator);
#endif
}