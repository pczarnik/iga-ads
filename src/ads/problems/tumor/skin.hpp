#ifndef ADS_PROBLEMS_TUMOR_SKIN_HPP_
#define ADS_PROBLEMS_TUMOR_SKIN_HPP_


namespace ads {
namespace tumor {


class skin_model {
public:
    enum layer {
        stratum_corneum,
        stratum_spinosum,
        basement_membrame,
        dermis,
        hypodermis,
        count
    };

    static constexpr std::size_t layer_count = static_cast<std::size_t>(layer::count);

    double diffusion_coefficient[layer_count] = {
        2.5, // stratum corneum
        15,  // stratum spinosum
        0.1, // basement membrame
        7.5, // dermis
        2.5, // hypodermis
    };

    double layer_size[layer_count] = {
        0.1,  // stratum corneum
        0.15, // stratum spinosum
        0.05, // basement membrame
        0.5,  // dermis
        0.2,  // hypodermis
    };

    double stratum_corneum_top = 1;
    double stratum_spinosum_top = 1 - layer_size[stratum_corneum];
    double basement_membrame_top = stratum_spinosum_top - layer_size[stratum_spinosum];
    double dermis_top = basement_membrame_top - layer_size[basement_membrame];
    double hypodermis_top = dermis_top - layer_size[dermis];

    double top[layer_count] = {
        stratum_corneum_top,
        stratum_spinosum_top,
        basement_membrame_top,
        dermis_top,
        hypodermis_top
    };


    double diffusion(double /*x*/, double /*y*/, double z) const {
        for (std::size_t i = 1; i < layer_count; ++ i) {
            if (z > top[i]) {
                return diffusion_coefficient[i - 1];
            }
        }
        return diffusion_coefficient[hypodermis];
    }
};


}
}


#endif /* ADS_PROBLEMS_TUMOR_SKIN_HPP_ */
