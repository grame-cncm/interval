#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "faust/dsp/dsp.h"
#include "faust/gui/UI.h"
#include "faust/gui/meta.h"

<<includeIntrinsic>>
<<includeclass>>

// Integration driver: initialized controls, deterministic nonzero audio input,
// and all output channels printed for comparison with a trusted Faust compiler.
int main()
{
    mydsp processor;
    processor.init(48000);
    constexpr int frames = 4096;
    std::vector<std::vector<FAUSTFLOAT>> inputs(processor.getNumInputs(),
                                               std::vector<FAUSTFLOAT>(frames, 0.25));
    std::vector<std::vector<FAUSTFLOAT>> outputs(processor.getNumOutputs(),
                                                std::vector<FAUSTFLOAT>(frames));
    std::vector<FAUSTFLOAT*> in, out;
    for (auto& channel : inputs) in.push_back(channel.data());
    for (auto& channel : outputs) out.push_back(channel.data());
    processor.compute(frames, in.data(), out.data());
    std::cout << std::setprecision(17);
    for (int sample = 0; sample < frames; ++sample) {
        for (const auto& channel : outputs) {
            if (!std::isfinite(channel[sample])) return 1;
            std::cout << channel[sample] << ' ';
        }
        std::cout << '\n';
    }
}
