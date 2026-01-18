# Green's Function Total Field Inversion (gTFI) for Quantitative Susceptibility Mapping (QSM) in MRI
![gTFI](image/gTFI.png?v=1&type=image)
## Introduction

Quantitative susceptibility mapping (QSM) is a magnetic resonance imaging technique that quantifies tissue magnetic susceptibility by deconvolving the measured signal phase data. Accurate background field removal is essential for QSM, particularly in surface regions of the brain, such as the cerebral cortex, where background field interference can be substantial.

To address this, Green's Function Total Field Inversion (gTFI) method is proposed. The gTFI method models the background field using integral equations composed of Green's function and boundary conditions, eliminating the need for traditional filtering, assumptions, or regularization. gTFI simultaneously determines the background field at the boundary and the tissue susceptibility from the measured phase data, allowing for effective background field modeling, ensuring consistency across the data.

## Features

- Accurate modeling of the background field using Green's function and boundary conditions.
- Total field inversion approach without explicit background field removal, mitigating errors from boundary assumptions.
- Incorporates Dirichlet boundary conditions directly into the total field fitting, enhancing accuracy of QSM results in regions near boundaries.

## License

This project is licensed under the BSD 3-Clause License. See the [LICENSE](LICENSE) file for details.

### Third-Party Libraries
The `lib` folder contains third-party libraries, which are included for convenience. These libraries are licensed under their own terms. These libraries are not covered by this project's BSD 3-Clause license.

## Citation

If you use this code in your research, please cite the corresponding paper:

H. Zhong, G. Li, Y. Wang and J. Li, "Green’s Function Total Field Inversion for Quantitative Susceptibility Mapping," in IEEE Transactions on Medical Imaging, doi: 10.1109/TMI.2025.3639776.



