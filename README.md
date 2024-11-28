# Green's Function Total Field Inversion (gTFI) for Quantitative Susceptibility Mapping
![输入图片描述](README_md_files/17ce1000-ad89-11ef-8795-ffc8b7b71db1.jpeg?v=1&type=image)
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


