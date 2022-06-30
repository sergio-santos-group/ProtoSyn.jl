# iGBR-NN

All details regarding the usage and inner workings of iGBR-NN are explained in the model's paper: 
https://academic.oup.com/bioinformatics/article/36/6/1757/5613804?login=true

**iGBR-NN is a machine learning model to predict Born radii from local environment descriptors.** All models were trained using the Keras framework, in Python. As such, the Fogolari et al. researchers make available the `.h5` and `.json` files for the pre-trained models.

For faster implementation in ProtoSyn, these files were converted to `.onnx` format using the `keras2onnx` tool:
https://github.com/Cuda-Chen/keras2onnx-example

A modified version of ONNX.jl package is used to load these files in ProtoSyn, available from:
https://github.com/JosePereiraUA/ONNX.jl/tree/ops-fix