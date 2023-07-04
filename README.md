# SHECS PIR
In the src folder there are two codes: the main file is implementation of SHECS-PIR 2.0 i.e., implementation of SHECS-PIR wih homomorphic evaluation of trace function and the other code is implementation of SHECS-PIR with conversion algorithm using gadget decomposition and scalar product.

The  implementation is written in the Rust programming language. Please refer to https://www.rust-lang.org/tools/install to install the required tools (e.g., rustup, cargo and rustc). If after installaion of alll the tools, while building the code any crate is said to be missing, install that crate manually. 

The main library we use is concrete-core, but this is imported automatically with cargo. The code uses Panacea as a dependency. Implementation of Panacea can be refered to here: https://github.com/KULeuven-COSIC/Panacea 

To compile the code run the command: ```cargo build```

To run the code run the command: ```cargo run```

Please note that if the slower SHECS-PIR is run, it will take a long time (maybe about an hour or above) to run completely. 
