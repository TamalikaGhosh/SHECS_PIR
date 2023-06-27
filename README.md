# SHECS PIR
In the src file there are two codes: the main file is implementation of SHECS-PIR with faster conversion algorith (with homomorphic evaluation of trace function) and the other uses slower conversion technique. 

The  implementation is written in the Rust programming language. Please refer to https://www.rust-lang.org/tools/install to install the required tools (e.g., rustup, cargo and rustc). The main library we use is concrete-core, but this is imported automatically with cargo. The code uses Panacea as a dependency. Implementation of Panacea can be refered to here: https://github.com/KULeuven-COSIC/Panacea 

To compile the code run the command: ```cargo build```

To run the code run the command: ```cargo run```
