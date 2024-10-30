# Trellis_QEC_decoder

This repository implements the decoder from the semester project "Tensor network decoding for stabilizer codes based on the syndrome trellis representation". The code can be split into two parts

The first part is a python code, which constructs the trellis as a 2d array from a set of stabilizers. The constructor has an inbuilt function to generate the stabilizers for surface and color codes. The trellis is saved as a .txt file along with a full stabilizer tableau, used by the decoder. All pauli operators in denoted as vectors in the symplectic vector space.

The decoder implements the TN decoder from the paper in the function "trellis_mps_constructor_function.jl". The TN is implemented as with the Tensor datatype from the ITensor package in the julia programming language. Once the resulting Tn is constructed and contracted into an MPS-form as described in the paper, the tensor network is turend into a MPS datatype. An option can be given to truncate the MPS for a given maximum bond dimension. The actual decoding process is done in the "MPS_decoder_function.jl" function, where the resulting MPS contraction done for a given syndrome, returning a probability distrubution over logical classes.

As of right now the decoder has some numerical instability found for surface code of distances d>11. This issue has not yet been resolved
