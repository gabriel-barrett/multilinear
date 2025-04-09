The witness will be a row-major trace.
Indices will be subdivided into row indices and column indices.
That means padding is needed both for rows and columns, as they
must be a power of 2.

To implement lookup arguments, the trace columns will have to be
divided into two. The pre-random trace and the post-random trace.

Let's say the pre-random trace goes from `0b0000...0b1010`. Suppose
the total trace is of width 2^10. Then, we can compute the commitment
of the pre-random trace by the first 10 elements which is a combination
of the merkle tree of size 8 and a merkle tree of size 2 (well ignoring
the blowup factor for now). This depends on the commitment scheme of
course.
