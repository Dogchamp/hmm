USAGE:
    -i: <filename (in cwd)>
    -n: <P(N,N)>
    -g: <P(G,G)>

HMM
===
states : { N (Noncoding), P (Coding) }
transition probabilities : {  N->N: x, N->P: 1-x, P->P: y, P->N: 1-y }
emission probabilities : {
    N: 1/64 for all codons,
    P: row-wise probabilities of given tables
    }

