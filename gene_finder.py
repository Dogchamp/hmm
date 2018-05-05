import sys, getopt, os   # For file and cmdline argIO
import math              # For logspace functions

# ====================================== MODEL

class hmm(object):
    """
    For this particular project:
    states: Either noncoding or coding
    t_probs: p(N,N), p(N,G)=1-p(N,N), p(G,G), p(G,N)=1-p(G,G)
    e_probs: {from table}
    start_p: {noncoding: 1}
    """

    def __init__(self, states, observations, start_p, t_probs, e_probs):
        self.states = states            # Hidden States
        self.observations = observations
        self.start_p = start_p          # Always {N:1.0} for this project
        self.t_probs = t_probs          # Transition probabilty table
        self.e_probs = e_probs          # Emission probability table

    # Finds most likely sequence of HMM states that results in sequence of observations
    def viterbi(self, obs, states, start_p, trans_p, emit_p ):
        V = [{}]
        # Initialize base case
        for st in states:
            V[0][st] = {"prob": start_p[st] * emit_p[st][obs[0]], "prev": None}

        # Viterbi
        for t in range(1, len(obs)):
            V.append({})
            for st in states:
                max_tr_prob = max(V[t-1][prev_st]["prob"]*trans_p[prev_st][st] for prev_st in states)
                for prev_st in states:
                    if V[t-1][prev_st]["prob"] * trans_p[prev_st][st] == max_tr_prob:
                        max_prob = max_tr_prob * emit_p[st][obs[t]]
                        V[t][st] = {"prob": max_prob, "prev": prev_st}
                        break
        opt = []
        max_prob = max(value["prob"] for value in V[-1].values())
        previous = None
        for st, data in V[-1].items():
            if data["prob"] == max_prob:
                opt.append(st)
                previous = st
                break
        for t in range(len(V)-2, -1, -1):
            opt.insert(0, V[t+1][previous]["prev"])
            previous = V[t+1][previous]["prev"]
        print('Steps: ' + ' '.join(opt) + 'with probability of %s' % max_prob)

    # Finds most likely sequence of HMM states that results in sequence of observations
    def viterbi_logspace(self, obs, states, start_p, trans_p, emit_p ):
        for s in start_p:
            start_p[s] = math.log2(start_p[s])
        for from_s in trans_p:
            for to_s in trans_p[from_s]:
                trans_p[from_s][to_s] = math.log2(trans_p[from_s][to_s])
        for state in emit_p:
            for seq in emit_p[state]:
                emit_p[state][seq] = math.log2(emit_p[state][seq])
        V = [{}]
        # Initialize base case
        for st in states:
            V[0][st] = {"prob": start_p[st] + emit_p[st][obs[0]], "prev": None}

        # Viterbi
        for t in range(1, len(obs)):
            V.append({})
            for st in states:
                max_tr_prob = max(V[t-1][prev_st]["prob"] + trans_p[prev_st][st] for prev_st in states)
                for prev_st in states:
                    if V[t-1][prev_st]["prob"] + trans_p[prev_st][st] == max_tr_prob:
                        max_prob = max_tr_prob + emit_p[st][obs[t]]
                        V[t][st] = {"prob": max_prob, "prev": prev_st}
                        break
        opt = []
        max_prob = max(value["prob"] for value in V[-1].values())
        previous = None
        for st, data in V[-1].items():
            if data["prob"] == max_prob:
                opt.append(st)
                previous = st
                break
        for t in range(len(V)-2, -1, -1):
            opt.insert(0, V[t+1][previous]["prev"])
            previous = V[t+1][previous]["prev"]
        print('Steps: ' + ' '.join(opt) + 'with log score of %s' % max_prob)
        return opt

    def joint_prob_log(self, init_p, path, trans_p, emit_s, emit_p):
        total = init_p["Noncoding"]
        for i in range(len(path)):
            total += trans_p[path[i]][path[i]]
        for codon in emit_s:
            for state in emit_p:
                total += emit_p[state][codon]
        return total

# ===================================== END MODE
def load_file(gene_sequence_file):
    with open(gene_sequence_file) as file:
        data = file.read()
        desc, close_paren, sequence = data.partition(")")
        return (desc, close_paren, sequence)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# =========================================== MAIN ===========================================
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def main(argv):
    pn = 0.0
    pg = 0.0

    try:
        opts, args = getopt.getopt(argv, "i:n:g:")
    except getopt.GetoptError:
        print("USAGE: \ngene_finder.py -i <inputfile> -n <PNN probabilty> -g <PGG probability>")
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-i':
            filename = arg
        if opt == '-n':
            pn = arg
        if opt == '-g':
            pg = arg

    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
    filepath = os.path.join(__location__, filename)

    load_file(filepath)
    gene_sequence = load_file(filepath)[2]
    gene_sequence = gene_sequence.replace('\n', '')
    codon_sequence = [gene_sequence[i:i+3] for i in range(0,len(gene_sequence),3)]

    print(codon_sequence)
    print("pn: %f, pg: %f\n" % (float(pn), float(pg)))

    states = ('Noncoding', 'Coding')
    start_p = {'Noncoding': 0.99999, 'Coding': 0.00001}
    trans_p = {
            'Coding' : {'Coding': float(pg), 'Noncoding': 1.0-float(pg)},
            'Noncoding': {'Coding': 1.0-float(pn), 'Noncoding': float(pn)}
            }
    emit_p = {
            'Noncoding' : {"AAA": 1/64, "AAC": 1/64, "AAG": 1/64, "AAT": 1/64,  # AA
                           "ACA": 1/64, "ACC": 1/64, "ACG": 1/64, "ACT": 1/64,  # AC
                           "AGA": 1/64, "AGC": 1/64, "AGG": 1/64, "AGT": 1/64,  # AG
                           "ATA": 1/64, "ATC": 1/64, "ATG": 1/64, "ATT": 1/64,  # AT
                           "CAA": 1/64, "CAC": 1/64, "CAG": 1/64, "CAT": 1/64,  # CA
                           "CCA": 1/64, "CCC": 1/64, "CCG": 1/64, "CCT": 1/64,  # CC
                           "CGA": 1/64, "CGC": 1/64, "CGG": 1/64, "CGT": 1/64,  # CG
                           "CTA": 1/64, "CTC": 1/64, "CTG": 1/64, "CTT": 1/64,  # CT
                           "GAA": 1/64, "GAC": 1/64, "GAG": 1/64, "GAT": 1/64,  # GA
                           "GCA": 1/64, "GCC": 1/64, "GCG": 1/64, "GCT": 1/64,  # GC
                           "GGA": 1/64, "GGC": 1/64, "GGG": 1/64, "GGT": 1/64,  # GG
                           "GTA": 1/64, "GTC": 1/64, "GTG": 1/64, "GTT": 1/64,  # GT
                           "TAA": 1/64, "TAC": 1/64, "TAG": 1/64, "TAT": 1/64,  # TA
                           "TCA": 1/64, "TCC": 1/64, "TCG": 1/64, "TCT": 1/64,  # TC
                           "TGA": 1/64, "TGC": 1/64, "TGG": 1/64, "TGT": 1/64,  # TG
                           "TTA": 1/64, "TTC": 1/64, "TTG": 1/64, "TTT": 1/64   # TT

                        },

                           #A              C              G               T
            'Coding'    : {"AAA": 77/184,  "AAC": 66/184, "AAG": 37/184,  "AAT": 4/184,  # AA
                           "ACA": 3/98,    "ACC": 63/98,  "ACG": 13/98,   "ACT": 19/98,  # AC
                           "AGA": 1/28,    "AGC": 23/28,  "AGG": 1/28,    "AGT": 3/28,   # AG
                           "ATA": 1/188,   "ATC": 98/188, "ATG": 60/188,  "ATT": 29/188, # AT
                           "CAA": 15/116,  "CAC": 23/116, "CAG": 73/116,  "CAT": 5/116,  # CA
                           "CCA": 11/76,   "CCC": 1/76,   "CCG": 55/76,   "CCT": 9/76,   # CC
                           "CGA": 1/137,   "CGC": 46/137, "CGG": 1/137,   "CGT": 89/137, # CG
                           "CTA": 1/190,   "CTC": 18/190, "CTG": 141/190, "CTT": 11/190, # CT
                           "GAA": 147/338, "GAC": 85/338, "GAG": 46/338,  "GAT": 60/338, # GA
                           "GCA": 30/128,  "GCC": 19/128, "GCG": 49/128,  "GCT": 30/128, # GC
                           "GGA": 1/131,   "GGC": 47/131, "GGG": 5/131,   "GGT": 78/131, # GG
                           "GTA": 34/144,  "GTC": 21/144, "GTG": 34/144,  "GTT": 55/144, # GT
                           "TAA": 1/56,    "TAC": 37/56,  "TAG": 1/56,    "TAT": 17/56,  # TA
                           "TCA": 2/77,    "TCC": 38/77,  "TCG": 5/77,    "TCT": 32/77,  # TC
                           "TGA": 1/18,    "TGC": 4/18,   "TGG": 8/18,    "TGT": 5/18,   # TG
                           "TTA": 2/69,    "TTC": 44/69,  "TTG": 8/69,    "TTT": 15/69   # TT
                        }
            }
    h = hmm(codon_sequence, states, start_p, trans_p, emit_p)

    # Viterbi and Joint Prob
    path = h.viterbi_logspace(codon_sequence, states, start_p, trans_p, emit_p)
    print("joint prob: %f" % h.joint_prob_log(start_p, path, trans_p, codon_sequence, emit_p))

    is_gene = False
    genes = []
    start = 0
    end = 0

    # Find the genes from path
    for i in range(0, len(path)):
        if path[i] == 'Coding' and not is_gene:
            is_gene = True
            start = i
        if path[i] == 'Noncoding' and is_gene:
            is_gene = False
            end = i-1
            genes.append((start,end))
    print(genes)

if __name__ == "__main__":
    main(sys.argv[1:])
