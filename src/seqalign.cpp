#include <iostream>
#include <getopt.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>


/*

Function Definitions found towards the bottom of the file. :)

*/


////////////////////////////////////////////////////////////////////
// Object in Dynamic Programming Matrix
struct State {
    
public:

    int M;    // Mis/Match
    int D;    // Delete (gap in seq1, L)
    int I;    // Insert (gap in seq2, U)
    int max;  // Max Value

    // Constructors 
    State() : M(0), D(0), I(0) {}
    State(int _a) : M(_a), D(_a), I(_a), max(_a) {}
    State(int _m, int _d, int _i) : M(_m), D(_d), I(_i) {}

    void set_max() {
        max = std::max(std::max(M, D), I);
    }

    // Get Possible Paths
    std::vector<char> get_trace() {
        std::vector<char> max_states;
        if (M == max) { max_states.push_back('M'); }
        if (D == max) { max_states.push_back('D'); }
        if (I == max) { max_states.push_back('I'); }
        return max_states;
    }

    // Overload output operator for debugging
    friend std::ostream& operator<<(std::ostream& os, const State& state) {
        os << state.max;
        return os;
    }

};


////////////////////////////////////////////////////////////////////
// Alignment Class
class Alignment {

private:

    int n, m;
    int score = 0;

    std::string seq1, seq2;
    std::vector<std::string> tb;
    std::vector<std::vector<State>> dp;

    void traceback(const int &i, const int &j, const int x, std::vector<std::string> &alns);
    int align_itr(const int &match, const int &mismatch, const int &gap);
    
public:

    Alignment() {}

    Alignment(const std::string &s1, const std::string &s2) : seq1(s1), seq2(s2) {
        n = seq1.size();
        m = seq2.size();
        dp = std::vector<std::vector<State>>(n + 1, std::vector<State>(m + 1));
    }

    void reset() {
        score = 0;
        tb.clear();
        dp = std::vector<std::vector<State>>(n + 1, std::vector<State>(m + 1));
    }

    int get_score() const { return score; }
    int set_max(const int &a, const int &b, const int &c) {
        return std::max(std::max(a, b), c);
    }

    int align(const std::string &method, const std::string &algo);
    void print_alignment() const;
    void print_dp_matrix() const;
};


/////////////////////////////////////////////////////////////
// Utility Functions
void print_help();
std::vector<std::string> read_fasta(const std::string &file);

/////////////////////////////////////////////////////////////
// Main
int main(int argc, char **argv) {

    std::string type, mode;

    // Define long options
    static struct option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"type", required_argument, nullptr, 't'},
        {"mode", required_argument, nullptr, 'm'},
        {nullptr, 0, nullptr, 0}
    };

    int _t, _m;
    int opt;
    int option_index = 0;

    // Parse optional arguments
    while ((opt = getopt_long(argc, argv, "t:m:", long_options, &option_index)) != -1) {
        
        switch (opt) {
            case 'h':
                print_help();
                exit(1);

            case 't':
                _t = std::atoi(optarg);
                if (_t < 0 || _t > 3) {
                    std::cerr << "Invalid --type value. Must be 0, 1, 2, or 3.\n";
                    print_help();
                    return 1;
                } else if (_t == 0) {
                    type = "nw";
                } else if (_t == 1) {
                    type = "sg";
                } else {
                    type = "sw";
                }
                break;

            case 'm':
                _m = std::atoi(optarg);
                if (_m < 0 || _m > 1) {
                    std::cerr << "Invalid --mode value. Must be 0 or 1.\n";
                    print_help();
                    return 1;
                } else if (_m == 0) {
                    mode = "itr";
                } else {
                    mode = "rec";
                }
                break;

            default:
                print_help();
                return 1;
        }
    }

    // Remaining arguments are fasta1 and fasta2
    if (argc - optind != 2) {
        std::cerr << "Error: Two fasta files must be specified.\n";
        print_help();
        return 1;
    }

    std::vector<std::string> seq1 = read_fasta(argv[optind]);
    std::vector<std::string> seq2 = read_fasta(argv[optind + 1]);

    Alignment aln(seq1[0], seq2[0]);
    int score = aln.align(type, mode);
        
    aln.print_alignment();

    return 0;
}


/////////////////////////////////////////////////////////////
// Help Message
void print_help() {
    std::cout << "usage: seqalign [-h] [--type {0,1,2,3}] [--mode {0,1}] fasta1 fasta2\n\n"
              << "positional arguments:\n"
              << "  fasta1            fasta file 1 (X Sequence)\n"
              << "  fasta2            fasta file 2 (Y Sequence)\n\n"
              << "options:\n"
              << "  -h, --help        show this help message and exit\n"
              << "  --type {0,1,2,3}  the type of alignment: 0=global (default), 1=semiglobal_y, 2=semiglobal_x, 3=local\n"
              << "  --mode {0,1}      the mode of alignment: 0=iterative (default), 1=recursive\n";
}


// Process FASTA file into vector of sequences
std::vector<std::string> read_fasta(const std::string &file) {

    std::string seq, id, line;
    std::vector<std::string> sequences;

    std::ifstream fa(file);
    if (fa.is_open()) {
   
        // Iterate through file
        while (std::getline(fa, line)) {
            if (line[0] == '>') {
                if (id != "") {
                    sequences.push_back(seq);
                    seq = "";
                }
                id = line;
            } else {
                seq += line;
            }
        }
        sequences.push_back(seq);
        fa.close(); // Close the file

    } else {
        std::cerr << "ERROR: Could not read fasta file: " << file << "\n";
        throw "ERROR: Make sure fasta file exists.";
    }

    return sequences;
}



/////////////////////////////////////////////////////////////
// Perform Traceback
void Alignment::traceback(const int &i, const int &j, const int x, std::vector<std::string> &alns) {

    // Initilize Alignments Vector
    if (alns.empty()) {
        alns.push_back("");
    }

    // End of Alignment Cases
    if (i <= 0 && j <= 0) {
        return;
    
    } else if (i == 0 && j > 0) {
        alns[x] += "D";
        traceback(i, j - 1, x, alns);
        return;
    
    } else if (i > 0 && j == 0) {
        alns[x] += "I";
        traceback(i - 1, j, x, alns);
        return;
    } 


    int t;
    std::string t_str = alns[x]; 
    std::vector<char> states = dp[i][j].get_trace();

    // Recurse for all maximum states
    for (int y = 0; y < states.size(); y++) {

        t = x;  

        // If multiple paths
        if (y > 0) {
            alns.push_back(t_str);
            t = alns.size() - 1;
        }
        
        if (states[y] == 'M') {
            alns[t] += "M";
            traceback(i - 1, j - 1, t, alns);
        
        } else if (states[y] == 'I') {
            alns[t] += "I";
            traceback(i - 1, j, t, alns);
        
        } else if (states[y] == 'D') {
            alns[t] += "D";
            traceback(i, j - 1, t, alns);
        }

    }   
}


/////////////////////////////////////////////////////////////
// Generic Dynamic Programming Implementation
int Alignment::align_itr(const int &match, const int &mismatch, const int &gap) {

    int m;

    for (int i = 1; i <= seq1.size(); i++) {    
        for (int j = 1; j <= seq2.size(); j++) {

            m = match;
            if (seq1[i - 1] != seq2[j - 1]) {
                m = mismatch;
            }

            dp[i][j].M = Alignment::set_max(
                            dp[i-1][j-1].M + m,
                            dp[i-1][j-1].D + m,
                            dp[i-1][j-1].I + m
                          );
            dp[i][j].D = dp[i][j-1].max + gap;
            dp[i][j].I = dp[i-1][j].max + gap;
            dp[i][j].set_max();
        }
    }

    return dp[seq1.size()][seq2.size()].max;
}


/////////////////////////////////////////////////////////////
// General Alignment Call
int Alignment::align(const std::string &method, const std::string &algo) {
        
    int match;
    int mismatch;
    int gap;

    // Initialize Scoring Penalities and Boundary Conditions
    if (method == "nw") {

        // Needleman Wunsch
        match = 1;
        mismatch = -1;
        gap = -2;

        for (int i = 0; i <= n; i++) { dp[i][0] = State(gap * i); }
        for (int j = 0; j <= m; j++) { dp[0][j] = State(gap * j); }

    } else if (method == "sg") {

        // Semi-global
        match = 1;
        mismatch = -1;
        gap = -2;

        // TO DO
        if (seq1.size() > seq2.size()) {
            for (int j = 0; j <= seq2.size(); j++) { dp[0][j] = State(gap * j); }
        } else {
            for (int i = 0; i <= seq1.size(); i++) { dp[i][0] = State(gap * i); } 
        }

    } else {
        std::cerr << "Error: Unknown alignment method " << method << "\n";
        exit(1);
    }

    // Perform Alignment
    if (algo == "rec") {
        ;
    } else {
        score = align_itr(match, mismatch, gap);
    }
    

    // print_dp_matrix();

    // Traceback Alignments
    int num_alignments = 0;
    traceback(n, m, num_alignments, this -> tb);

    return score;
}


/////////////////////////////////////////////////////////////
// Print Dynamic Programming Matrix
void Alignment::print_dp_matrix() const {
    std::cerr << "X\tX\t";
    for (int x = 0; x < seq2.size(); x++) { std::cerr << seq2[x] << "\t"; }
    std::cerr << "\n";
    for (int i = 0; i <= seq1.size(); i++) {
        if (i == 0) {
            std::cerr << "X\t";
        } else {
            std::cerr << seq1[i - 1] << "\t";
        }
        for (int j = 0; j <= seq2.size(); j++) {
            if (i == 0) {
                std::cerr << dp[i][j].D << "\t";
            } else if (j == 0) {
                std::cerr << dp[i][j].I << "\t";
            } else {
                std::cerr << dp[i][j] << "\t";
            }
        }
        std::cerr << "\n";
    }
}


/////////////////////////////////////////////////////////////
// Print Alignment
void Alignment::print_alignment() const {

    int i, j;
    int x = 0;
    std::string aln_seq1;
    std::string aln_seq2;

    for (auto &t : tb) {
        
        i = 0; j = 0;
        aln_seq1 = ""; aln_seq2 = "";

        std::string str = t;
        std::reverse(str.begin(), str.end());

        for (const auto &c : str) {

            if (c == 'M') {
                aln_seq1 += seq1[i];
                aln_seq2 += seq2[j];
                i += 1; j += 1;

            } else if (c == 'I') {
                aln_seq1 += seq1[i];
                aln_seq2 += "_";
                i += 1;

            } else if (c == 'D') {
                aln_seq1 += "_";
                aln_seq2 += seq2[j];
                j += 1;
            }

        }

        std::cout << "Alignment #" << x + 1 << " (Score: " << this -> score << "):\n";
        std::cout << aln_seq1 << "\n" << aln_seq2 << "\n\n";
        x += 1;
    }
}