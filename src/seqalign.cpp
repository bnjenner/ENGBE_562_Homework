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
// Object in Dynamic Programming F Matrix
struct State {
    
public:

    int M;           // Mis/Match
    int D;           // Delete (gap in seq1, L)
    int I;           // Insert (gap in seq2, U)
    int max;         // Max Value
    bool visited;    // Visited

    // Constructors 
    State() : M(0), D(0), I(0), visited(false) {}
    State(int _a) : M(_a), D(_a), I(_a), max(_a), visited(true) {}
    State(int _m, int _d, int _i) : M(_m), D(_d), I(_i), visited(false) {}

    // Make State Decision
    void set_max() {
        max = std::max(std::max(M, D), I);
        visited = true;
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

    int n, m;                                  // Length of sequences

    int match = 1;                             // Match Reward
    int mismatch = -1;                         // Mismatch Penality
    int gap = -2;                              // Gap Penality

    int max_score = -1;                        // Max Alignment Score
    int max_i = -1, max_j = -1;                // Max Position
    int tb_i = -1, tb_j = -1;                  // Traceback Start Position

    bool flip = false;                         // Flip Sequences in Output (not necessary)

    std::string seq1, seq2;                    // Input Sequences
    std::vector<std::string> tb;               // Traceback Paths
    std::vector<std::vector<State>> dp;        // F Matrix

    //////////////////////////////////
    // Private Alignment Methods
    void align_itr(const std::string &method);
    void align_rec(const int i, const int j, const std::string &method);
    void traceback(const int i, const int j, const int x, std::vector<std::string> &alns, const std::string &method);
    
    // Private Scoring Utilities
    int get_max(const std::vector<int> vec, const std::string &method);
    int score_funct(const int &i, const int &j, const char c, const std::string &method);
    void record_max(const int &i, const int &j, const std::string &method) {
        // Keep track of max scores for sw and sg
        if (method == "sw" && max_score < dp[i][j].max) {
            max_score = dp[i][j].max;
            max_i = i, max_j = j;   
        } else if (method == "sg" && j == this -> m && max_score < dp[i][j].max) {
            max_score = dp[i][j].max;
            max_i = i;
        }
    }

public:

    //////////////////////////////////
    // Constructors 
    Alignment() {}
    Alignment(const std::string &s1, const std::string &s2) : seq1(s1), seq2(s2) {
        n = seq1.size(); m = seq2.size();
        if (m > n) { // Flip to keep seq2 smallest
            n = seq2.size(); m = seq1.size();
            seq1 = s2; seq2 = s1;
            flip = true;
        }
        dp = std::vector<std::vector<State>>(n + 1, std::vector<State>(m + 1));
    }

    //////////////////////////////////
    // Public Methods
    int align(const std::string &method, const std::string &algo);
    int get_score() const { return max_score; }
    void print_alignment() const;
    void reset() {
        max_score = -1;
        max_i = -1, max_j = -1;
        tb_i = -1, tb_j = -1;
        tb.clear();
        dp = std::vector<std::vector<State>>(n + 1, std::vector<State>(m + 1));
    }
};


/////////////////////////////////////////////////////////////
// Utility Functions (Defined Below)
void print_help();
std::vector<std::string> read_fasta(const std::string &file);


/////////////////////////////////////////////////////////////
// Main
int main(int argc, char **argv) {

    // Option Variables
    int _t, _m, opt;
    int option_index = 0;
    std::string type = "nw", mode = "itr";

    // Define long options
    static struct option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"type", required_argument, nullptr, 't'},
        {"mode", required_argument, nullptr, 'm'},
        {nullptr, 0, nullptr, 0}
    };

    // Parse optional arguments
    while ((opt = getopt_long(argc, argv, "t:m:", long_options, &option_index)) != -1) {   
        switch (opt) {
            case 'h':
                print_help();
                return 1;
            case 't':
                _t = std::atoi(optarg);
                if (_t < 0 || _t > 3) {
                    std::cerr << "Invalid --type value. Must be 0, 1, 2, or 3.\n";
                    print_help();
                    return 1;
                } else if (_t == 0) {
                    type = "nw";
                } else if (_t == 1 || _t == 2) {
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

    // Read in Fasta Files
    std::vector<std::string> seq1 = read_fasta(argv[optind]);
    std::vector<std::string> seq2 = read_fasta(argv[optind + 1]);

    // Create Alignment Object
    Alignment aln(seq1[0], seq2[0]);

    // Align and Print Alignments
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


////////////////////////////////////////////////////////////////////
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
// Get Max of Vector Utility
int Alignment::get_max(const std::vector<int> vec, const std::string &method) {
    int res = vec[0];
    if (method == "sw") { res = 0; } // Maximum Score in SW is 0.
    for (const auto &v: vec) { res = std::max(res, v); }
    return res;
}


/////////////////////////////////////////////////////////////
// F Matrix Scoring Function
int Alignment::score_funct(const int &i, const int &j, const char c, const std::string &method) {
    int score, s;
    if (c == 'M') {
        s = match;
        if (seq1[i - 1] != seq2[j - 1]) { s = mismatch; }   
        score = Alignment::get_max({
                                    dp[i-1][j-1].M + s,
                                    dp[i-1][j-1].D + s,
                                    dp[i-1][j-1].I + s
                                    },
                                   method);
    } else if (c == 'D') {
        score = Alignment::get_max({dp[i][j-1].max + gap}, method);
    } else {
        score = Alignment::get_max({dp[i-1][j].max + gap}, method);
    }
    return score;
}



/////////////////////////////////////////////////////////////
// Perform Traceback
void Alignment::traceback(const int i, const int j, const int x, 
                          std::vector<std::string> &alns, const std::string &method) {

    // Initilize Alignments Vector
    if (alns.empty()) { alns.push_back(""); }

    // Termination / Special Traceback Cases
    if (i <= 0 && j <= 0) {
        return; 
    } else if (i == 0 && j > 0) {
        alns[x] += "D";
        traceback(i, j - 1, x, alns, method);
        return;
    } else if (i > 0 && j == 0 || (max_i < i && method == "sg")) {
        if (method == "sw") { return; }
        alns[x] += "I";
        traceback(i - 1, j, x, alns, method);
        return;
    }

    int t;
    std::string t_str = alns[x]; 
    std::vector<char> states = dp[i][j].get_trace();

    // Recurse for all maximum states (possible paths)
    for (int y = 0; y < states.size(); y++) {

        // If multiple paths
        t = x;  
        if (y > 0) {
            alns.push_back(t_str);
            t = alns.size() - 1; // New aln string index
        }

        // Traceback other routes        
        if (states[y] == 'M') {
            alns[t] += "M";
            traceback(i - 1, j - 1, t, alns, method);
        } else if (states[y] == 'I') {
            alns[t] += "I";
            traceback(i - 1, j, t, alns, method);
        } else if (states[y] == 'D') {
            alns[t] += "D";
            traceback(i, j - 1, t, alns, method);
        }

    }   
}


/////////////////////////////////////////////////////////////
// Iterative Dynamic Programming Implementation
void Alignment::align_itr(const std::string &method) {

    for (int i = 1; i <= n; i++) {    
        for (int j = 1; j <= m; j++) {

            dp[i][j].M = Alignment::score_funct(i, j, 'M', method);
            dp[i][j].D = Alignment::score_funct(i, j, 'D', method);
            dp[i][j].I = Alignment::score_funct(i, j, 'I', method);    
           
            dp[i][j].set_max();
            Alignment::record_max(i, j, method);
        }
    }
}


/////////////////////////////////////////////////////////////
// Recursive Dynamic Programming Implementation
void Alignment::align_rec(const int i, const int j, const std::string &method) {

    int s;
    if (dp[i][j].visited == true) { return; }

    align_rec(i-1, j-1, method);
    align_rec(i, j-1, method);
    align_rec(i-1, j, method);

    dp[i][j].M = Alignment::score_funct(i, j, 'M', method);
    dp[i][j].D = Alignment::score_funct(i, j, 'D', method);
    dp[i][j].I = Alignment::score_funct(i, j, 'I', method);

    dp[i][j].set_max();
    Alignment::record_max(i, j, method);
}


/////////////////////////////////////////////////////////////
// General Alignment Call
int Alignment::align(const std::string &method, const std::string &algo) {
        
    tb_i = n;
    tb_j = m;

    // Initialize Scoring Penalities and Boundary Conditions
    if (method == "nw") {
        for (int i = 0; i <= n; i++) { dp[i][0] = State(gap * i); }
        for (int j = 0; j <= m; j++) { dp[0][j] = State(gap * j); }
        max_i = n;
        max_j = m;
    } else if (method == "sg") {
        for (int i = 0; i <= n; i++) { dp[i][0] = State(0); }
        for (int j = 0; j <= m; j++) { dp[0][j] = State(gap * j); }
        max_j = m;
    } else if (method == "sw") {
        for (int i = 0; i <= n; i++) { dp[i][0] = State(0); }
        for (int j = 0; j <= m; j++) { dp[0][j] = State(0); }
    } else {
        std::cerr << "Error: Unknown alignment method " << method << "\n";
        exit(1);
    }

    // Perform Alignment and get score
    if (algo == "rec") {
        align_rec(n, m, method);
    } else {
        align_itr(method);
    }
    max_score = dp[max_i][max_j].max;

    // Adjust Traceback Starts if needed
    if (method == "sw") { tb_i = max_i; tb_j = max_j; }
    if (method == "sg") { tb_i = n; }

    // Traceback Alignments
    int num_alignments = 0;
    traceback(tb_i, tb_j, num_alignments, this -> tb, method);

    return Alignment::get_score();
}

/////////////////////////////////////////////////////////////
// Print Alignment
void Alignment::print_alignment() const {

    int i, j;
    int x = 0;
    std::string cons;
    std::string aln_seq1;
    std::string aln_seq2;

    // For each possible alignment
    for (auto &t : tb) {
        
        cons = "";
        aln_seq1 = "";
        aln_seq2 = "";
        i = tb_i; j = tb_j;

        // Iterate through traceback strings
        for (const auto &c : t) {

            if (c == 'M') {
                aln_seq1 = seq1[i - 1] + aln_seq1;
                aln_seq2 = seq2[j - 1] + aln_seq2;
                if (seq1[i - 1] == seq2[j - 1]) {
                    cons = '.' + cons;
                } else {
                    cons = seq2[j - 1] + cons;
                    if (flip) { cons[0] = seq1[i - 1]; }
                }
                i -= 1; j -= 1;

            } else if (c == 'I') {
                aln_seq1 = seq1[i - 1] + aln_seq1;
                aln_seq2 = "_" + aln_seq2;
                cons = '-' + cons;    
                if (flip) { cons[0] = seq1[i - 1]; } // Replace is output is flipped
                i -= 1;

            } else if (c == 'D') {
                aln_seq1 = "_" + aln_seq1;
                aln_seq2 = seq2[j - 1] + aln_seq2;
                cons = seq2[j - 1] + cons;
                if (flip) { cons[0] = '-'; } // Replace is output is flipped
                j -= 1;
            }
        }

        // Report Score and Alignments
        std::cout << "Alignment #" << x + 1 << " (Score: " << this -> max_score << "):\n";
        if (flip) {
            std::cout << aln_seq2 << "\n" << aln_seq1 << "\n" << cons << "\n\n";
        } else {
            std::cout << aln_seq1 << "\n" << aln_seq2 << "\n" << cons << "\n\n";
        }
        x += 1;
    }
}