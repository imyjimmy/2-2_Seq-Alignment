/* Data Structures */
import java.util.ArrayList;
import java.util.List;
//may not need above.
import java.util.Stack;

import java.util.Hashtable;
import java.util.Set;

/* For file input and outputs */
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

public class SequenceAligner{
    public Hashtable<String, String> sequences = new Hashtable<String, String>();
    public Hashtable<String, Integer> pam100 = new Hashtable<String, Integer>();

    //MSA
    public Hashtable<String, String> msa_pairs = new Hashtable<String, String>(); 
    //key1:key2 => "AATTCCGG"
    public Hashtable<String, Integer> msa_pairs_scores = new Hashtable<String, Integer>();
    //key1:key2 => 34, etc
    public Hashtable<String, Double> msa_pairs_distances = new Hashtable<String, Double>();
    //{(key1:key2, 3.04), (key1:key3, 3.13), etc
    
    public String[] msa_alignment; //lol

    public double avgOutgroupDistance = 0.0;
    public String outGroup;

    //some default values
    public int openGapPenalty = -5;
    public int gapExtensionPenalty = -1;
    public int matchScore = 1;
    public int mismatchPenalty = -1;

    public String seq_type;
    public int numSolutions = 0;
    public boolean multiseq;
    public String pamFile;

    //matrix, cell_origin, seq1, seq2 (?)

    public enum Direction {
        STRT, LEFT, TOP, DG, LD, TD, LT, ALL;
    }
    
    /*
     Todo: as the user for scoring parameters. As noted in the assignment.
     Particularly in regards to the affine gap penalty: A+(B⋅L)
     params: --interactive or -i
     --type or -t (n|N or p|P)
     --file or -f [filepath]
     -m or --match
     -n or --notmatch
     -o or --ogap
     -e or --egap
     -h or --help 
    */
    public static void main(String []args){
        Hashtable<String, String> parameters = new Hashtable<String, String>();
        SequenceAligner seq_aligner = new SequenceAligner();
        if (args.length > 0 && (args[0].equals("-h") || args[0].equals("--help"))) {
            System.out.println("params: -h or --help to print this screen.");
            System.out.println("--type or -t (n|N or p|P) for nucleotide or protein match");
            System.out.println("--file or -f [filepath]");
            System.out.println("-m or --match [score]");
            System.out.println("-n or --nonmatch [score]");
            System.out.println("-o or --ogap open gap penalty [score]");
            System.out.println("-e or --egap gap extension penalty [score]");
            System.exit(0);
        } else {
            for (int i = 0; i < args.length-1; i+=2) {
                parameters.put(args[i],args[i+1]);
            }    
            System.out.println(parameters);
            seq_aligner.start(parameters);    
        }
        
        Hashtable<String, String> s = seq_aligner.sequences;
        String[] keys = s.keySet().toArray(new String[s.size()]);

        // System.out.println(seq_aligner.seq_type);
        if (seq_aligner.seq_type.equals("P") || seq_aligner.seq_type.equals("p")) {
            seq_aligner.parsePAM();
        }
        
        if (keys.length == 2) {
            seq_aligner.multiseq = false;
            String seq1 = s.get(keys[1]);
            String seq2 = s.get(keys[0]);

            int[][] matrix = seq_aligner.createMatrix(seq1, seq2);
            Direction[][] cell_origin_matrix = new Direction[seq1.length()+1][seq2.length()+1];

            // seq_aligner.printMatrix(matrix, cell_origin_matrix, seq1, seq2);
            seq_aligner.initializeMatrix(matrix, cell_origin_matrix);

            //Begin the algorithm.
            seq_aligner.needlemanWunsch(matrix, cell_origin_matrix, seq1, seq2);

            //return the alignment.
            seq_aligner.getAlignment(matrix, cell_origin_matrix, keys[1], keys[0], seq1, seq2);
        } else { //Starting MSA Alignment
            seq_aligner.multiseq = true; //not sure why we need this...to tell prettyprint something?
            String consensus = "";

            keys = s.keySet().toArray(new String[s.size()]);
            while( keys.length > 1) {
                System.out.println("keys.length: " + keys.length);

                for (int i = 0; i < keys.length-1; i++) {
                    for (int j = i+1; j < keys.length; j++) {
                        String seq1 = s.get(keys[i]);
                        String seq2 = s.get(keys[j]);

                        int[][] matrix = seq_aligner.createMatrix(seq1, seq2);
                        Direction[][] cell_origin = new Direction[seq1.length()+1][seq2.length()+1];
                        System.out.println("aligning: " + keys[i] + " with " + keys[j]);
                        seq_aligner.initializeMatrix(matrix, cell_origin);
                        seq_aligner.needlemanWunsch(matrix, cell_origin, seq1, seq2);
                        seq_aligner.getAlignment(matrix, cell_origin, keys[i], keys[j], seq1, seq2);
                    }
                }

                String[] keypairs = seq_aligner.msa_pairs_scores.keySet().toArray(new String[seq_aligner.msa_pairs_scores.size()]);
               
                // System.out.println("");
                System.out.println("keypairs.length: " + keypairs.length);
                //find the max (highest scoring) keypair.
                int max = -10000;
                int max_i = -1;
                for (int i = 0; i < keypairs.length; i++) {
                    System.out.println("pairs: " + keypairs[i]); 
                    int score = seq_aligner.msa_pairs_scores.get(keypairs[i]).intValue();
                    int distance = seq_aligner.msa_pairs_distances.get(keypairs[i]).intValue();

                    System.out.println("score: " + score);
                    System.out.println("distance: " + distance);

                    //"In the MSA, you can use simple alignment scores.""
                    if (score > max) {
                        max = score;   
                        max_i = i;
                    }
                }

                System.out.println("max_i: " + max_i);
                String maxkeys = keypairs[max_i];
                String[] key = maxkeys.split(":");

                String alignment = seq_aligner.msa_pairs.get(maxkeys);
                String distanceStr = seq_aligner.msa_pairs_distances.get(maxkeys).toString();
                System.out.println("closest alignments: " + alignment + " from: " + maxkeys);

                consensus = seq_aligner.consensus(alignment);

                String newKey = "(";
                
                //removing old keys, adding new ones. System.out.printf("Value: %.2f", value);
                //DecimalFormat df = new DecimalFormat("#.00");
                //df.format(...);
                for (int i = 0; i < key.length; i++) {
                    s.remove(key[i]);

                    for (int j = 0; j < keypairs.length; j++) {
                        if (keypairs[j].contains(key[i])) {
                            seq_aligner.msa_pairs_scores.remove(keypairs[j]);
                            //seq_aligner.msa_pairs_distances.remove(keypairs[j]);
                            //probably remove, then make new distances.
                            seq_aligner.msa_pairs.remove(keypairs[j]);
                        }
                    }

                    newKey += key[i];
                    if (i != key.length -1){
                        newKey += ",";
                    } else if (i == key.length -1) {
                        newKey += " " + distanceStr + " )";
                    }
                }
                s.put(newKey, consensus);

                keys = s.keySet().toArray(new String[s.size()]);
            }

            String consensusKey = s.keySet().toArray(new String[s.size()])[0];
            System.out.println("Newick Tree Format: " + consensusKey + "\nConsensus: " + s.get(consensusKey));
        }
     }

    /* helper method
    */
    public void parseFile(String filename) {
        // System.out.println("public void parseFile(String filename) {");
        File file = new File(filename);
        try { 
            Scanner inputFile = new Scanner(file); 
            String key = "";
            String val = "";
            do {
                String line = inputFile.nextLine();
                // System.out.println(line);
                if (line.startsWith(">")) {
                    if (!"".equals(val) && !"".equals(key)) { //puts the previous key val pair into sequences.

                        if (this.seq_type.equals("n") || this.seq_type.equals("N")) {
                            if ( (val.contains("T") && val.contains("U")) || (!val.matches("[ATCGU]+")) ) {
                                System.out.println("sequence either contains both T and U, or contains illegal characters.");
                                System.exit(1);
                            } 
                            else {
                                sequences.put(key, val);    
                            }
                        } else { 
                            if (!val.matches("[ARNDCQEGHILKMFPSTWYV]+")) {
                                System.out.println("sequence contains illegal amino acid letters. Sorry.");
                                System.exit(1);
                            } else {
                                sequences.put(key, val);
                            }
                        }
                        
                        // System.out.println("putting into sequences: key is " + key + " val is " + val);
                        key = "";
                        val = "";
                    }
                    //new key
                    key = line.substring(1);
                    // System.out.println("key substring: " + key);
                } else {
                    val += line;
                }
            } while (inputFile.hasNextLine());

            if (!"".equals(val) && !"".equals(key)) { //avoids null pointer exceptions.
                sequences.put(key, val);
                // System.out.println("putting into sequences: key is " + key + " val is " + val);
                key = "";
                val = "";
            }
            /* 
            * Assumption: the sequence is now loaded.
            */

        } catch (FileNotFoundException fnfe) {
            fnfe.printStackTrace();
        }
    }

    public void parsePAM() {
        System.out.println("inside parsePam");
        File pam = new File(this.pamFile);
        String column = "";
        String line;
        //char[] column = new char[](20); //20 amino acids. hardcoded by nature since apx billions of years ago.
        try {
            Scanner inputFile = new Scanner(pam);
            int row = 0;
            do {
                line = inputFile.nextLine();
                // System.out.println("inside parsePam");
                // System.out.println(line);
                if (line.startsWith(" ")) {    //first line
                    line = line.replaceAll("\\s+","");
                    column = line;
                    // System.out.println(line);
                    // System.out.println(line.length());
                    // System.out.println(column);
                } else { // not the first line
                    char aa = line.charAt(0);
                    for (int i = row-1; i < column.length(); i++) {
                        String number = line.substring(3*i+2, 3*i+4);
                        number = number.replaceAll("\\s+", "");
                        
                        Integer value = Integer.valueOf(number);
                        String key = "" + aa + column.charAt(i);

                        // System.out.println("key: " + key + " value: " + number);
                        pam100.put(key, value);
                    }
                }
                row++;
            } while (inputFile.hasNextLine());
        } catch (FileNotFoundException fnfe) {
            fnfe.printStackTrace();
        }
    }

    /*  parses command line params.
    */
    public void start(Hashtable<String, String> parameters) {
        //assumption: all the params are given.
        String openGapPenalty = parameters.get("-o");
        if (openGapPenalty == null) {
            openGapPenalty = parameters.get("--ogap");
            if (openGapPenalty == null) {
                System.out.println("open gap param is missing.");
                System.exit(1);
            }
        }
        this.openGapPenalty = Integer.parseInt(openGapPenalty);

        String gapExtensionPenalty = parameters.get("-e");
        if (gapExtensionPenalty == null) {
            gapExtensionPenalty = parameters.get("--egap");
            if (gapExtensionPenalty == null) {
                System.out.println("gap extension param is missing");
                System.exit(1);
            }
        }
        this.gapExtensionPenalty = Integer.parseInt(gapExtensionPenalty);

        String matchScore = parameters.get("-m");
        if (matchScore == null) {
            matchScore = parameters.get("--match");
            if (matchScore == null) {
                System.out.println("match score param is missing");
                System.exit(1);
            }
        }
        this.matchScore = Integer.parseInt(matchScore);

        String nonmatch = parameters.get("-n");
        if (nonmatch == null) {
            nonmatch = parameters.get("--nonmatch");
            if (nonmatch == null) {
                System.out.println("mismatch score param is missing");
                System.exit(1);
            }
        } 
        this.mismatchPenalty = Integer.parseInt(nonmatch);

        String type = parameters.get("-t");
        if (type == null) {
            type = parameters.get("--type");
            if (type == null) {
                System.out.println("type param is missing");
                System.exit(1);
            } 
            if (!type.matches("[NnPp]+")) {
                System.out.println("illegal value entered for type parameter.");
                System.exit(1);
            }
        }
        this.seq_type = type;

        String filename = parameters.get("-f"); 
        if (filename == null) {
            filename = parameters.get("--file");
            if (filename == null) {
                System.out.println("file name param is missing");
                System.exit(1);
            }
        }
        this.parseFile(filename);

        if (this.seq_type.equals("p") ||  this.seq_type.equals("P")) {
        //another parameter specifying location of pam matrix? -p?
            String pamFile = parameters.get("-p");
            if (pamFile == null) {
                pamFile = parameters.get("--pamfile");
                if (pamFile == null) {
                    System.out.println("pam file param is missing");
                    System.exit(1);
                }
            }
            this.pamFile = pamFile;
        }
        // System.out.println("got to end of start method");
    }

    public int[][] createMatrix(String seq1, String seq2) {
        // System.out.println("seq1: " + seq1);
        // System.out.println("seq2: " + seq2);
        return new int[seq1.length()+1][seq2.length()+1];
    }

    //matrix dimensions have already been determined.
    public void initializeMatrix(int[][] matrix, Direction[][] cell_origin) {
        //-5 = opening penalty
        //-1 = extension penalty
        // System.out.println("public void initializeMatrix(int[][] matrix) ");
        for (int i = 0; i < matrix.length; i++) {
            if (i == 0) {
                for (int j = 0; j < matrix[i].length; j++) {
                    matrix[i][j] = 0;
                } 
            } else { //i = 1, 2, etc.
                matrix[i][0] = 0;
            }
        }

        for (int i = 0; i < cell_origin.length; i++) {
            if (i == 0) {
                for (int j = 0; j < cell_origin[i].length; j++) {
                    if (j == 0) {
                        cell_origin[i][j] = Direction.STRT;
                    } else {
                        cell_origin[i][j] = Direction.LEFT;
                    }
                }
            } else {
                cell_origin[i][0] = Direction.TOP;
            }
        }

        // this.printMatrix(matrix, cell_origin, "", "");
     }

    /* The global gap alignment algorithm begins here.
    */
    public void needlemanWunsch(int[][] matrix, Direction[][] cell_origin, String seq1, String seq2) {
        // System.out.println("public void needlemanWunsch(int[][] matrix, Direction[][] cell_origin, String seq1, String seq2)");
        // this.printMatrix(matrix, cell_origin, seq1, seq2);
        //important that both i and j start with 1, 
        //because matrix 1st row and column has already been initialized!!
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix[0].length; j++) {
                matrix[i][j] = this.scoreEntry(matrix, cell_origin, i, j, seq1, seq2);
            }
        }

        System.out.println("Ran the al-gore-rhythm, here is the matrix:");
        this.printMatrix(matrix, cell_origin, seq1, seq2);
    }

    /* Assigns max scores to matrix, assigns direction to cell_origin matrix */
    public int scoreEntry(int[][] matrix, Direction[][] cell_origin, int i, int j, String seq1, String seq2) {
        // System.out.println("public void scoreEntry(int[][] matrix, int i, int j)");
        // this.printMatrix(matrix, cell_origin, seq1, seq2);

        //from the left
        int left = matrix[i][j-1] + this.decideGapPenalty(matrix, cell_origin, i, j, i, j-1);
        int max = left; //max score so far. 
        cell_origin[i][j] = Direction.LEFT;

        //from top therfore decrement i.
        int top = matrix[i-1][j] + this.decideGapPenalty(matrix, cell_origin, i, j, i-1, j);
        
        if (top > max) {
            max = top;
            cell_origin[i][j] = Direction.TOP;
        } else if (top == left) {
          //tie score, what to do? max is still same value, but cell_origin needs to take note
            cell_origin[i][j] = Direction.LT;
        }

        //DG, LD, TD, ALL;
        int diag = matrix[i-1][j-1] + this.matchOrMismatch(i-1, j-1, seq1, seq2);
        if (diag > max) {
            max = diag;
            cell_origin[i][j] = Direction.DG;
        } else if (diag == max && cell_origin[i][j] == Direction.LEFT) {
            cell_origin[i][j] = Direction.LD;
        } else if (diag == max && cell_origin[i][j] == Direction.TOP) {
            cell_origin[i][j] = Direction.TD;
        } else if (diag == max && cell_origin[i][j] == Direction.LT) {
            cell_origin[i][j] = Direction.ALL;
        }

        // System.out.println("left: " + left + " top: " + top + " diag: " + diag + " max: " + max);
        return max;
    }

    /* 
    * decides open / extend gap penalty.
    */
    public int decideGapPenalty(int[][] matrix, Direction[][] cell_origin, int i, int j, int prev_i, int prev_j) {
        if (prev_i < i) { //from top
            if (prev_i != 0 && ( 
                matrix[prev_i][prev_j] - matrix[prev_i-1][prev_j] == gapExtensionPenalty ||
                matrix[prev_i][prev_j] - matrix[prev_i-1][prev_j] == openGapPenalty))
            { //this.cellOrigin(matrix, prev_i, prev_j, seq1, seq2).equals("TOP")
                //the previous cell either started the gap from the same direction or it extended it.
                //prev_i can't be 0, if it were we can't extend the gap from top.
                return gapExtensionPenalty;
            } else { 
                return openGapPenalty;
            }
        } else { //prev_j < j
            //from left
            //analogous logic
            if (prev_j != 0 && (
                matrix[prev_i][prev_j] - matrix[prev_i][prev_j-1] == gapExtensionPenalty ||
                matrix[prev_i][prev_j] - matrix[prev_i][prev_j-1] == openGapPenalty)) 
            { //this.cellOrigin(matrix, prev_i, prev_j, seq1, seq2).equals("LEFT")
                return gapExtensionPenalty;
            } else {
                return openGapPenalty;
            }
        }
    }

    public int matchOrMismatch(int i, int j, String seq1, String seq2) {
        // System.out.println("sequence type: " + seq_type);
        if (seq_type.equals("N") || seq_type.equals("n")) { 
            if (seq1.charAt(i) == seq2.charAt(j)) {
                return matchScore;
            } else {
                return mismatchPenalty;
            }
        } else { //protein comparison. use that pam matrix.
            Integer score = new Integer(0);
            char x = seq1.charAt(i);
            char y = seq2.charAt(j);

            // System.out.println("Comparing: " + x + " and " + y + "...");
            
            score = this.pam100.get("" + x + y);
            if (score == null) {
                score = this.pam100.get("" + y + x);
            }

            // System.out.println("protein score: " + score);
            return score.intValue();
        }
    }

    /* BULK OF THE WORK HERE 
    *  j is the column, i is the row.
    */
    public void getAlignment(int[][] matrix, Direction[][] cell_origin, String key1, String key2, String seq1, String seq2) {
        //first, get the max value of the bottom row of 'matrix'
        int max = 0;
        int max_j_index = 0;
        int max_i_index = 0;
        int i = matrix.length-1;

        List<Integer>max_i_indices = new ArrayList<Integer>();
        List<Integer>max_j_indices = new ArrayList<Integer>();

        for (int j = 0; j < matrix[i].length; j++) {
            if (j == 0) {
                max = matrix[i][j];
                max_j_index = 0;

                max_j_indices.add(new Integer(max_j_index));
                max_i_indices.add(new Integer(max_i_index));
            } else {
                if (matrix[i][j] > max) {
                    max = matrix[i][j];
                    // max_j_index = j;
                    // max_i_index = i;
                    max_j_indices.clear();
                    max_i_indices.clear();

                    // System.out.println("max_j_indices.size(); " + max_j_indices.size());
                    // System.out.println("clearing prev, then adding: " + i + " , " + j);

                    max_j_indices.add(new Integer(j));
                    max_i_indices.add(new Integer(i));
                } else if (matrix[i][j] == max) {
                    //tie. add without clearing.
                    // System.out.println("max_j_indices.size(); " + max_j_indices.size());
                    // System.out.println("adding without clearing: " + i + " , " + j);
                    

                    max_j_indices.add(new Integer(j));
                    max_i_indices.add(new Integer(i));
                }
            }
        }
        //now we have the max along the bottom row.
        //now going along the rightmost column.
        int j = matrix[0].length-1;
        for (int k = 0; k < matrix.length-1; k++) {
            if (matrix[k][j] > max) {
                // System.out.println("along the rightmost column, coordinates of max score: " + i + " , " + j);
                max_j_indices.clear();
                max_i_indices.clear();

                max_j_indices.add(new Integer(j));
                max_i_indices.add(new Integer(k));
                max = matrix[k][j];
            } else if (matrix[k][j] == max) {
                // System.out.println("max_j_indices.size(); " + max_j_indices.size());
                // System.out.println("adding without clearing: " + k + " , " + j);
                
                max_j_indices.add(new Integer(j));
                max_i_indices.add(new Integer(k));
            }
        }

        //assuming one max value for now...
        // System.out.println("Max: " + max + " at j: " + max_j_index + " at i: " + max_i_index);

        //@todo: before traverse, pretty print terminal gap
        List<String> terminalGaps = new ArrayList<String>();

        // System.out.println("max_j_indices.size(); " + max_j_indices.size());
        for (int k = 0; k < max_j_indices.size(); k++) {
            String terminalGap = "";
            i = max_i_indices.get(k).intValue();
            j = max_j_indices.get(k).intValue();
            
            // System.out.println("coordinates of max score: " + i + " , " + j);
            int m = matrix.length - 1;
            int n = matrix[0].length - 1;
            // System.out.println("matrix.length-1: " +  m + " , maxtrix[0].length -1: " + n);
            if (i != matrix.length-1 || j != matrix[0].length -1 ) {
                // System.out.println("adding a terminal gap");
                terminalGaps.add(k, prettyPrintTerminal(i, j, matrix, cell_origin, seq1, seq2));
            } else {
                terminalGaps.add(k, "");
            }

            this.traverse(i, j, max, matrix, cell_origin, key1, key2, seq1, seq2, terminalGaps.get(k));
            this.numSolutions = 0;
        }
    }

    /* STRT, LEFT, TOP, DG, LD, TD, LT, ALL;
    */
    public void traverse(int i, int j, int max, int[][] matrix, Direction[][] cell_origin, String key1, String key2, String seq1, String seq2, String output) { //should it return String
        if (this.numSolutions == 5 && !this.multiseq) {
            System.out.println("Multiple Solutions were found...the top 5 solutions are printed above.");
            System.out.println("This algorithm traverses multiple solutions by doing diagonals (matches)"); 
            System.out.println("first. Therefore the best solution is most likely printed.");
            System.exit(0);
        }
        if (cell_origin[i][j] == Direction.STRT) {
            //base case
            this.numSolutions++;
            // System.out.println("number of solutions: " + this.numSolutions);
            // System.out.println("base case reached. output: " + output);
            if (this.multiseq) {
                String key = key1 + ":" + key2;
                if (this.msa_pairs.get(key) == null) {
                    //updates msa_pairs!
                    System.out.println("adding pair " + key + " into msa hashes");
                    this.msa_pairs.put(key, output);  
                    this.msa_pairs_scores.put(key, new Integer(max));
                    System.out.println("putting into msa_pairs_distances key: " + key + " distance: " + this.distance(output));
                    
                    this.msa_pairs_distances.put(key, this.distance(output));
                } else {
                    // since traverse method favors going diagonal in case of ties, I am inclined to comment out the following
                    //String entry = this.msa_pairs.get(key);
                    // entry += "\n" + output;
                    // this.alignment.put(key, entry);
                }
            } else {
                prettyPrintAlignment(key1, key2, max, output);    
            }
        } else if (cell_origin[i][j] == Direction.DG && (this.numSolutions < 5 || this.multiseq)) { // diagonal
            output += String.valueOf(seq1.charAt(i-1)) + String.valueOf(seq2.charAt(j-1));
            // System.out.println("diag case reached. output: " + output);
            this.traverse(i-1, j-1, max, matrix, cell_origin, key1, key2, seq1, seq2, output);  
        
        } else if (cell_origin[i][j] == Direction.LEFT && (this.numSolutions < 5 || this.multiseq)) {
            // System.out.println("left case reached.");    
            output += "-" + String.valueOf(seq2.charAt(j-1));
            this.traverse(i, j-1, max, matrix, cell_origin, key1, key2, seq1, seq2, output);
        
        } else if (cell_origin[i][j] == Direction.TOP && (this.numSolutions < 5 || this.multiseq)) {
            // System.out.println("top case reached.");    
            output += String.valueOf(seq1.charAt(i-1)) + "-";
            this.traverse(i-1, j, max, matrix, cell_origin, key1, key2, seq1, seq2, output);
        
        } else if (cell_origin[i][j] == Direction.LD && (this.numSolutions < 5 || this.multiseq)) {
            // System.out.println("tie between left, diag");
            String temp = output;
            output += String.valueOf(seq1.charAt(i-1)) + String.valueOf(seq2.charAt(j-1));
            this.traverse(i-1, j-1, max, matrix, cell_origin, key1, key2, seq1, seq2, output);  //traversing diag
            
            // System.out.println("done traversing diag, now going left.");
            temp += "-" + String.valueOf(seq2.charAt(j-1));
            this.traverse(i , j-1, max, matrix, cell_origin, key1, key2, seq1, seq2, temp);

        } else if (cell_origin[i][j] == Direction.TD && (this.numSolutions < 5 || this.multiseq)) {
            // System.out.println("tie between top, diag");
            String temp = output;
            output += String.valueOf(seq1.charAt(i-1)) + String.valueOf(seq2.charAt(j-1));
            this.traverse(i-1, j-1, max, matrix, cell_origin, key1, key2, seq1, seq2, output);  //traversing diag
            
            // System.out.println("done traversing diag, now going top.");
            temp += String.valueOf(seq1.charAt(i-1)) + "-";
            this.traverse(i-1, j, max, matrix, cell_origin, key1, key2, seq1, seq2, temp);
        
        } else if (cell_origin[i][j] == Direction.LT && (this.numSolutions < 5 || this.multiseq)) {
            // System.out.println("tie between top, left. go left first");
            String temp = output;
            
            output += "-" + String.valueOf(seq2.charAt(j-1));
            this.traverse(i, j-1, max, matrix, cell_origin, key1, key2, seq1, seq2, output);

            // System.out.println("done going left, now going top");
            temp += String.valueOf(seq1.charAt(i-1)) + "-";
            this.traverse(i-1, j, max, matrix, cell_origin, key1, key2, seq1, seq2, temp);

        } else if (cell_origin[i][j] == Direction.ALL && (this.numSolutions < 5 || this.multiseq)) {
            // System.out.println("3 way tie");
            String temp = output;
            String temp2 = output;

            //diag
            output += String.valueOf(seq1.charAt(i-1)) + String.valueOf(seq2.charAt(j-1));
            this.traverse(i-1, j-1, max, matrix, cell_origin, key1, key2, seq1, seq2, output);  
            
            //left 
            temp += "-" + String.valueOf(seq2.charAt(j-1));
            this.traverse(i, j-1, max, matrix, cell_origin, key1, key2, seq1, seq2, temp);

            //top
            temp2 += String.valueOf(seq1.charAt(i-1)) + "-";
            this.traverse(i-1, j, max, matrix, cell_origin, key1, key2, seq1, seq2, temp2);
        }

        // System.out.println("got to end of this traverse method");
    }

    public String consensus(String input) {
        if (input.length() % 2 != 0) {
            System.out.println("alignment string should have an even length!");
            System.exit(1);
        }

        String reverse = "";
        for (int i = input.length()-1; i>=0; i--) {
            reverse += String.valueOf(input.charAt(i));
        }
        
        String top = "";
        String bottom = "";
        for (int i = 0; i < input.length(); i++) {
            if (i % 2 == 0) {
                bottom += String.valueOf(reverse.charAt(i));
            } else {
                top += String.valueOf(reverse.charAt(i));
            }
        }

        String toReturn = "";

        for (int i = 0; i< top.length(); i++) {
            if (top.charAt(i) == bottom.charAt(i)) {
                toReturn += "" + top.charAt(i);
            } else { //more on this later(?)
                toReturn += "-";
            }
        }

        return toReturn;
    }

    public double distance(String key1, String key2) {
        String alignment = this.msa_pairs.get(key1 + ":" + key2);

        if (alignment == null) {
            alignment = this.msa_pairs.get(key2 + ":" + key1);
        }

        // String seq2 = this.msa_pairs.get(key2);
        return this.distance(alignment);
    }

    //returns pairwise distance of the alignment. 
    //kimura distances returned NaN in some cases due to negative logs, which were not useful.
    public Double distance(String alignment) {
        int distance = 0;
        for (int i = 0; i < alignment.length()-1; i++) {
            if (alignment.charAt(i) != alignment.charAt(i+1)) {
                distance++;
            }
        }
        return (double) distance;

        // int transitions = 0;
        // int transversions = 0;
        // int length = 0;

        // if (alignment.length() % 2 != 0) {
        //     System.out.println("error, alignment string is of odd length");
        //     System.exit(1);
        // } else {
        //     length = alignment.length() / 2;
        //     System.out.println("computing kimura distance, length: " + length);
        //     for (int i = 0; i < alignment.length()-1; i++) {
        //         if (alignment.charAt(i) == 'A' && alignment.charAt(i+1) == 'G' || 
        //             alignment.charAt(i) == 'G' && alignment.charAt(i+1) == 'A' ||
        //             alignment.charAt(i) == 'C' && alignment.charAt(i+1) == 'T' ||
        //             alignment.charAt(i) == 'T' && alignment.charAt(i+1) == 'C' ||
        //             alignment.charAt(i) == 'C' && alignment.charAt(i+1) == 'U' ||
        //             alignment.charAt(i) == 'U' && alignment.charAt(i+1) == 'C' ) {
        //             transitions++;
        //         } else if (
        //             alignment.charAt(i) == 'A' && alignment.charAt(i+1) == 'T' ||
        //             alignment.charAt(i) == 'T' && alignment.charAt(i+1) == 'A' ||
        //             alignment.charAt(i) == 'A' && alignment.charAt(i+1) == 'U' || 
        //             alignment.charAt(i) == 'U' && alignment.charAt(i+1) == 'A' || 

        //             alignment.charAt(i) == 'G' && alignment.charAt(i+1) == 'T' || 
        //             alignment.charAt(i) == 'G' && alignment.charAt(i+1) == 'U' || 
        //             alignment.charAt(i) == 'U' && alignment.charAt(i+1) == 'G' ||  
        //             alignment.charAt(i) == 'T' && alignment.charAt(i+1) == 'G' ||

        //             alignment.charAt(i) == 'A' && alignment.charAt(i+1) == 'C' || 
        //             alignment.charAt(i) == 'C' && alignment.charAt(i+1) == 'A' ||

        //             alignment.charAt(i) == 'G' && alignment.charAt(i+1) == 'C' || 
        //             alignment.charAt(i) == 'C' && alignment.charAt(i+1) == 'G' ) {
        //             transversions++;
        //         }
        //     }
        // }

        // double p = (double) transitions / (double) length;
        // double q = (double) transversions / (double) length;

        // System.out.println("p: " + p);
        // System.out.println("q: " + q);
        // System.out.println("(1/(1 - 2*p -q)) : " + 1/(1-2*p-q));
        // System.out.println("0.25*Math.log(1/(1 - 2*q)) : " + 0.25*Math.log(1/(1 - 2*q)));
        // return 0.5*Math.log(1/(1-2*p-q)) + 0.25*Math.log(1/(1 - 2*q));
    }

    //
    public double computeAvgOutgroupDistance(String outGroup) {
        String[] keys = this.msa_pairs_distances.keySet().toArray(new String[this.msa_pairs_distances.size()]);
        int sum = 0;
        // System.out.println("computing avg distance, sum: " + sum );
        for (int i = 0; i < keys.length; i++) {
            if (keys[i].contains(outGroup)) {
                sum+=this.msa_pairs_distances.get(keys[i]);
                // System.out.println("computing avg distance, i: " + i + " sum: " + sum ); 
            }
        }
        return (double) sum / keys.length;
    }

    public double transformDistance(String speciesA, String speciesB, String outGroup) {
        System.out.println("in transformDistance");
        System.out.println("getting msa_pairs_distances, key: " + speciesA + ":" + speciesB);
        Double ab_distance = this.msa_pairs_distances.get(speciesA + ":" + speciesB);
        if (ab_distance == null) {
            System.out.println("ab_distance was null, reversing the key");
            System.out.println("getting msa_pairs_distances, key: " + speciesB + ":" + speciesA);
            ab_distance = this.msa_pairs_distances.get(speciesB + ":" + speciesA);   
            if (ab_distance == null) {
                System.out.println("paired distance is still null...crashing");
                System.exit(1);
            }
        }  
        double ab = (double) ab_distance.intValue();

        System.out.println("getting msa_pairs_distances, key: " + speciesA + ":" + outGroup);
        Double ao_distance = this.msa_pairs_distances.get(speciesA + ":" + outGroup); 
        if (ao_distance == null) {
            System.out.println("ao_distance was null, reversing the key");
            System.out.println("getting msa_pairs_distances, key: " + outGroup + ":" + speciesA);
            ao_distance = this.msa_pairs_distances.get(outGroup + ":" + speciesA);
            if (ao_distance == null) {
                System.out.println("ao_distance is still null...crashing");
                System.exit(1);
            }
        }
        double ao = (double) ao_distance.intValue();
        
        System.out.println("getting msa_pairs_distances, key: " + speciesB + ":" + outGroup);    
        Double bo_distance = this.msa_pairs_distances.get(speciesB + ":" + outGroup); 
        if (bo_distance == null) {
            System.out.println("bo_distance was null, reversing the key");
            System.out.println("getting msa_pairs_distances, key: " + outGroup + ":" + speciesB);
            bo_distance = this.msa_pairs_distances.get(outGroup + ":" + speciesB);
            if (bo_distance == null) {
                System.out.println("bo_distance is still null...crashing");
                System.exit(1);
            }
        }
        double bo = (double) bo_distance.intValue();
    
        return (ab - ao - bo) / 2.0 + this.avgOutgroupDistance;
    }

    public void prettyPrintAlignment(String key1, String key2, int max, String alignment) {
        String reverse = "";
        for (int i = alignment.length()-1; i>=0; i--) {
            reverse += String.valueOf(alignment.charAt(i));
        }
        String top = "";
        String bottom = "";
        for (int i = 0; i < alignment.length(); i++) {
            if (i % 2 == 0) {
                bottom += String.valueOf(reverse.charAt(i));
            } else {
                top += String.valueOf(reverse.charAt(i));
            }
        }
         //todo: print things out better-- 
        //24 chars to a line, and print alignment # (eg 1st possible alignment etc, total number of alignments)
        System.out.println("Alignment Score: " + max);
        // System.out.println("saved to alignment scores: " + this.alignment_scores.get(key1 + ":" + key2));
        int bottom_index = 0;
        for (int i = 0; i < top.length(); i++) {             
            System.out.print(top.charAt(i));    
            if (((i % 60 == 0) && i > 0) || (i == top.length()-1)) {
                System.out.print("\n");
                for (int j = bottom_index; j<=i; j++) {
                    System.out.print("|");
                    if (j == i) {
                        System.out.print("\n");
                    }
                }
                for (int k = bottom_index; k<=i; k++) {
                    System.out.print(bottom.charAt(k));
                    bottom_index++;
                }
                System.out.print("\n\n");
            }
        }

        // System.out.println(top);
        // System.out.println(bottom);

        // System.out.println("saved to alignment hash: " + this.alignment.get(key1 + ":" + key2));
    }

    public String prettyPrintTerminal(int i, int j, int[][] matrix, Direction[][] cell_origin, String seq1, String seq2) {
        String terminal = "";
        if (i < matrix.length-1) { //not the bottom row
            // System.out.println("not the bottom row");
            for (int k = i; k < matrix.length-1; k++) {
                terminal += "" + "-" + seq1.charAt(k);
                // System.out.println("terminal: " + terminal);
            }
            //example: "-A-A-A"
        } else { //not the rightmost column
            // System.out.println("not the rightmost column");
            for (int l = j; l < matrix[0].length-1; l++) {
                terminal += seq2.charAt(l) + "-";
            }
        }

        String reverse = "";
        for (int k = terminal.length()-1; k>=0; k--) {
            reverse += String.valueOf(terminal.charAt(k));
        }
        
        return reverse;
    }

    public void printMatrix(int[][] matrix, Direction[][] cell_origin, String seq1, String seq2) {
        System.out.println("here's the alignment matrix:");
        System.out.println("seq1:" + seq1);
        System.out.println("seq2:" + seq2);
        
        System.out.print("\t\t");
        for (int j = 0; j < seq2.length(); j++) {
            System.out.print(seq2.charAt(j) + "\t");
        }
        System.out.println();

        for (int i = 0; i < matrix.length; i++) {
            if (i == 0) {
                System.out.print("\t");
            }
            if (i >= 1) {
                System.out.print(seq1.charAt(i-1) + "\t");
            }
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.print(matrix[i][j] + "\t");
            }   
            System.out.println();
        }

        System.out.println("here's the cell origin matrix"); 
        System.out.print("\t\t");
        for (int j = 0; j < seq2.length(); j++) {
            System.out.print(seq2.charAt(j) + "\t");
        }
        System.out.println();

        for (int i = 0; i < cell_origin.length; i++) {
            if (i == 0) {
                System.out.print("\t");
            }
            if (i >= 1) {
                System.out.print(seq1.charAt(i-1) + "\t");
            }
            for (int j = 0; j < cell_origin[i].length; j++) {
                System.out.print(cell_origin[i][j] + "\t");
            }   
            System.out.println();
        }
    }
}