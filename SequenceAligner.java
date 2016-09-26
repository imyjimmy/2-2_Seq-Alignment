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
    
    // public Hashtable<String, String> parameters = new Hashtable<String, String>();

    //some default values
    public int openGapPenalty = -5;
    public int gapExtensionPenalty = -1;
    public int matchScore = 1;
    public int mismatchPenalty = -1;

    public String seq_type;

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
        if (args.length > 0 && (args[0].equals("-i") || args[0].equals("--interactive"))) {
            //--interactive or -i
            seq_aligner.startInteractive();
            // parameters.put(args[0], "True");
        } else if (args.length > 0 && (args[0].equals("-h") || args[0].equals("--help"))) {
            System.out.println("params: --interactive or -i");
            System.out.println("-h or --help to print this screen.");
            System.out.println("--type or -t (n|N or p|P) for nucleotide or protein match");
            System.out.println("--file or -f [filepath]");
            System.out.println("-m or --match [score]");
            System.out.println("-n or --notmatch [score]");
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

        if (keys.length == 2) {
            String seq1 = s.get(keys[0]);
            String seq2 = s.get(keys[1]);

            int[][] matrix = seq_aligner.createMatrix(seq1, seq2);
            Direction[][] cell_origin_matrix = new Direction[seq1.length()+1][seq2.length()+1];

            seq_aligner.printMatrix(matrix, cell_origin_matrix, seq1, seq2);
            seq_aligner.initializeMatrix(matrix, cell_origin_matrix);

            //Begin the algorithm.
            seq_aligner.needlemanWunsch(matrix, cell_origin_matrix, seq1, seq2);

            //return the alignment.
            seq_aligner.getAlignment(matrix, cell_origin_matrix);

        } else { //multiple pair alignment. todo.
            //nothing for now.
        }
     }

     /*   public int gapExtensionPenalty = -1;
    public int matchScore = 1;
    public int mismatchPenalty = -1;
    */
    public void start() {
        //assumption: all the params are given.
        this.openGapPenalty = parameters.get("-o");
        if (this.openGapPenalty == null) {
            System.out.println("open gap param is missing.");
            this.openGapPenalty = parameters.get("-ogap");
        }

        this.gapExtensionPenalty = parameters.get("-e");
    }

     /*
     * The interactive interface to this program.
     * @todo: ask for scoring params.
     */
    public void startInteractive() {
        System.out.println("Welcome to Sequence Alignment. \nWe will help you align two sequences. \nLet's begin. \nProtein (p|P) or Nucleic Acid (n|N)?");
        //store response.
        Scanner scan = new Scanner(System.in);
        seq_type = scan.nextLine();

        System.out.println("Please give a complete file name. \nInclude the relative path to file \nfrom this directory or \nan absolute file path:");

        String filename = scan.nextLine();
        File file = new File(filename);
        
        try { 
            Scanner inputFile = new Scanner(file); 
            String key = "";
            String val = "";
            do {
                String line = inputFile.nextLine();
                System.out.println(line);
                if (line.startsWith(">")) {
                    key = line.substring(1);
                    System.out.println("substring: " + key);
                } else {
                    val = line;
                }

                if (!"".equals(val) && !"".equals(key)) { //avoids null pointer exceptions.
                    sequences.put(key, val);
                    key = "";
                    val = "";
                }
            } while (inputFile.hasNextLine());

            /* 
            * Assumption: the sequence is now loaded.
            */

        } catch (FileNotFoundException fnfe) {
            fnfe.printStackTrace();
        }

        System.out.println("configure the opening gap penalty: (-5)");
        System.out.println("configure the gap extension penalty: (-1)");
     }

    public int[][] createMatrix(String seq1, String seq2) {
        System.out.println("seq1: " + seq1);
        System.out.println("seq2: " + seq2);
        return new int[seq1.length()+1][seq2.length()+1];
    }

    public void initializeMatrix(int[][] matrix, Direction[][] cell_origin) {
        //-5 = opening penalty
        //-1 = extension penalty
        System.out.println("public void initializeMatrix(int[][] matrix) ");
        for (int i = 0; i < matrix.length; i++) {
            if (i == 0) {
                for (int j = 0; j < matrix[i].length; j++) {
                    if (j == 0) {
                        matrix[i][j] = 0;
                    } else {
                        matrix[i][j] = -5+(-1*(j-1));
                    }
                } 
            } else { //i = 1, 2, etc.
                matrix[i][0] = -5+(-1*(i-1));
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
        System.out.println("public void needlemanWunsch(int[][] matrix, Direction[][] cell_origin, String seq1, String seq2)");
        this.printMatrix(matrix, cell_origin, seq1, seq2);
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

    /* Assigns max scores to matrix */
    public int scoreEntry(int[][] matrix, Direction[][] cell_origin, int i, int j, String seq1, String seq2) {
        System.out.println("public void scoreEntry(int[][] matrix, int i, int j)");
        this.printMatrix(matrix, cell_origin, seq1, seq2);

        //from the left
        int left = matrix[i][j-1] + this.decideGapPenalty(matrix, i, j, i, j-1);
        int max = left; //max score so far. 
        cell_origin[i][j] = Direction.LEFT;

        //from top therfore decrement i.
        int top = matrix[i-1][j] + this.decideGapPenalty(matrix, i, j, i-1, j);
        
        if (top > max) {
            max = top;
            cell_origin[i][j] = Direction.TOP;
        } else if (top == left) {
          //tie score, what to do? max is still same value, but cell_origin needs to take note
            cell_origin[i][j] = Direction.LT;
        }

        //DG, LD, TD, LT, ALL;
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

        System.out.println("left: " + left + " top: " + top + " diag: " + diag + " max: " + max);
        return max;
    }

    /* 
    * LOTS OF ASSUMPTIONS
    * when u assuuuume
    * THIS WILL BE FIXED TOO
    */
    public int decideGapPenalty(int[][] matrix, int i, int j, int prev_i, int prev_j) {
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
            return 0;
        }
    }

    //is this method even worth it?
    public Direction cellOrigin(Direction[][] cell_origin, int i, int j) {
        return cell_origin[i][j];
    }

    /* BULK OF THE WORK HERE 
    *  
    */
    public void getAlignment(int[][] matrix, Direction[][] cell_origin) {
        //first, get the max value of the bottom row of 'matrix'
        int max = 0;
        int max_j_index = 0;
        int i = matrix.length-1;
        for (int j = 0; j < matrix[i].length; j++) {
            if (j == 0) {
                max = matrix[i][j];
                max_j_index = 0;
            } else {
                if (matrix[i][j] > max) {
                    max = matrix[i][j];
                    max_j_index = j;
                }
            }
        }
        //now we have the max.
        System.out.println("Max: " + max + " at: " + max_j_index);

        // Stack<Direction[]> path = new Stack<Direction[]>();
        // System.out.println(cell_origin[i][max_j_index]);
        // path.push()
        //@todo: depth first traversal.
        //to construct
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
        for (int i = 0; i < cell_origin.length; i++) {
            for (int j = 0; j < cell_origin[i].length; j++) {
                System.out.print(cell_origin[i][j] + "\t");
            }   
            System.out.println();
        }
    }
}