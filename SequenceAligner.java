/* Data Structures */
import java.util.ArrayList;
import java.util.List;
import java.util.Hashtable;
import java.util.Set;

/* For file input and outputs */
import java.util.Scanner;
import java.io.File;
import java.io.FileNotFoundException;

public class SequenceAligner{
    public Hashtable<String, String> sequences = new Hashtable<String, String>();
    public int openGapPenalty = -5;
    public int gapExtensionPenalty = -1;
    public int matchScore = 1;
    public int mismatchPenalty = -1;

    /*
     Todo: as the user for scoring parameters. As noted in the assignment.
     Particularly in regards to the affine gap penalty: A+(B⋅L)
    */
    public static void main(String []args){
        SequenceAligner seq_aligner = new SequenceAligner();
        seq_aligner.start();

        Hashtable<String, String> s = seq_aligner.sequences;
        String[] keys = s.keySet().toArray(new String[s.size()]);

        if (keys.length == 2) {
            String seq1 = s.get(keys[0]);
            String seq2 = s.get(keys[1]);
            int[][] matrix = seq_aligner.createMatrix(seq1, seq2);

            seq_aligner.printMatrix(matrix);
            seq_aligner.initializeMatrix(matrix);

            //Begin the algorithm.
            seq_aligner.needlemanWunsch(matrix, seq1, seq2);

        } else { //multiple pair alignment. todo.
            //nothing for now.
        }
     }

     /*
     * The interface to this program.
     */
    public void start() {
        System.out.println("Welcome to Sequence Alignment. \nWe will help you align two sequences. \nLet's begin. \nProtein (p|P) or Nucleic Acid (n|N)?");
        //store response.
        Scanner scan = new Scanner(System.in);
        String seq_type = scan.nextLine();

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
        return new int[seq1.length()+1][seq2.length()+1];
    }

    public void initializeMatrix(int[][] matrix) {
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
            }
            else { //i = 1, 2, etc.
                matrix[i][0] = -5+(-1*(i-1));
            }
        }

        this.printMatrix(matrix);
     }

    /* The global gap alignment algorithm begins here.
    */
    public void needlemanWunsch(int[][] matrix, String seq1, String seq2) {
        System.out.println("public void needlemanWunsch(int[][] matrix, String seq1, String seq2)");
        this.printMatrix(matrix);
        //important that both i and j start with 1, 
        //because matrix 1st row and column has already been initialized!!
        for (int i = 1; i < matrix.length; i++) {
            for (int j = 1; j < matrix.length; j++) {
                matrix[i][j] = this.findMaxScore(matrix, i, j, seq1, seq2);
            }
        }

        System.out.println("Ran the al-gore-rhythm, here is the matrix:");
        this.printMatrix(matrix);
    }

    public int findMaxScore(int[][] matrix, int i, int j, String seq1, String seq2) {
        System.out.println("public void findMaxScore(int[][] matrix, int i, int j)");
        this.printMatrix(matrix);

        //from the left
        int left = matrix[i][j-1] + this.decideGapPenalty(matrix, i, j, i, j-1);
        int max = left; //max score so far.

        //from top therfore decrement i.
        int top = matrix[i-1][j] + this.decideGapPenalty(matrix, i, j, i-1, j);
        
        if (top > max) {
            max = top;
        } else if (top == left) {
          //tie score, what to do? @todo
        }

        int diag = matrix[i-1][j-1] + this.matchOrMismatch(i-1, j-1, seq1, seq2);
        if (diag > max) {
            max = diag;
        }

        return max;
    }

    /* 
    * LOTS OF ASSUMPTIONS
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
/*
    THIS NEEDS WORK RIGHT NOW.
    4   10
    4   5

*/
    // public String cellOrigin(int[][] matrix, int i, int j, String seq1, String seq2) {
    //     String toReturn;
    //     if () {
    //         return "LEFT";
    //     } else if () {
    //         return "TOP";
    //     } else if (seq1.charAt(i) == seq2.charAt(j) && matrix[i][j] - matrix[i-1][j-1] == matchScore) {
    //         //assume the match score will be awarded from the diagonal when there's a match. 
    //         //(highest possible value)
    //         toReturn = "DIAG";
    //     } else if (seq1.charAt(i) != seq2.charAt(j) && matrix[i][j] - matrix[i-1][j-1] == mismatchPenalty) {
    //         //assume that mismatch penalty will be higher priority than gaps. 
    //         toReturn = "DIAG";
    //     }
    // }

    public int matchOrMismatch(int i, int j, String seq1, String seq2) {
        if (seq1.charAt(i) == seq2.charAt(j)) {
            return matchScore;
        } else {
            return mismatchPenalty;
        }
    }

    public void printMatrix(int[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.print(matrix[i][j] + "\t");
            }   
            System.out.println();
        }
     }
}