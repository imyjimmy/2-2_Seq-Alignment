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

    //some default values
    public int openGapPenalty = -5;
    public int gapExtensionPenalty = -1;
    public int matchScore = 1;
    public int mismatchPenalty = -1;

    public String seq_type;

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
            // Integer i = seq_aligner.pam100.get("IV");
            // System.out.println("getting match IV: " + i.toString());
        }
        
        if (keys.length == 2) {
            // System.out.println("inside keys.length == 2");
            String seq1 = s.get(keys[1]);
            String seq2 = s.get(keys[0]);

            int[][] matrix = seq_aligner.createMatrix(seq1, seq2);
            Direction[][] cell_origin_matrix = new Direction[seq1.length()+1][seq2.length()+1];

            // seq_aligner.printMatrix(matrix, cell_origin_matrix, seq1, seq2);
            seq_aligner.initializeMatrix(matrix, cell_origin_matrix);

            //Begin the algorithm.
            seq_aligner.needlemanWunsch(matrix, cell_origin_matrix, seq1, seq2);

            //return the alignment.
            seq_aligner.getAlignment(matrix, cell_origin_matrix, seq1, seq2);
        } else { //multiple pair alignment. todo.
            //nothing for now.
            System.out.println("starting multi pair alignment");

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
                    if (!"".equals(val) && !"".equals(key)) { //avoids null pointer exceptions.
                        sequences.put(key, val);
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
        File pam = new File("100pam.txt");
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

    /*  
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

        //another parameter specifying location of pam matrix? -p?

        // System.out.println("got to end of start method");
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
        
        this.parseFile(filename);

        System.out.println("configure the opening gap penalty: (-5)");
        System.out.println("configure the gap extension penalty: (-1)");
     }

    public int[][] createMatrix(String seq1, String seq2) {
        System.out.println("seq1: " + seq1);
        System.out.println("seq2: " + seq2);
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

        // System.out.println("Ran the al-gore-rhythm, here is the matrix:");
        // this.printMatrix(matrix, cell_origin, seq1, seq2);
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
    * @todo: LOTS OF ASSUMPTIONS
    * when u assuuuume
    * THIS WILL BE FIXED TOO--
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

    //is this method even worth it?
    public Direction cellOrigin(Direction[][] cell_origin, int i, int j) {
        return cell_origin[i][j];
    }

    /* BULK OF THE WORK HERE 
    *  j is the column, i is the row.
    */
    public void getAlignment(int[][] matrix, Direction[][] cell_origin, String seq1, String seq2) {
        //first, get the max value of the bottom row of 'matrix'
        int max = 0;
        int max_j_index = 0;
        int max_i_index = 0;
        int i = matrix.length-1;

        for (int j = 0; j < matrix[i].length; j++) {
            if (j == 0) {
                max = matrix[i][j];
                max_j_index = 0;
            } else {
                if (matrix[i][j] > max) {
                    max = matrix[i][j];
                    max_j_index = j;
                    max_i_index = i;
                } else if (matrix[i][j] == max) {
                    //we have to worry about tiebreakers..?
                }
            }
        }
        //now we have the max along the bottom row.
        //now going along the rightmost column.
        int j = matrix[0].length-1;
        for (int k = 0; k < matrix.length-1; k++) {
            if (matrix[k][j] > max) {
                max_j_index = j;
                max_i_index = k;
                max = matrix[k][j];
            } else if (matrix[k][j] == max) {
                //tiebreakers.
            }
        }

        //assuming one max value for now...
        // System.out.println("Max: " + max + " at j: " + max_j_index + " at i: " + max_i_index);

        //@todo: before traverse, pretty print terminal gap
        String terminalGap = "";
        if (max_i_index != matrix.length-1 || max_j_index != matrix[0].length-1) {
            terminalGap = prettyPrintTerminal(max_i_index, max_j_index, matrix, cell_origin, seq1, seq2);
        }

        this.traverse(max_i_index, max_j_index, matrix, cell_origin, seq1, seq2, terminalGap);
    }

    /* STRT, LEFT, TOP, DG, LD, TD, LT, ALL;
    */
    public void traverse(int i, int j, int[][] matrix, Direction[][] cell_origin, String seq1, String seq2, String output) { //should it return String
        if (cell_origin[i][j] == Direction.STRT) {
            //base case
            // System.out.println("base case reached. output: " + output);
            //@todo: prettyPrint needs to handle terminal gaps. but it cant be in the base case
            prettyPrintAlignment(output);
        } else if (cell_origin[i][j] == Direction.DG) { // diagonal
            output += String.valueOf(seq1.charAt(i-1)) + String.valueOf(seq2.charAt(j-1));
            // System.out.println("diag case reached. output: " + output);
            this.traverse(i-1, j-1, matrix, cell_origin, seq1, seq2, output);  
        
        } else if (cell_origin[i][j] == Direction.LEFT) {
            // System.out.println("left case reached.");    
            output += "-" + String.valueOf(seq2.charAt(j-1));
            this.traverse(i, j-1, matrix, cell_origin, seq1, seq2, output);
        
        } else if (cell_origin[i][j] == Direction.TOP) {
            // System.out.println("top case reached.");    
            output += String.valueOf(seq1.charAt(i-1)) + "-";
            this.traverse(i-1, j, matrix, cell_origin, seq1, seq2, output);
        
        } else if (cell_origin[i][j] == Direction.LD) {
            // System.out.println("tie between left, diag");
            String temp = output;
            output += String.valueOf(seq1.charAt(i-1)) + String.valueOf(seq2.charAt(j-1));
            this.traverse(i-1, j-1, matrix, cell_origin, seq1, seq2, output);  //traversing diag
            
            // System.out.println("done traversing diag, now going left.");
            temp += "-" + String.valueOf(seq2.charAt(j-1));
            this.traverse(i , j-1, matrix, cell_origin, seq1, seq2, temp);

        } else if (cell_origin[i][j] == Direction.TD) {
            // System.out.println("tie between top, diag");
            String temp = output;
            output += String.valueOf(seq1.charAt(i-1)) + String.valueOf(seq2.charAt(j-1));
            this.traverse(i-1, j-1, matrix, cell_origin, seq1, seq2, output);  //traversing diag
            
            // System.out.println("done traversing diag, now going top.");
            temp += String.valueOf(seq1.charAt(i-1)) + "-";
            this.traverse(i-1, j, matrix, cell_origin, seq1, seq2, temp);
        
        } else if (cell_origin[i][j] == Direction.LT) {
            // System.out.println("tie between top, left. go left first");
            String temp = output;
            
            output += "-" + String.valueOf(seq2.charAt(j-1));
            this.traverse(i, j-1, matrix, cell_origin, seq1, seq2, output);

            // System.out.println("done going left, now going top");
            temp += String.valueOf(seq1.charAt(i-1)) + "-";
            this.traverse(i-1, j, matrix, cell_origin, seq1, seq2, temp);

        } else if (cell_origin[i][j] == Direction.ALL) {
            // System.out.println("3 way tie");
            String temp = output;
            String temp2 = output;

            //diag
            output += String.valueOf(seq1.charAt(i-1)) + String.valueOf(seq2.charAt(j-1));
            this.traverse(i-1, j-1, matrix, cell_origin, seq1, seq2, output);  
            
            //left 
            temp += "-" + String.valueOf(seq2.charAt(j-1));
            this.traverse(i, j-1, matrix, cell_origin, seq1, seq2, temp);

            //top
            temp2 += String.valueOf(seq1.charAt(i-1)) + "-";
            this.traverse(i-1, j, matrix, cell_origin, seq1, seq2, temp2);
        }

        // System.out.println("got to end of this traverse method");
    }

    public void prettyPrintAlignment(String alignment) {
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
        System.out.println(top);
        System.out.println(bottom);
    }

    public String prettyPrintTerminal(int i, int j, int[][] matrix, Direction[][] cell_origin, String seq1, String seq2) {
        String toReturn = "";
        if (i < matrix.length-1) { //not the bottom row
            for (int k = i; k < matrix.length-1; k++) {
                toReturn += "" + seq1.charAt(i-1) + "-";
            }
            //example: "-A-A-A"
        } else { //not the rightmost column
            for (int l = j; l < matrix[0].length-1; l++) {
                toReturn += "-" + seq2.charAt(l-1);
            }
        }
        return toReturn;
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