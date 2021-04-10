import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;


public class Celilmete_Caglasen {
    public static String bases="ATGC";
    public static File file;
    public static String[] gens = new String[10];
    public static String[] motifs = new String[10];
    public static double[][] profile = new double[4][10];
    public static double[][] resMatrix = new double[10][491];

    public static void main(String args[]){
        int k = 10;
        createFile();
//        writeRandomBasesToFile();
        readGens();
        randomizedMotifSearch(k);


    }

    public static void randomizedMotifSearch(int k) {
        getRandomMotif(k);
        int lastScore = 0;
        int tempScore = 0;
        String[] lastMotif= new String[10];
        int counter = 0;
        while (true) {
            while (true) {
                int score = getMotifs(k);
                if (tempScore == score)
                    break;
                else tempScore = score;
                counter++;
            }
            getRandomMotif(k);
            if (lastScore == tempScore) {
                counter += 1;
                break;
            }else if (lastScore > tempScore){
                lastScore = tempScore;
                counter = 0;
            }
        }

        for (int i = 0; i < 10; i++) {
            System.out.println(motifs[i]);
        }
        System.out.println("with total score of " + lastScore);

    }

    /********************************
     * this function reads the gens */
    public static void readGens(){
        try {
            Scanner reader = new Scanner(file);
            int line = 0;
            while (reader.hasNextLine()) {
                gens[line] = reader.nextLine();
                line++;
            }
            reader.close();
        }catch (Exception e){
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    /*************************************
    * this function calculates resMatrix */
    public static void getResMatrix() {
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j <491 ; j++) {
                double temp = 1;
                String tempString = gens[i].substring(j, j + 10);
                for (int k = 0; k < 10; k++) {
                    switch (tempString.charAt(k)) {
                        case 'A' -> temp *= profile[0][k];
                        case 'C' -> temp *= profile[1][k];
                        case 'G' -> temp *= profile[2][k];
                        case 'T' -> temp *= profile[3][k];
                    }
                }
                resMatrix[i][j] = temp;
            }
        }
    }

    /***********************************************
     * This function calculates the profile matrix */
    public static void getProfile() {
        profile = new double[4][10];
        for (int i = 0; i < 10; i++) {
            for (int j = 0; j < 10; j++) {
                switch (motifs[j].charAt(i)) {
                    case 'A' -> profile[0][i] += 0.1;
                    case 'C' -> profile[1][i] += 0.1;
                    case 'G' -> profile[2][i] += 0.1;
                    case 'T' -> profile[3][i] += 0.1;
                }

            }
        }

    }

    /***************************************************************
    * This function calculates profile matrix and gets next motifs */
    public static int getMotifs(int k) {
        getProfile();
        getResMatrix();
        int score = 0;
        double debug = 0;
        for (int i = 0; i < 10; i++) {
            int lastIndex = 0;
            double lastVal = 0;
            for (int j = 0; j < 491; j++) {
                if (resMatrix[i][j] > lastVal) {
                    lastVal = resMatrix[i][j];
                    lastIndex = j;
                }
            }
            debug += lastVal;
            motifs[i] = gens[i].substring(lastIndex, lastIndex + k);
        }
//        System.out.println(debug);
        debug = 0;
        // calculation of the score
        ArrayList<Integer> nums = new ArrayList<>();
        for (int i = 0; i < 10; i++) {
            int numa = 0, numc = 0, numg = 0, numt = 0;
            for (int j = 0; j < 10; j++) {
                switch (motifs[i].charAt(j)) {
                    case 'A' -> numa++;
                    case 'C' -> numc++;
                    case 'G' -> numg++;
                    case 'T' -> numt++;
                }
            }
            nums.add(numa);nums.add(numc);nums.add(numg);nums.add(numt);
            Collections.sort(nums);
            score += 10 - nums.get(3);
            nums.clear();
        }

        return score;
    }

    /************************************************************************
    * This function gets motifs starting from random positions in each line */
    public static void getRandomMotif(int k) {
        int range = 500 - k;
        Random random = new Random();
        int position;

        for (int i = 0; i < 10; i++) {
            position = random.nextInt(range);
            motifs[i] = gens[i].substring(position,position+10);
        }
    }

    /****************************************************************************************
    * _4RandomPositions                                                                     *
    * This function returns an ArrayList of 4 distinct random numbers                       *
    * Distinctivity is important to not to mutate the already mutated position              *
    * Firstly we create temporary list of numbers from 1 to 10                              *
    * We shuffle it                                                                         *
    * Then we insert the first 4 of them into the randomPositions list                      *
    ****************************************************************************************/
    public static ArrayList<Integer> _4RandomPositions(){
        ArrayList<Integer> list = new ArrayList<Integer>();
        for (int i=0; i<10; i++) {
            list.add(new Integer(i));
        }
        Collections.shuffle(list);
        ArrayList<Integer> randomPositions = new ArrayList<Integer>();
        for (int i=0; i<4; i++) {
            randomPositions.add(list.get(i));
        }
        return randomPositions;
    }

    /****************************************************************************************
    * createFile                                                                            *
    * This function creates a txt file if there is no file                                  *
    ****************************************************************************************/
    public static void createFile(){
        try {
            file = new File("randomMotif.txt");
            if (file.createNewFile()) {
                System.out.println("File created: " + file.getName());
            } else {
                System.out.println("File already exists.");
            }
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }

    /****************************************************************************************
    * writeRandomBasesToFile                                                                *
    * This function is very important & has a key role to this program                      *
    * Firstly, we create an array of 10 strings to store the mutated versions of k-mers     *
    * In a for loop we create & store the same DNA strings that are identical to each other *
    * We store the original k-mer in the 0th index of mutatedKMersStrings array             *
    * Then we copy this original k-mer to the remaining indexes                             *
    * Strings are immutable in Java, so we will use StringBuilder class to mutate the bases *
    * We have 10 DNA strings to mutate(1 String for each line)-> for loop will iterate to 10*
    * We'll use _4RandomPositions function to get 4 distinct positions to mutate the DNA    *
    * We put the return value of _4RandomPositions function to mutatedPositions list        *
    * 4 mutations are needed so the for loop will iterate through 4                         *
    * We will get a random index from bases string then we'll get the char at this position *
    * We get a random base to replace with, put the return value to randomBase              *
    * Now we have to find a random DNA index to mutate                                      *
    * We get the current index element from mutatedPositions list                           *
    * That is the base that we will replace with randomBase                                 *
    * But this base may be same with the randomBase, we have to control before mutating it  *
    * If they're the same, we change the randomBase with the next base in the bases string  *
    * Then we mutate the stringBuilder's base at randomPosition with randomBase             *
    * Then we update the element of mutatedKMerStrings by converting stringBuilder to string*
    * Then we start to write to file                                                        *
    * We have 10 lines and in each line we have 500 bases                                   *
    * We find a random starting index from 0 to 490 to insert the mutated DNA               *
    * We assign the return value to randStartingIndex                                       *
    * We iterate through the loop and insert randomly selected bases one by one             *
    * But when we came to the randStartingIndex, we insert the mutated DNA here             *
    * We put a newline character at the end of 500 bases.                                   *
    * We close the file and catch exceptions                                                *
    *
    *
    *
    ****************************************************************************************/
    public static void writeRandomBasesToFile(){
        int kmer=10;
        Random r = new Random();
        // There will be 10 mutated kmers in mutatedKMerStrings array
        String[] mutatedKMerStrings= new String[10];
        for (int i = 0; i < kmer ; i++) {
            mutatedKMerStrings[i]="";
            if(i!=0){
                mutatedKMerStrings[i]=mutatedKMerStrings[0];        // Copy the original k-mer to k-mers from 2nd to 10th
                continue;
            }
            for (int j = 0; j <10 ; j++) {                          //Generate the original k-mer
                mutatedKMerStrings[i]=mutatedKMerStrings[i].concat(""+bases.charAt(r.nextInt(bases.length())));
            }
        }
        System.out.println("o: "+mutatedKMerStrings[0]);
        StringBuilder stringBuilder = null;

        for (int i = 0; i <10 ; i++) {
            List<Integer> mutatedPositions=_4RandomPositions(); //DNAdaki önceden mutate edilen yerleri tekrar mutate etmemek için mutasyonlu pozisyonları burada saklıyoruz

            stringBuilder = new StringBuilder(mutatedKMerStrings[i]);
            for (int j = 0; j < 4 ; j++) {      //4 mutations needed
                int randomBaseIndex=r.nextInt(bases.length());
                char randomBase=bases.charAt(randomBaseIndex);
                int randomPosition=mutatedPositions.get(j);
                if(stringBuilder.charAt(randomPosition)==randomBase){
                    //System.out.println("--randombase: "+randomBase+" origBase: "+stringBuilder.charAt(randomPosition));
                    randomBase= bases.charAt((randomBaseIndex+1)%4);
                    //System.out.println(" -randombase: "+randomBase+" origBase: "+stringBuilder.charAt(randomPosition));
                }

                stringBuilder.setCharAt(randomPosition, randomBase);
            }
            mutatedKMerStrings[i]=stringBuilder.toString();
        }


        try {
            FileWriter myWriter = new FileWriter(file);
            for (int i = 0; i < kmer; i++) {                    //random base file will consist of 10 lines
                int randStartingIndex = r.nextInt(500-kmer);         //find a random starting index for 10-mer in the current line
                for (int j = 0; j < 500-kmer; j++) {            //each line will consist of 500 characters
                    if(j==randStartingIndex){                              //if the for iterator comes to specified random starting position
                        myWriter.write(mutatedKMerStrings[i]);  //insert the 10-mer at the specified random position
                    }
                    myWriter.write(bases.charAt(r.nextInt(bases.length())));
                }
                myWriter.write("\n");
            }

            myWriter.close();
            System.out.println("Successfully wrote to the file.");
        } catch (IOException e) {
            System.out.println("An error occurred.");
            e.printStackTrace();
        }
    }
}
