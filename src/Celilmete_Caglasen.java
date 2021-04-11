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

    public static void main(String args[]){
        int k = 10;
        createFile();
//        writeRandomBasesToFile();
        readGens();
        System.out.println(randomizedMotifSearch(k));
        printMotifs();



    }

    public static void printProfile() {
        for (int i = 0; i < profile.length; i++) {
            for (int j = 0; j < profile[i].length; j++) {
                System.out.printf("%.1f ", profile[i][j]);
            }
            System.out.println();
        }
    }

    public static void printMotifs() {
        System.out.println("---------");
        for (int i = 0; i < motifs.length; i++) {
            System.out.println(motifs[i]);
        }
        System.out.println("---------");
    }

    public static int randomizedMotifSearch(int k) {
        getRandomMotif(k);
        String[] bestMotifs = motifs.clone();
        int bestScore = score();
        int score;
        int count = 1;
        while (true){
            getRandomMotif(k);
            String[] tempMotifs = motifs.clone();
            int tempScore = 99999;
            while (true) {
                getProfile();
                getBestKMers(k);
                score = score();
                if (score < tempScore) {
                    tempMotifs = motifs.clone();
                    tempScore = score;
                }
                else {
                    break;
                }
            }
            if(tempScore < bestScore) {
                bestScore = tempScore;
                bestMotifs = tempMotifs.clone();
                count = 0;
            }
            else if (count % 50 == 0) {
                return bestScore;
            }
            else {
                motifs = bestMotifs.clone();
                count++;
            }
        }

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

    public static int baseToIndex(char base) {
        return "ACGT".indexOf(base);
    }

    public static char indexToBase(int i) {
        return "ACGT".charAt(i);
    }

    public static double getProb(int k, String kmer) {
        double prob = 1;
        for (int i = 0; i < kmer.length(); i++) {
            prob *= profile[baseToIndex(kmer.charAt(i))][i];
        }
        return prob;
    }

    public static void  getBestKMers(int k) {
        for (int i = 0; i < gens.length; i++) {
            double bestVal = -1;
            String bestKmer = "";
            for (int j = 0; j < gens[i].length() - k + 1; j++) {
                String kmer = gens[i].substring(j, j + k);
                double prob = getProb(k, kmer);
                if (bestVal < getProb(k, kmer)) {
                    bestKmer = kmer;
                    bestVal = prob;
                }
            }
            motifs[i] = bestKmer;
        }
    }

    /************************************************************************
     * This function gets motifs starting from random positions in each line */
    public static void getRandomMotif(int k) {
        int range = 500 - k;
        Random random = new Random();
        int position;

        for (int i = 0; i < 10; i++) {
            position = random.nextInt(range);
            motifs[i] = gens[i].substring(position,position+k);
        }
    }

    public static int score() {
        int score = 0;
        for (int i = 0; i < motifs.length; i++) {
            ArrayList<Integer> counts = new ArrayList<>();
            int numa = 0, numc = 0, numg = 0, numt = 0;
            for (int j = 0; j < motifs[i].length(); j++) {
                switch (motifs[j].charAt(i)) {
                    case 'A' -> numa++;
                    case 'C' -> numc++;
                    case 'G' -> numg++;
                    case 'T' -> numt++;
                }
            }
            counts.add(numa);counts.add(numc);counts.add(numg);counts.add(numt);
            Collections.sort(counts);
            score += 10 - counts.get(3);
            counts.clear();
        }

        return score;
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
