package org.panda.resource;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;


public class KinaseLibrary extends FileServer {

    /*
    A map which can be used to determine the probability that a kinase has some
    amino acid at some position
     */
    private Map<String, Map<Integer, Map<String, Double>>> probabilityMap;

    private int maximumSiteVal;

    private int minimumSiteVal;

    public static final int validPeptideLength = 10;

    public static final int phosphoAcceptorIndex = 5;

    @Override
    public String[] getLocalFilenames() {
        return new String[]{"KinaseLibraryMatrices.tsv"};
    }

    // Returns the minimum location value
    public int minimumLocationValue() {
        return minimumSiteVal;
    }

    // Returns the maximum location value
    public int maximumLocationValue() {
        return maximumSiteVal;
    }


    // Returns all kinases within the map
    public Set<String> kinaseSet() {
        return probabilityMap.keySet();
    }

    private double calculateScore(String kinase, HashMap<Integer, Character> locationAminoAcid) {
        double probability = 1;
        for (int location : locationAminoAcid.keySet()) {
            if (location != 0) {
                String s = Character.toString(locationAminoAcid.get(location));
                double value = probabilityMap.get(kinase).get(location).get(Character.toString(locationAminoAcid.get(location)));
                probability = probability * value;
            }
        }

        return probability;
    }


    private HashMap<Integer, Character> processPeptide(String aminoAcidSequence) {

        // Produce a new hashmap that will contain the location as the key and the corresponding
        // amino acid that that location as the value
        HashMap<Integer, Character> locationAminoAcid = new HashMap<>();

        // Obtain the index of the phosphorylated residue
        int phosphoResidueIndex = findPhosphorylatedResidue(aminoAcidSequence);

        // The phosphorylated residue is at the 0th location in the peptide
        locationAminoAcid.put(0, aminoAcidSequence.charAt(phosphoResidueIndex));


        // Splitting into negative and positive locations
        String negativeLocations = aminoAcidSequence.substring(0, phosphoResidueIndex);
        String positiveLocations = aminoAcidSequence.substring(phosphoResidueIndex + 1);


        // Process the negative locations into hashmap
        for (int i = 0; i < negativeLocations.length(); i++) {
            locationAminoAcid.put(i - negativeLocations.length(), negativeLocations.charAt(i));
        }
        // Process the positive locations into hashmap
        for (int i = 0; i < positiveLocations.length(); i++) {
            locationAminoAcid.put(i + 1, positiveLocations.charAt(i));
        }

        return locationAminoAcid;
    }


    /* Method to access the value for each amino acid from the map
     */
    private double aminoAcidProbability(String kinase, int location, String aminoAcid) {

        // Check for user requesting information that does not exist in the map
        if (!probabilityMap.containsKey(kinase)) {
            throw new IllegalArgumentException("No such kinase available.");
        } else if (!probabilityMap.get(kinase).containsKey(location)) {
            throw new IllegalArgumentException("No such location available.");
        } else if (!probabilityMap.get(kinase).get(location).containsKey(aminoAcid)) {
            throw new IllegalArgumentException("No such amino acid available");
        }

        return probabilityMap.get(kinase).get(location).get(aminoAcid);

    }

    /*
    Method reads from KinaseLibraryMatrices file and inserts data into
    probabilityMap instance variable
     */
    @Override
    public boolean load() throws IOException {
        // Get the path to the file
        Path p = (Paths.get(this.locateInBase(this.getLocalFilenames()[0])));

        // Read all the lines to begin processing
        List<String> lines = Files.readAllLines(p);

        // Process the first line in order to get column headers (do not hardcode these)
        String[] firstLineSplit = lines.get(0).split("\t");

        // This next part of the code will find the locations/amino acids from the column headers
        // We will store these in two arrays

        int[] locations = new int[lines.size()];
        String[] aminoAcids = new String[lines.size()];

        splitColumnHeaders(firstLineSplit, locations, aminoAcids);

        // Assign maximum location value and minimum location value
        maximumSiteVal = findMax(locations);
        minimumSiteVal = findMin(locations);

        probabilityMap = new HashMap<>();

        /* Process file row-by-row to fill out map
           Note we have already processed 0th row,
           for column header information
         */

        for (int row = 1; row < lines.size(); row++) {
            // Split the current row or line by the delimiter tab
            String[] splitLine = lines.get(row).split("\t");

            // Determine the kinase for curr row, and put it into the kinase-Location map
            String currKinase = splitLine[0];

            probabilityMap.put(currKinase, new HashMap<>());

            for (int cell = 1; cell < splitLine.length; cell++) {

                // Use the column headers to determine the location/amino acid for current cell
                int currLocation = locations[cell];
                String currAminoAcid = aminoAcids[cell];
                double measurement = Double.parseDouble(splitLine[cell]);


                // Put the measurement into the map appropriately
                if (!probabilityMap.get(currKinase).containsKey(currLocation)) {
                    probabilityMap.get(currKinase).put(currLocation, new HashMap<>());
                }

                probabilityMap.get(currKinase).get(currLocation).put(currAminoAcid, measurement);

            }


        }

        return true;
    }


    private int findMax(int[] arr) {
        int max = arr[0];
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] > max) {
                max = arr[i];
            }
        }
        return max;
    }

    private int findMin(int[] arr) {
        int min = arr[0];
        for (int i = 0; i < arr.length; i++) {
            if (arr[i] < min) {
                min = arr[i];
            }
        }
        return min;
    }

    /*
    Method to process first row of data file. Will read all the column headers and deduce
    the locations and amino acid reffered to by each column header
     */
    private void splitColumnHeaders(String[] firstLine, int[] sites, String[] aminoAcids) {
        for (int i = 0; i < firstLine.length; i++) {
            int startIndex = 0;
            boolean foundDivider = false;
            while (startIndex < firstLine[i].length() && !foundDivider)
                if (firstLine[i].charAt(startIndex) == '-' || Character.isDigit(firstLine[i].charAt(startIndex))) {
                    startIndex++;
                } else {
                    int siteValue = Integer.parseInt(firstLine[i].substring(0, startIndex));
                    String aminoAcidValue = firstLine[i].substring(startIndex);
                    sites[i] = siteValue;
                    aminoAcids[i] = aminoAcidValue;
                    foundDivider = true;
                }

        }
    }


    // Method that validates the peptide length
    private boolean validatePeptideLength(String aminoAcidSequence) {
        return aminoAcidSequence.length() == validPeptideLength;
    }


    /* Method that returns the index of the phosphorylated residue
       Note that in the future, the phosphoAcceptorIndex may change
       assuming KinaseLibraryMatrices is updated
     */
    private int findPhosphorylatedResidue(String aminoAcidSequence) {
        return phosphoAcceptorIndex;
    }

    /*
    Method which takes an amino acid sequence as the argument and will
    return a map of probabilities for putative kinases.

    The input sequence should be in the format:
    XXXXXsXXXX or XXXXXtXXXX where s or t are the phosphorylated residue

    In future, it may become possible to score more varied amino acid sequences.
    I.e. Sequences with length of 12. In this case both validatePeptideLength and
    findPhosphorylatedResidue must be updated appropriately.
     */
    public HashMap<String, Double> peptideScore(String aminoAcidSequence) {
        // Process the amino acid sequence by splitting into negative and positive locations

        if (!validatePeptideLength(aminoAcidSequence)) {
            throw new IllegalArgumentException("Invalid sequence: sequence must be length 10");
        }


        // A hashmap which contains each amino acid in amino acid sequence as values, and the
        // corresponding key will be the location of that amino acid in the sequence passed as argument
        HashMap<Integer, Character> locationAminoAcid = processPeptide(aminoAcidSequence);

        // Create a hashmap to store the scores for each kinase
        HashMap<String, Double> kinaseScore = new HashMap<String, Double>();

        // Calculate a score for each kinase and put it into the map
        for (String kinase : probabilityMap.keySet()) {
            double calculatedScore = calculateScore(kinase, locationAminoAcid);
            kinaseScore.put(kinase, calculatedScore);
        }

        return kinaseScore;


    }


}


