package org.panda.causalpath.resource;

import org.panda.resource.FileServer;

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

    public static final int validPeptideLength = 11;

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
        String negativeLocations = aminoAcidSequence.substring(0, phosphoResidueIndex - 1);
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


    /* Public method for user to access the probability that if kinase phosphorlyated
     a given peptide, what are the odds the peptide has given aminoAcid at specified location
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


    public static void main(String[] args) {

        // Load

        KinaseLibrary kN = new KinaseLibrary();
        try {
            kN.load();
        } catch (Exception e) {

        }


        System.out.println(kN.testHandler());




    }


    private boolean testProbabilityMap(String kinase, String[] tableData, String[] correspondingColumnHeads) {

        boolean test = true;

        int[] sites = new int[correspondingColumnHeads.length];
        String[] aminoAcids = new String[correspondingColumnHeads.length];
        splitColumnHeaders(correspondingColumnHeads, sites, aminoAcids);

        for (int i = 0; i < tableData.length; i++) {
            double observedData = Double.parseDouble(tableData[i]);
            int correspondingLocation = sites[i];
            String correspondingAmino = aminoAcids[i];

            if (observedData != aminoAcidProbability(kinase, correspondingLocation, correspondingAmino)) {
                System.out.println("Test failed at kinase " + kinase + "observed/recorded " + observedData + " " +
                        probabilityMap.get(kinase).get(correspondingLocation).get(correspondingAmino));
                test = false;
            }

        }

        return test;

    }


    private boolean testHandler() {

        // Test #1
        String[] ds1 = {"0.7382", "1.1795", "1.0834", "1.0338", "0.9636", "0.9636", "0.7587", "0.8737", "0.7526", "0.7223", "0.8296", "0.9636", "0.9062", "1.0807", "1.3753", "1.6725", "0.8025", "1.2005", "1.0043", "1.0561", "1.04", "1.04", "1.3041", "0.905", "1.2345", "1.0772", "0.9535", "0.9518", "0.9518", "0.6988", "0.7062", "0.7658", "0.8184", "0.7252", "0.7825", "0.966", "0.9518", "1.3531", "1.7077", "0.9085", "1.0403", "1.3963", "0.9627", "0.9723", "0.9723", "0.9622"};
        String[] ch1 = {"-5P", "-5G", "-5A", "-5C", "-5S", "-5T", "-5V", "-5I", "-5L", "-5M", "-5F", "-5Y", "-5W", "-5H", "-5K", "-5R", "-5Q", "-5N", "-5D", "-5E", "-5s", "-5t", "-5y", "-4P", "-4G", "-4A", "-4C", "-4S", "-4T", "-4V", "-4I", "-4L", "-4M", "-4F", "-4Y", "-4W", "-4H", "-4K", "-4R", "-4Q", "-4N", "-4D", "-4E", "-4s", "-4t", "-4y", "-3P", "-3G", "-3A", "-3C", "-3S", "-3T", "-3V", "-3I", "-3L", "-3M", "-3F", "-3Y", "-3W", "-3H", "-3K", "-3R", "-3Q", "-3N", "-3D", "-3E", "-3s", "-3t", "-3y", "-2P", "-2G", "-2A", "-2C", "-2S", "-2T", "-2V", "-2I", "-2L", "-2M", "-2F", "-2Y", "-2W", "-2H", "-2K", "-2R", "-2Q", "-2N", "-2D", "-2E", "-2s", "-2t", "-2y", "-1P", "-1G", "-1A", "-1C", "-1S", "-1T", "-1V", "-1I", "-1L", "-1M", "-1F", "-1Y", "-1W", "-1H", "-1K", "-1R", "-1Q", "-1N", "-1D", "-1E", "-1s", "-1t", "-1y", "1P", "1G", "1A", "1C", "1S", "1T", "1V", "1I", "1L", "1M", "1F", "1Y", "1W", "1H", "1K", "1R", "1Q", "1N", "1D", "1E", "1s", "1t", "1y", "2P", "2G", "2A", "2C", "2S", "2T", "2V", "2I", "2L", "2M", "2F", "2Y", "2W", "2H", "2K", "2R", "2Q", "2N", "2D", "2E", "2s", "2t", "2y", "3P", "3G", "3A", "3C", "3S", "3T", "3V", "3I", "3L", "3M", "3F", "3Y", "3W", "3H", "3K", "3R", "3Q", "3N", "3D", "3E", "3s", "3t", "3y", "4P", "4G", "4A", "4C", "4S", "4T", "4V", "4I", "4L", "4M", "4F", "4Y", "4W", "4H", "4K", "4R", "4Q", "4N", "4D", "4E", "4s", "4t", "4y"};
        String k1 = "AURA";
        boolean t1Result = testProbabilityMap(k1, ds1, ch1);


        // Test #2
        String[] ds2 = {"1.3084", "1.1278", "1.1916", "0.9525", "0.8012", "0.8012", "0.6998", "0.5772", "0.6924", "0.9278", "0.6697", "0.7127", "0.6657", "1.2586", "2.2192", "2.3084", "0.8742", "0.8012", "0.4195", "0.5458", "0.6436", "0.6436", "0.475"};
        String[] ch2 = {"-4P", "-4G", "-4A", "-4C", "-4S", "-4T", "-4V", "-4I", "-4L", "-4M", "-4F", "-4Y", "-4W", "-4H", "-4K", "-4R", "-4Q", "-4N", "-4D", "-4E", "-4s", "-4t", "-4y"};
        String k2 = "AMPKA2";
        boolean t2Result = testProbabilityMap(k2, ds2, ch2);


        return t1Result && t2Result;

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

    // Method that validates that the phosphorylated site is in the appropriate position
    private boolean validatePhosphorylationSite(String aminoAcidSequence) {
        return (findPhosphorylatedResidue(aminoAcidSequence) == 6);
    }

    // Method that returns the index of the phosphorylated residue
    private int findPhosphorylatedResidue(String aminoAcidSequence) {

        // Variable to store index of the phosphorylated residue in the sequence
        int phosphoResidueIndex = -1;

        // Iterate through sequence to find phosphorylated residue pX
        for (int i = 0; i < aminoAcidSequence.length(); i++) {
            if (aminoAcidSequence.charAt(i) == 'p') {
                phosphoResidueIndex = ++i;
                break;
            }
        }

        return phosphoResidueIndex;
    }

    /*
    Method which takes an amino acid sequence as the argument and will
    return a map of probabilities for putative kinases.

    The input sequence should be in the format:
    XXXXXpSXXXX or XXXXXpTXXXX where pS or pT are the phosphorylated residue

    In future, it may become possible to score more varied amino acid sequences.
    In this case, alteration of validatePeptideLength necessary.
     */
    public HashMap<String, Double> peptideScore(String aminoAcidSequence) {
        // Process the amino acid sequence by splitting into negative and positive locations

        if (!validatePeptideLength(aminoAcidSequence)) {
            throw new IllegalArgumentException("Invalid sequence: XXXXXpSXXXX or XXXXXpTXXXX format required");
        }

        if (!validatePhosphorylationSite(aminoAcidSequence)) {
            throw new IllegalArgumentException("Phosphorylated residue must be preceded by 'p'. Phosphorylated residue must occur" +
                    " at the 6th index");
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


