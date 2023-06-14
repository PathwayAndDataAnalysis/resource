//
// Source code recreated from a .class file by IntelliJ IDEA
// (powered by FernFlower decompiler)
//

package org.panda.resource;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.panda.utility.TermCounter;

public class UniProtSequence extends FileServer {
    private static UniProtSequence instance;
    private Map<String, String> nameToID;
    private Map<String, String> nameToSymbol;
    private Map<String, Map<String, String>> symbolToNames;
    private Map<String, String> idToSeq;

    public UniProtSequence() {
    }

    public static synchronized UniProtSequence get() {
        if (instance == null) {
            instance = new UniProtSequence();
        }

        return instance;
    }

    public String[] getLocalFilenames() {
        return new String[]{"uniprot-sequence.txt"};
    }

    public String[] getDistantURLs() {
        return new String[]{"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"};
    }

    public int getStartLocation(String nameOrID, String peptide) {
        String id = (String) this.nameToID.get(nameOrID);
        if (id == null && this.idToSeq.containsKey(nameOrID)) {
            id = nameOrID;
        }

        if (id != null) {
            String seq = (String) this.idToSeq.get(id);

            assert seq != null;

            return seq.indexOf(peptide) + 1;
        } else {
            return -1;
        }
    }

    public String getAminoacidAt(String nameOrID, int loc) {
        if (loc < 1) {
            throw new IllegalArgumentException("Location cannot be smaller than 1. loc = " + loc);
        } else {
            String id = (String) this.nameToID.get(nameOrID);
            if (id == null && this.idToSeq.containsKey(nameOrID)) {
                id = nameOrID;
            }

            if (id != null) {
                String seq = (String) this.idToSeq.get(id);

                assert seq != null;

                if (loc <= seq.length()) {
                    return seq.substring(loc - 1, loc);
                }
            }

            return null;
        }
    }

    public String getSeqAround(String nameOrID, int loc, int width) {
        if (loc < 1) {
            throw new IllegalArgumentException("Location cannot be smaller than 1. loc = " + loc);
        } else if (width % 2 != 1) {
            throw new IllegalArgumentException("Sequence width has to be an odd number. widht = " + width);
        } else if (loc <= width / 2) {
            throw new IllegalArgumentException("The location has to be greater than width/2. width = " + width + ", loc = " + loc);
        } else {
            String id = (String) this.nameToID.get(nameOrID);
            if (id == null && this.idToSeq.containsKey(nameOrID)) {
                id = nameOrID;
            }

            if (id != null) {
                String seq = (String) this.idToSeq.get(id);

                assert seq != null;

                int halfW = width / 2;
                if (loc <= seq.length() - halfW) {
                    return seq.substring(loc - halfW - 1, loc + halfW);
                }
            }

            return null;
        }
    }

    public String getIDOfName(String uniprotName) {
        return (String) this.nameToID.get(uniprotName);
    }

    /**
     * Method will return an aminoacid sequence from the desired protein at the
     * specificed location, with prefixLen amino acids preceding the location and
     * suffixLen amino acids following the location
     *
     * @param idOrName  name of protein
     * @param prefixLen number of preceding amino acids in sequence
     * @param suffixLen number of following amino acids in sequence
     * @param location  position for which user desires sequence with prefixLen amino acids prior to that location and suffixLen amino acids forward of that location
     * @return
     */
    public String getSeqAround(String idOrName, int prefixLen, int suffixLen, int location) {

		/* Obtain the longer of prefix or suffix, and then multiply by
		   2 and add one to get minimum length for getSeqAround call
		   Note that max will always be odd, as it is 2k + 1
		 */
        int max = ((prefixLen > suffixLen ? prefixLen : suffixLen) * 2) + 1;

        // Will take an odd-length as argument, and return sequence of that length
        // centered around location

        String seq = getSeqAround(idOrName, location, max);

		/*
		This occurs, for example, when the location is at the beggining(or end) of the amino
		acid sequence for the protein, such that there is not enough amino acids to left(or right)
		to recover the sequence length the user desires
		 */
        if (seq == null) {
            return null;
        }

        // Obtain the index of location
        int midPoint = max / 2;

        String prefixSeq = seq.substring(midPoint - prefixLen, prefixLen);
        String suffixSeq = seq.substring(midPoint + 1, midPoint + 1 + suffixLen);

        return prefixSeq + seq.charAt(midPoint) + suffixSeq;

    }

    public String getSymbolOfName(String name) {
        return (String) this.nameToSymbol.get(name);
    }

    public Map<String, String> getNamesOfSymbol(String symbol) {
        return (Map) this.symbolToNames.get(symbol);
    }

    public String getNameOfSymbol(String symbol, String organism) {
        return this.symbolToNames.containsKey(symbol) ? (String) ((Map) this.symbolToNames.get(symbol)).getOrDefault(organism, (Object) null) : null;
    }

    public boolean load() throws IOException {
        this.nameToID = new HashMap();
        this.nameToSymbol = new HashMap();
        this.symbolToNames = new HashMap();
        this.idToSeq = new HashMap();
        String name = null;
        String id = null;
        StringBuilder sequence = null;
        BufferedReader reader = this.getResourceReader(this.getLocalFilenames()[0]);

        for (String line = reader.readLine(); line != null; line = reader.readLine()) {
            if (line.startsWith(">")) {
                if (sequence != null) {
                    assert id != null && name != null;

                    this.idToSeq.put(id, sequence.toString());
                }

                String[] t = line.split("\\|| ");
                id = t[1];
                name = t[2];
                this.nameToID.put(name, id);
                sequence = new StringBuilder();
                int oInd = line.indexOf(" OX=");
                String organism = line.substring(oInd + 4, line.indexOf(" ", oInd + 4));
                int sInd = line.indexOf(" GN=");
                if (sInd > 0) {
                    String symbol = line.substring(sInd + 4, line.indexOf(" ", sInd + 4));
                    this.nameToSymbol.put(name, symbol);
                    if (!this.symbolToNames.containsKey(symbol)) {
                        this.symbolToNames.put(symbol, new HashMap());
                    }

                    ((Map) this.symbolToNames.get(symbol)).put(organism, name);
                }
            } else {
                sequence.append(line);
            }
        }

        reader.close();
        return true;
    }

    public static void main(String[] args) {
        String sym = "UL38";
        Map<String, String> names = get().getNamesOfSymbol(sym);
        System.out.println("names = " + names);
    }

    private static void countAAs() {
        TermCounter tc = new TermCounter();
        get().idToSeq.values().stream().forEach((s) -> {
            for (int i = 0; i < s.length(); ++i) {
                String aa = s.substring(i, i + 1);
                tc.addTerm(aa);
            }

        });
        tc.print();
    }
}
