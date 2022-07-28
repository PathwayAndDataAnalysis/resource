package org.panda.resource;

import org.panda.utility.TermCounter;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.*;

/**
 * This class provides a mapping between HGNC IDs and approved gene symbols.
 *
 * @author Ozgun Babur
 */
public class UniProtSequence extends FileServer
{
	private static UniProtSequence instance;

	private Map<String, String> nameToID;
	private Map<String, String> nameToSymbol;
	//The map is from organism to UP name
	private Map<String, Map<String, String>> symbolToNames;
	private Map<String, String> idToSeq;

	public static synchronized UniProtSequence get()
	{
		if (instance == null) instance = new UniProtSequence();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"uniprot-sequence.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/" +
			"uniprot_sprot.fasta.gz"};
	}

	/**
	 * Starts from 1.
	 */
	public int getStartLocation(String nameOrID, String peptide)
	{
		String id = nameToID.get(nameOrID);

		if (id == null && idToSeq.containsKey(nameOrID))
		{
			id = nameOrID;
		}

		if (id != null)
		{
			String seq = idToSeq.get(id);
			assert seq != null;

			return seq.indexOf(peptide) + 1;
		}

		return -1;
	}

	public String getAminoacidAt(String nameOrID, int loc)
	{
		if (loc < 1) throw new IllegalArgumentException("Location cannot be smaller than 1. loc = " + loc);

		String id = nameToID.get(nameOrID);

		if (id == null && idToSeq.containsKey(nameOrID))
		{
			id = nameOrID;
		}

		if (id != null)
		{
			String seq = idToSeq.get(id);
			assert seq != null;

			if (loc <= seq.length())
			{
				return seq.substring(loc - 1, loc);
			}
		}

		return null;
	}

	public String getSeqAround(String nameOrID, int loc, int width)
	{
		if (loc < 1) throw new IllegalArgumentException("Location cannot be smaller than 1. loc = " + loc);
		if (width % 2 != 1) throw new IllegalArgumentException("Sequence width has to be an odd number. widht = " + width);
		if (loc <= width / 2) throw new IllegalArgumentException("The location has to be greater than width/2. width = " + width + ", loc = " + loc);

		String id = nameToID.get(nameOrID);

		if (id == null && idToSeq.containsKey(nameOrID))
		{
			id = nameOrID;
		}

		if (id != null)
		{
			String seq = idToSeq.get(id);
			assert seq != null;

			int halfW = width / 2;
			if (loc <= seq.length() - halfW)
			{
				return seq.substring(loc - halfW - 1, loc + halfW);
			}
		}
		return null;
	}

	public String getIDOfName(String uniprotName)
	{
		return nameToID.get(uniprotName);
	}

	public String getSymbolOfName(String name)
	{
		return nameToSymbol.get(name);
	}

	public Map<String, String> getNamesOfSymbol(String symbol)
	{
		return symbolToNames.get(symbol);
	}

	public String getNameOfSymbol(String symbol, String organism)
	{
		return symbolToNames.containsKey(symbol) ? symbolToNames.get(symbol).getOrDefault(organism, null) : null;
	}

	/**
	 * Does not read the last sequence, which belongs to Whitewater Arroyo mammarenavirus.
	 */
	@Override
	public boolean load() throws IOException
	{
		nameToID = new HashMap<>();
		nameToSymbol = new HashMap<>();
		symbolToNames = new HashMap<>();
		idToSeq = new HashMap<>();

		String name = null;
		String id = null;
		StringBuilder sequence = null;

		BufferedReader reader = getResourceReader(getLocalFilenames()[0]);

		String line = reader.readLine();

		while (line != null)
		{
			if (line.startsWith(">"))
			{
				if (sequence != null)
				{
					assert id != null && name != null;
					idToSeq.put(id, sequence.toString());
				}

				String[] t = line.split("\\|| ");
				id = t[1];
				name = t[2];
				nameToID.put(name, id);
				sequence = new StringBuilder();

				int oInd = line.indexOf(" OX=");
				String organism = line.substring(oInd + 4, line.indexOf(" ", oInd + 4));

				int sInd = line.indexOf(" GN=");
				if (sInd > 0)
				{
					String symbol = line.substring(sInd + 4, line.indexOf(" ", sInd + 4));
					nameToSymbol.put(name, symbol);
					if (!symbolToNames.containsKey(symbol)) symbolToNames.put(symbol, new HashMap<>());
					symbolToNames.get(symbol).put(organism, name);
				}
			}
			else
			{
				sequence.append(line);
			}

			line = reader.readLine();
		}


		reader.close();

		return true;
	}


	public static void main(String[] args)
	{
//		int startLocation = get().getStartLocation("P0DP23", "vFDk".toUpperCase());
//		System.out.println("startLocation = " + startLocation);

		String sym = "UL38";
		Map<String, String> names = get().getNamesOfSymbol(sym);
		System.out.println("names = " + names);
//		System.out.println(get().getSeqAround("P0DP23", 80, 5));

//		countAAs();

//		System.out.println(get().getSymbolOfName("P53_HUMAN"));
//		System.out.println(get().getNamesOfSymbol("TP53"));
	}

	private static void countAAs()
	{
		TermCounter tc = new TermCounter();
		get().idToSeq.values().stream().forEach(s ->
		{
			for (int i = 0; i < s.length(); i++)
			{
				String aa = s.substring(i, i + 1);
				tc.addTerm(aa);
			}
		});
		tc.print();
	}
}
