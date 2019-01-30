package org.panda.resource;

import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.Graph;
import org.panda.utility.graph.UndirectedGraph;

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

	public String getIDOfName(String uniprotName)
	{
		return nameToID.get(uniprotName);
	}

	@Override
	public boolean load() throws IOException
	{
		nameToID = new HashMap<>();
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
		int startLocation = get().getStartLocation("P62158", "vFDk".toUpperCase());
		System.out.println("startLocation = " + startLocation);
	}
}
