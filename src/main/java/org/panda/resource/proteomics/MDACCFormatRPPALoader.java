package org.panda.resource.proteomics;

import org.panda.resource.HGNC;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * For loading a standard format MDACC RPPA dataset. THe data comes multi-sheet XSLX file. The last sheet
 * (NormLog2_MedianCentered) should be saved as tab-delimited text file (choose .cvs). Then it can be loaded using this
 * class.
 */
public class MDACCFormatRPPALoader
{
	public static Map<String, Map<String, Double>> load(String filename) throws IOException
	{
		String[] abs = Files.lines(Paths.get(filename)).filter(l -> l.contains("\tAntibody Name\t")).findFirst().get().split("\t");
		String[] geneNames = Files.lines(Paths.get(filename)).filter(l -> l.contains("\tGene Name\t")).findFirst().get().split("\t");
		List<String> colNames = Arrays.asList(Files.lines(Paths.get(filename)).filter(l -> l.contains("\tSample description\t")).findFirst().get().split("\t"));

		int sampleIndex = colNames.indexOf("Sample description");
		if (sampleIndex < 0) throw new RuntimeException("Sample description not found.");

		int startCol = Arrays.asList(abs).indexOf("Antibody Name") + 1;
		if (sampleIndex < 1) throw new RuntimeException("Antibody Name not found.");

		boolean halt = false;
		// check if there is missing ab
		List<Integer> inds = new ArrayList<>();
		for (int i = startCol; i < abs.length; i++)
		{
			String preferredID = RPPAIDMapper.get().getPreferredID(abs[i]);
			if (preferredID == null)
			{
				halt = true;
				if (abs[i].contains("_p") || HGNC.get().getSymbol(geneNames[i]) == null)
					System.err.println("Not found: " + abs[i] + "   gene = " + geneNames[i]);
				else inds.add(i);
			}
		}
		if (!inds.isEmpty()) System.out.println("\nAdd these to the ID mapper file:\n");
		inds.forEach(i -> System.out.println(abs[i] + "\t" + abs[i] + "\t" + geneNames[i]));

		halt = halt || !inds.isEmpty();

		if (halt)
		{
			System.out.println();
			throw new RuntimeException("There we unrecognized antibodies. See above lines.");
		}
		//--------

		Map<String, Map<String, Double>> map = new HashMap<>();

		Files.lines(Paths.get(filename)).skip(11).map(l -> l.split("\t")).forEach(t ->
		{
			String sample = t[sampleIndex];

			for (int i = startCol; i < t.length; i++)
			{
				String ab = RPPAIDMapper.get().getPreferredID(abs[i]);
				Double val = Double.valueOf(t[i]);

				if (!map.containsKey(ab)) map.put(ab, new HashMap<>());
				map.get(ab).put(sample, val);
			}
		});

		return map;
	}

	public static void main(String[] args) throws IOException
	{
		Map<String, Map<String, Double>> map = load("/home/ozgun/Data/LINCS/repA.csv");
	}
}

