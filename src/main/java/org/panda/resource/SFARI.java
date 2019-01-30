package org.panda.resource;

import java.io.IOException;
import java.sql.SQLOutput;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Provides genes in SFARI for Autism.
 *
 * @author Ozgun Babur
 */
public class SFARI extends FileServer
{
	private static SFARI instance;

	private Map<String, Entry> data;

	public static synchronized SFARI get()
	{
		if (instance == null) instance = new SFARI();
		return instance;
	}

	public Set<String> getAllGenes()
	{
		return data.keySet();
	}

	public boolean isAutismGene(String sym)
	{
		return data.containsKey(sym);
	}

	/**
	 * Low score is better.
	 */
	public Set<String> getGenesWithMaxScore(int scoreCutoff)
	{
		return data.values().stream().filter(e -> e.scoreSatisfy(scoreCutoff)).map(e -> e.gene).collect(Collectors.toSet());
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"SFARI.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"https://gene.sfari.org//wp-content/themes/sfari-gene/utilities/download-csv.php?api-endpoint=genes"};
	}

	@Override
	public boolean load() throws IOException
	{
		data = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0]).skip(1).map(l ->
			l.split(",(?=([^\\\"]*\\\"[^\\\"]*\\\")*[^\\\"]*$)")).forEach(t ->
				data.put(t[1], new Entry(t[1], t[3], t[5].isEmpty() ? null : Integer.valueOf(t[5]),
					t[6].equals("1"), Integer.valueOf(t[7]))));

		return true;
	}

	public void plot()
	{
		for (Entry entry : data.values())
		{
//			if (entry.hasScore())
				System.out.println((entry.score == null? 0 : entry.score) + "\t" + entry.numOfReports);
		}
	}

	class Entry
	{
		String gene;
		String chromosome;
		Integer score;
		boolean syndromic;
		int numOfReports;

		public Entry(String gene, String chromosome, Integer score, boolean syndromic, int numOfReports)
		{
			this.gene = gene;
			this.chromosome = chromosome;
			this.score = score;
			this.syndromic = syndromic;
			this.numOfReports = numOfReports;
		}

		boolean scoreSatisfy(int thr)
		{
			return hasScore() && score <= thr;
		}

		public boolean hasScore()
		{
			return score != null;
		}
	}

	public static void main(String[] args)
	{
//		get().plot();
//		get().getGenesWithMaxScore(4).stream().sorted().forEach(System.out::println);
		get().getAllGenes().forEach(System.out::println);
	}

}
