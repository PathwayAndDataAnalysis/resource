package org.panda.resource.autismdatasets;

import org.panda.resource.CancerGeneCensus;
import org.panda.resource.FileServer;
import org.panda.resource.OncoKB;

import java.io.IOException;
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

	public boolean isSyndromic(String gene)
	{
		Entry entry = data.get(gene);
		if (entry != null)
		{
			return entry.syndromic;
		}

		return false;
	}

	/**
	 * Low score is better.
	 */
	public Set<String> getGenesWithMaxScore(int scoreCutoff)
	{
		return data.values().stream().filter(e -> e.scoreSatisfy(scoreCutoff)).map(e -> e.gene).collect(Collectors.toSet());
	}

	public Set<String> getGenesWithExactScore(int score)
	{
		return data.values().stream().filter(e -> score == 0 ? !e.hasScore() : e.hasScore() && e.score == score)
			.map(e -> e.gene).collect(Collectors.toSet());
	}

	public Set<String> getGenesSyndromic()
	{
		return data.values().stream().filter(e -> e.syndromic).map(e -> e.gene).collect(Collectors.toSet());
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

	public String getClassification(String gene)
	{
		Entry entry = data.get(gene);
		if (entry == null) return "";

		String s = entry.syndromic ? "S" : "";
		if (entry.score != null)
		{
			s = entry.score + s;
		}

		if (s.isEmpty()) s = "?";
		return s;
	}

	public Map<String, Set<String>> getGeneSetsUpToGivenRanks(String... ranks)
	{
		Map<String, Set<String>> genesMap = new HashMap<>();

		Arrays.stream(ranks).forEach(rank ->
		{
			Set<String> genes = rank.equals("all") ? SFARI.get().getAllGenes() :
				rank.equals("S") ? SFARI.get().getGenesSyndromic() :
					rank.startsWith("S") ? SFARI.get().getGenesWithMaxScore(Integer.valueOf(rank.substring(1))) :
						SFARI.get().getGenesWithMaxScore(Integer.valueOf(rank));

			if (rank.startsWith("S") && rank.length() > 1) genes.addAll(SFARI.get().getGenesSyndromic());

			genesMap.put(rank, genes);
		});
		return genesMap;
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
		printStats();
//		System.out.println(get().isAutismGene("RBFOX1"));
	}

	private static void printStats()
	{
		Set<String> cancer = new HashSet<>(CancerGeneCensus.get().getAllSymbols());
		cancer.addAll(OncoKB.get().getAllSymbols());

		for (int i = 0; i < 7; i++)
		{
			Set<String> genes = SFARI.get().getGenesWithExactScore(i);
			Set<String> synd = genes.stream().filter(g -> SFARI.get().isSyndromic(g)).collect(Collectors.toSet());
			genes.removeAll(synd);
			System.out.println(i + "\t" + genes.size() + "\t" + synd.size());
			System.out.println(synd.stream().filter(cancer::contains).sorted().collect(Collectors.toList()) + "\t" + genes.stream().filter(cancer::contains).sorted().collect(Collectors.toList()));
		}

		System.out.println("All genes in SFARI = " + SFARI.get().getAllGenes().size());
	}

}
