package org.panda.resource.autismdatasets;

import org.panda.resource.FileServer;
import org.panda.utility.AltMatrixUtil;
import org.panda.utility.statistics.Histogram;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Provides genomic alterations in denovo-db.
 *
 * @author Ozgun Babur
 */
public class IossifovDataset extends FileServer
{
	private static IossifovDataset instance;

	private List<Entry> data;

	public static final DataFilter CODING = e -> !e.gene.isEmpty() && !e.effectType.equals("synonymous") &&
		!e.effectType.equals("intron");

	public static synchronized IossifovDataset get()
	{
		if (instance == null) instance = new IossifovDataset();
		return instance;
	}

	public Set<String> getAllGenes()
	{
		return data.stream().map(d -> d.gene).filter(g -> !g.equals("NA") && !g.isEmpty()).collect(Collectors.toSet());
	}

	public List<Entry> getData()
	{
		return data;
	}

	public Stream<Entry> getDataStream(DataFilter filter)
	{
		if (filter != null) return data.stream().filter(filter::select);
		else return data.stream();
	}

	public int getProbandSize()
	{
		return 2508;
	}

	public int getControlSize()
	{
		return 1911;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"Iossifov.csv"};
	}

	@Override
	public boolean load() throws IOException
	{
		data = new ArrayList<>();

		getResourceAsStream(getLocalFilenames()[0]).skip(2).forEach(l -> data.add(new Entry(l)));

		return true;
	}

	public class Entry
	{
		public String familyID;
		public boolean inProband;
		public boolean inControl;
		public String gene;
		public String effectType;

		public Entry(String line)
		{
			String[] t = line.split("\t");
			familyID = t[0];
			inProband = t[4].contains("p");
			inControl = t[4].contains("s");
			gene = t[6];
			effectType = t[7];
		}
	}

	public interface DataFilter
	{
		boolean select(Entry e);
	}

	public Map<String, Set<String>> getHitMap(DataFilter filter)
	{
		Map<String, Set<String>> hitMap = new HashMap<>();

		getDataStream(filter).filter(e -> !e.gene.isEmpty()).forEach(e ->
		{
			if (!hitMap.containsKey(e.gene)) hitMap.put(e.gene, new HashSet<>());
			hitMap.get(e.gene).add(e.familyID);
		});

		return hitMap;
	}

	public Map<String, boolean[]> getProbandAlterationMatrix(Set<String> genes, DataFilter filter)
	{
		Map<String, Set<String>> hitMap = getHitMap(filter);

		List<String> families = getDataStream(e -> e.inProband).map(e -> e.familyID).distinct().sorted().collect(Collectors.toList());

		Map<String, boolean[]> matrix = new HashMap<>();

		for (String gene : genes)
		{
			if (hitMap.containsKey(gene))
			{
				Set<String> hits = hitMap.get(gene);
				boolean[] b = new boolean[getProbandSize()];
				for (int i = 0; i < families.size(); i++)
				{
					b[i] = hits.contains(families.get(i));
				}
				matrix.put(gene, b);
			}
		}

		return matrix;
	}

	public static void main(String[] args)
	{
//		checkMutualExclusivity();
		printRowsPerFamilyHistogram();
	}

	private static void checkMutualExclusivity()
	{
		String[] ranks = {"S", "1", "S1", "2", "S2", "3", "S3", "4", "S4", "5", "S5", "6", "S6", "all"};
		Map<String, Set<String>> geneMap = SFARI.get().getGeneSetsUpToGivenRanks(ranks);

		Map<String, Integer> covMap = new HashMap<>();
		Map<String, Integer> ovMap = new HashMap<>();
		Map<String, Double> pvals = new HashMap<>();

		geneMap.keySet().parallelStream().forEach(rank ->
		{
			Map<String, boolean[]> matrix = get().getProbandAlterationMatrix(geneMap.get(rank), CODING);
			List<boolean[]> m = AltMatrixUtil.toList(matrix);

			covMap.put(rank, AltMatrixUtil.getCoverage(m));
			ovMap.put(rank, AltMatrixUtil.getOverlap(m));

			double p = AltMatrixUtil.getOverallMutexnessPValue(m, 1000);

			pvals.put(rank, p);
		});

		System.out.println("SFARI Rank\tGenes size\tCoverage\tOverlap\tP-value");
		Arrays.stream(ranks).forEach(rank -> System.out.println(rank + "\t" + geneMap.get(rank).size() + "\t" +
			covMap.get(rank) + "\t" + ovMap.get(rank) + "\t" + pvals.get(rank)));
	}

	private static void printRowsPerFamilyHistogram()
	{
		Map<String, Long> map = get().getDataStream(e -> true).collect(Collectors.groupingBy(e -> e.familyID, Collectors.counting()));
		Histogram h = new Histogram(1);
		h.setBorderAtZero(true);
		h.setUseLowerBorderForPrinting(true);
		map.values().forEach(h::count);
		h.print();
	}
}
