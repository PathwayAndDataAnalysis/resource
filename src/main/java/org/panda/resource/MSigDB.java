package org.panda.resource;

import org.biopax.paxtools.controller.Cloner;
import org.biopax.paxtools.controller.Completer;
import org.biopax.paxtools.controller.SimpleEditorMap;
import org.biopax.paxtools.io.SimpleIOHandler;
import org.biopax.paxtools.model.BioPAXElement;
import org.biopax.paxtools.model.BioPAXLevel;
import org.biopax.paxtools.model.Model;
import org.biopax.paxtools.model.level3.EntityReference;
import org.biopax.paxtools.model.level3.Pathway;
import org.biopax.paxtools.model.level3.Xref;
import org.biopax.paxtools.pattern.util.Blacklist;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.Kronometre;
import org.panda.utility.statistics.FDR;
import org.panda.utility.statistics.FishersExactTest;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Provides gene sets from MSigDB.
 *
 * @author Ozgun Babur
 */
public class MSigDB extends FileServer
{
	private static MSigDB instance;
	private static final String FILE = "msigdb.v7.2.symbols.gmt";

	private Map<String, Set<String>> geneSets;
	private Map<String, String> urls;

	public static synchronized MSigDB get()
	{
		if (instance == null) instance = new MSigDB();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{FILE};
	}

	@Override
	public boolean load() throws IOException
	{
		geneSets = new HashMap<>();
		urls = new HashMap<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(line -> line.split("\t")).forEach(t ->
		{
			String name = t[0];
			String url = t[1];

			urls.put(name, url);

			Set<String> genes = new HashSet<>(Arrays.asList(t).subList(2, t.length));
			geneSets.put(name, genes);
		});

		return true;
	}

	public void crop(NameFilter filter)
	{
		for (String name : new HashSet<>(geneSets.keySet()))
		{
			if (!filter.select(name)) geneSets.remove(name);
		}
	}

	public void crop(Set<String> background, int minCountToKeepASet)
	{
		Set<String> remove = new HashSet<>();
		geneSets.forEach((name, set) ->
		{
			set.retainAll(background);
			if (set.size() < minCountToKeepASet) remove.add(name);
		});
		remove.forEach(geneSets::remove);
	}

	public Map<String, Set<String>> getSetsNameFiltered(NameFilter filter)
	{
		return geneSets.keySet().stream().filter(filter::select).collect(Collectors.toMap(Function.identity(), geneSets::get));
	}

	public void printSetsNameFiltered(NameFilter filter)
	{
		geneSets.keySet().stream().filter(filter::select).forEach(System.out::println);
	}

	public String getURL(String name)
	{
		return urls.get(name);
	}

	public int getSize(String name)
	{
		return geneSets.get(name).size();
	}

	public Set<String> getGeneSet(String name)
	{
		return geneSets.get(name);
	}

	/**
	 * Gets the enrichment pvals and pval limits.
	 * @return two maps, first is for pvals, second is for limits
	 */
	public Map<String, Double>[] getEnrichmentPvals(Set<String> selectedGenes)
	{
		// Filter out genes in the selected genes that are not part of any gene sets
		Set<String> background = geneSets.values().stream().flatMap(Set::stream).collect(Collectors.toSet());
		selectedGenes.retainAll(background);

		Map<String, Integer> selectionGeneCnt = count(selectedGenes);

		Map<String, Double> mapP = new HashMap<>();
		Map<String, Double> mapL = new HashMap<>();

		Set<Set<String>> memory = new HashSet<>();

		selectionGeneCnt.keySet().forEach(name ->
		{
			// If two gene sets are exactly the same, we don't want to test them multiple times
			Set<String> geneSet = geneSets.get(name);
			if (memory.contains(geneSet)) return;
			else memory.add(geneSet);

			int size = background.size();
			int featuredOverall = geneSets.get(name).size();
			int selected = selectedGenes.size();
			int featuredSelected = selectionGeneCnt.get(name);

			double pval = FishersExactTest.calcEnrichmentPval(size, featuredOverall, selected, featuredSelected);

			int maxPossibleHit = Math.min(featuredOverall, selected);

			double limit = FishersExactTest.calcEnrichmentPval(size, featuredOverall, selected, maxPossibleHit);

			mapP.put(name, pval);
			mapL.put(name, limit);
		});

		return new Map[]{mapP, mapL};
	}

	private Map<String, Integer> count(Set<String> selectedGenes)
	{
		Map<String, Integer> cnt = new HashMap<>();

		for (String name : geneSets.keySet())
		{
			int c = CollectionUtil.countOverlap(geneSets.get(name), selectedGenes);
			cnt.put(name, c);
		}
		return cnt;
	}

	public void writeEnrichmentResults(Set<String> selected, Set<String> background, int minMemberSize, String filename)
		throws IOException
	{
		if (background != null) crop(background, minMemberSize);

		final Map<String, Double>[] pVals = getEnrichmentPvals(selected);
		Map<String, Double> qVals = FDR.getQVals(pVals[0], pVals[1]);

		List<String> names = new ArrayList<>(qVals.keySet());
		names.sort(Comparator.comparing(o -> pVals[0].get(o)));

		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));

		writer.write("# Gene set name: Name of the gene set in MSigDB.\n");
		writer.write("# URL: The URL of the gene set in MSigDB.\n");
		writer.write("# P-value: Enrichment p-value of the gene set calculated by Fisher's exact test.\n");
		writer.write("# Q-value: Estimated FDR (false discovery rate) if this p-value is used as cutoff threshold.\n");
		writer.write("# Hit size: Number of query genes that overlaps with this gene set.\n");
		writer.write("# Effective gene set size: Intersection of genes in this gene set and the given background.\n");
		writer.write("# Molecules contributed to enrichment: Names of the selected genes that overlaps with this gene set.\n");
		writer.write("Gene set name\tURL\tP-value\tQ-value\tHit size\tEffective gene set size\tMolecules contributed to enrichment");
		for (String name : names)
		{
			if (pVals[0].get(name) > 0.05) break;

			Set<String> hitGenes = CollectionUtil.getIntersection(getGeneSet(name), selected);

			int geneSetSize = getGeneSet(name).size();

			FileUtil.lnwrite(name + "\t" + getURL(name) + "\t" + pVals[0].get(name) + "\t" + qVals.get(name) +
				"\t" + hitGenes.size() + "\t" + geneSetSize + "\t" + CollectionUtil.merge(hitGenes, " "), writer);
		}
		writer.close();
	}

	interface NameFilter
	{
		boolean select(String name);
	}

	public static void main(String[] args) throws IOException
	{
//		get().printSetsNameFiltered(name -> name.contains("TARGET_GENES"));
//		writeHighlightFiles();
		printMembers();
	}

	private static void writeHighlightFiles() throws IOException
	{
		String dir = "/Users/ozgun/Documents/Temp/highlight/";
		NameFilter filter = name -> name.startsWith("HALLMARK_");
		Map<String, Set<String>> sets = get().getSetsNameFiltered(filter);
		for (String name : sets.keySet())
		{
			BufferedWriter writer = Files.newBufferedWriter(Paths.get(dir + name + ".highlight"));
			Set<String> set = sets.get(name);
			set.stream().sorted().forEach(gene -> FileUtil.writeln("node\t" + gene, writer));
			writer.close();
		}
	}

	private static void printMembers()
	{
		String name = "HALLMARK_ADIPOGENESIS";
		Set<String> geneSet = get().getGeneSet(name);
		geneSet.stream().sorted().forEach(System.out::println);
	}
}
