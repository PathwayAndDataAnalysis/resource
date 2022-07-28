package org.panda.resource;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.TermCounter;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * A curated cancer gene list from Bushman lab.
 * From http://www.bushmanlab.org/links/genelists
 *
 * @author Ozgun Babur
 */
public class GWASCatalog extends FileServer
{
	private static GWASCatalog instance;

	private Map<String, Set<String>> geneSets;

	public static synchronized GWASCatalog get()
	{
		if (instance == null) instance = new GWASCatalog();
		return instance;
	}

	public Set<String> getAllSymbols()
	{
		return geneSets.values().stream().flatMap(Set::stream).collect(Collectors.toSet());
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"GWASCatalog.tsv"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"https://www.ebi.ac.uk/gwas/api/search/downloads/full"};
	}

	@Override
	public boolean load() throws IOException
	{
		TermCounter tc = new TermCounter();
		geneSets = new HashMap<>();

		String[] header = getResourceAsStream(getLocalFilenames()[0]).findFirst().get().split("\t");

		int nameIndex = ArrayUtil.indexOf(header, "DISEASE/TRAIT");
		int geneIndex = ArrayUtil.indexOf(header, "REPORTED GENE(S)");

		getResourceAsStream(getLocalFilenames()[0]).skip(1).map(l -> l.split("\t"))
			.forEach(token -> {
				String trait = token[nameIndex];
				Set<String> genes = new HashSet<>(Arrays.asList(token[geneIndex].split(", ")));

				genes = genes.stream().filter(g -> HGNC.get().getSymbol(g) != null).collect(Collectors.toSet());

				if (!genes.isEmpty())
				{
					if (!geneSets.containsKey(trait)) geneSets.put(trait, new HashSet<>());
					geneSets.get(trait).addAll(genes);
				}
			});

		tc.print();

		return true;
	}

	public Set<String> getMembersOf(String trait)
	{
		return geneSets.get(trait);
	}

	public Map<String, Set<String>> getGeneSets()
	{
		return geneSets;
	}

	public static void main(String[] args)
	{
		System.out.println(get().getMembersOf("Subjective well-being"));

		System.out.println(get().geneSets.keySet().size());
	}
}
