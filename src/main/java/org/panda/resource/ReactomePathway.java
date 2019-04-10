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
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Provides pathways and their contents in Pathway Commons. Also runs enrichment analyses.
 *
 * @author Ozgun Babur
 */
public class ReactomePathway extends FileServer
{
	private static ReactomePathway instance;

	private Map<String, Set<String>> gene2pathway;
	private Map<String, Set<String>> pathway2gene;
	private Map<String, String> pathway2name;

	public static synchronized ReactomePathway get()
	{
		if (instance == null) instance = new ReactomePathway();
		return instance;
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"Reactome-pathways.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"https://www.pathwaycommons.org/archives/PC2/v11/PathwayCommons11.reactome.hgnc.gmt.gz"};
	}

	@Override
	public boolean load() throws IOException
	{
		pathway2gene = new HashMap<>();
		gene2pathway = new HashMap<>();
		pathway2name = new HashMap<>();

		Set<Set<String>> groups = new HashSet<>();

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(line -> line.split("\t")).forEach(t ->
		{
			String id = t[0];
			String name = t[1].substring(t[1].indexOf(" ") + 1, t[1].indexOf(";"));
			pathway2name.put(id, name);

			Set<String> group = new HashSet<>(Arrays.asList(t).subList(2, t.length));
			if (groups.contains(group))
			{
//				System.err.println("Duplicate pathway: " + name);
				return;
			}
			groups.add(group);

			pathway2gene.put(id, group);

			group.stream().filter(gene -> !gene2pathway.containsKey(gene))
				.forEach(gene -> gene2pathway.put(gene, new HashSet<>()));

			group.forEach(gene -> gene2pathway.get(gene).add(id));
		});

		return true;
	}

	public Set<String> getPathways(String gene)
	{
		return gene2pathway.getOrDefault(gene, Collections.emptySet());
	}

	public Set<String> getGenes(String pathwayID)
	{
		return pathway2gene.getOrDefault(pathwayID, Collections.emptySet());
	}

	public Map<String, Set<String>> getCroppedPathways(Set<String> queryGenes)
	{
		Map<String, Set<String>> map = new HashMap<>();

		Set<Set<String>> memory = new HashSet<>();

		for (String id : pathway2gene.keySet())
		{
			Set<String> genes = pathway2gene.get(id);
			Set<String> overlap = CollectionUtil.getIntersection(genes, queryGenes);

			if (overlap.size() > 1 && !memory.contains(overlap))
			{
				map.put(id, overlap);
				memory.add(overlap);
			}
		}

		return map;
	}

	public Map<String, Set<String>> getAllPathways()
	{
		return pathway2gene;
	}

	public String getName(String id)
	{
		return pathway2name.get(id);
	}

	public String getCoverageStr(String id, Set<String> genes)
	{
		if (!pathway2gene.containsKey(id)) return null;
		Set<String> set = new HashSet<String>(genes);
		set.retainAll(pathway2gene.get(id));
		return set.size() + "/" + pathway2gene.get(id).size();
	}

	/**
	 * Gets pathways sorted to their containment of the query genes. Does not control the density,
	 * but only the count, so it is likely that first pathways will be big ones.
	 */
	public TreeMap<String, Integer> getSortedPathways(Collection<String> genes)
	{
		final Map<String, Integer> map = new HashMap<>();

		for (String gene : genes)
		{
			for (String pathway : getPathways(gene))
			{
				if (map.containsKey(pathway)) map.put(pathway, map.get(pathway) + 1);
				else map.put(pathway, 1);
			}
		}

		TreeMap<String, Integer> sorted = new TreeMap<>((o1, o2) -> map.get(o2).compareTo(map.get(o1)));

		sorted.putAll(map);
		return sorted;
	}

	public static void main(String[] args) throws IOException
	{
		Set<String> set1 = get().getGenes("http://identifiers.org/reactome/R-HSA-2219530");
		Set<String> set2 = get().getGenes("http://identifiers.org/reactome/R-HSA-1257604");
		Set<String> set3 = get().getGenes("http://identifiers.org/reactome/R-HSA-400253");

		CollectionUtil.printVennSets(set1, set2, set3);
	}
}
