package org.panda.resource;

import org.panda.utility.CollectionUtil;
import org.panda.utility.statistics.FishersExactTest;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Ozgun Babur
 */
public class PCPathwayHGNC extends FileServer
{
	Map<String, Set<String>> idToGenes;
	Map<Set<String>, Set<Pathway>> genesToPathwaySets;
	Set<String> allGenes;

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"pc-pathway-hgnc.txt"};
	}

	@Override
	public String[] getDistantURLs()
	{
		return new String[]{"http://www.pathwaycommons.org/archives/PC2/v9/PathwayCommons9.All.hgnc.gmt.gz"};
	}

	@Override
	public boolean load() throws IOException
	{
		idToGenes = new HashMap<>();
		genesToPathwaySets = new HashMap<>();

		getResourceAsStream(getLocalFilenames()[0]).forEach(line ->
		{
			Pathway p = new Pathway(line);
			if (p.isInvalid()) return;

			if (!genesToPathwaySets.containsKey(p.genes))
			{
				idToGenes.put(p.id, p.genes);
				genesToPathwaySets.put(p.genes, new HashSet<>());
			}

			genesToPathwaySets.get(p.genes).add(p);
		});

		allGenes = genesToPathwaySets.keySet().stream().flatMap(Collection::stream).collect(Collectors.toSet());

		return true;
	}

	public void cropToBackground(Set<String> background)
	{
		Map<String, Set<String>> idToGenes = new HashMap<>();
		Map<Set<String>, Set<Pathway>> genesToPathwaySets = new HashMap<>();

		genesToPathwaySets.values().stream().flatMap(Collection::stream).distinct().forEach(p ->
		{
			p.genes.retainAll(background);

			if (p.genes.size() < 2) return;

			if (!genesToPathwaySets.containsKey(p.genes))
			{
				idToGenes.put(p.id, p.genes);
				genesToPathwaySets.put(p.genes, new HashSet<>());
			}

			genesToPathwaySets.get(p.genes).add(p);
		});

		this.idToGenes = idToGenes;
		this.genesToPathwaySets = genesToPathwaySets;
		allGenes = genesToPathwaySets.keySet().stream().flatMap(Collection::stream).collect(Collectors.toSet());
	}

	public Map<String, Double> calculateEnrichment(Set<String> genes, int minSetSize, int maxSetSize)
	{
		int size = allGenes.size();
		return idToGenes.keySet().stream()
			.filter(id -> idToGenes.get(id).size() >= minSetSize && idToGenes.get(id).size() <= maxSetSize).collect(
				Collectors.toMap(id -> id, id ->
				{
					Set<String> geneSet = idToGenes.get(id);
					int featuredSelected = CollectionUtil.countOverlap(geneSet, genes);
					return FishersExactTest.calcEnrichmentPval(size, genes.size(), geneSet.size(), featuredSelected);
				}));
	}

	public Set<Pathway> getPathways(String id)
	{
		return genesToPathwaySets.get(idToGenes.get(id));
	}

	public Set<String> getGeneSet(String id)
	{
		return idToGenes.get(id);
	}

	public class Pathway
	{
		String name;
		String source;
		String id;
		Set<String> genes;

		public Pathway(String line)
		{
			String[] t = line.split("\t");
			if (t.length < 4) return;

			id = t[0];

			String[] tt = t[1].split(";");

			name = tt[0].substring(tt[0].indexOf(": ") + 2).trim();
			source = tt[1].substring(tt[1].indexOf(": ") + 2).trim();

			genes = new HashSet<>(Arrays.asList(t).subList(2, t.length));
		}

		boolean isInvalid()
		{
			return id == null;
		}

		public Set<String> getGenes()
		{
			return genes;
		}

		public String getName()
		{
			return name;
		}

		public String getSource()
		{
			return source;
		}
	}
}
