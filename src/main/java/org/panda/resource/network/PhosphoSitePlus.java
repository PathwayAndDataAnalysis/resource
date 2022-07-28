package org.panda.resource.network;

import org.panda.resource.FileServer;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.FileUtil;
import org.panda.utility.graph.DirectedGraph;
import org.panda.utility.graph.SiteSpecificGraph;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Serves the PhosphoSitePlus database relations.
 *
 * @author Ozgun Babur
 */
public class PhosphoSitePlus extends FileServer
{
	static PhosphoSitePlus instance;

	Map<SignedType, SiteSpecificGraph> graphs;

	public static synchronized PhosphoSitePlus get()
	{
		if (instance == null) instance = new PhosphoSitePlus();
		return instance;
	}

	public SiteSpecificGraph getGraph(SignedType type)
	{
		return graphs.get(type);
	}

	public Map<SignedType, DirectedGraph> getAllGraphs()
	{
		return new HashMap<>(graphs);
	}

	@Override
	public String[] getLocalFilenames()
	{
		return new String[]{"psp-regulatory-sites", "psp-prot-name-to-gene-symbol.txt", "psp-Kinase_Substrate_Dataset"};
	}

	@Override
	public boolean load() throws IOException
	{
		graphs = new HashMap<>();
		graphs.put(SignedType.PHOSPHORYLATES, new SiteSpecificGraph("PhosphoSitePlus phosphorylation", SignedType.PHOSPHORYLATES.getTag()));
		graphs.put(SignedType.DEPHOSPHORYLATES, new SiteSpecificGraph("PhosphoSitePlus dephosphorylation", SignedType.DEPHOSPHORYLATES.getTag()));
		graphs.put(SignedType.ACETYLATES, new SiteSpecificGraph("PhosphoSitePlus acetylation", SignedType.ACETYLATES.getTag()));
		graphs.put(SignedType.DEACETYLATES, new SiteSpecificGraph("PhosphoSitePlus deacetylation", SignedType.DEACETYLATES.getTag()));
		graphs.put(SignedType.METHYLATES, new SiteSpecificGraph("PhosphoSitePlus methylation", SignedType.METHYLATES.getTag()));
		graphs.put(SignedType.DEMETHYLATES, new SiteSpecificGraph("PhosphoSitePlus demethylation", SignedType.DEMETHYLATES.getTag()));

		Set<String> mods = new HashSet<>(Arrays.asList("p", "ac", "m1", "m2", "me"));

		Map<String, String> protNameToSym = Files.lines(Paths.get(locateInBase(getLocalFilenames()[1])))
			.map(l -> l.split("\\||:"))
			.collect(Collectors.toMap(t -> t[2].replaceAll("_", " "), t -> t[1]));

		Files.lines(Paths.get(locateInBase(getLocalFilenames()[0]))).map(l -> l.split("\t"))
			.filter(t -> t.length > 13)
			.filter(t -> t[6].equals("human"))
			.filter(t -> !t[13].isEmpty())
			.forEach(t ->
		{
			String site = t[7];

			String mod = site.substring(site.lastIndexOf("-") + 1);
			if (!mods.contains(mod)) return;

			site = site.substring(0, site.lastIndexOf("-"));
			String target = t[0];

			for (String sourceName : t[13].split("; "))
			{
				boolean effect = true;
				if (sourceName.endsWith("(DISRUPTS)"))
				{
					effect = false;
				}
				else if (!sourceName.endsWith("(INDUCES)"))
				{
					continue;
				}

				sourceName = sourceName.substring(0, sourceName.lastIndexOf("("));

				// Skip this if there is conflicting control
				String conflict = sourceName + "(" + (effect ? "DISRUPTS" : "INDUCES") + ")";
				if (t[13].contains(conflict)) continue;

				Set<String> pmids = Arrays.stream(t[15].split("; ")).map(s -> "PMID:" + s).collect(Collectors.toSet());

				SignedType st = effect ?
					(mod.equals("p") ? SignedType.PHOSPHORYLATES : mod.equals("ac") ? SignedType.ACETYLATES : mod.equals("m1") || mod.equals("m2") || mod.equals("me") ? SignedType.METHYLATES : null) :
					(mod.equals("p") ? SignedType.DEPHOSPHORYLATES : mod.equals("ac") ? SignedType.DEACETYLATES : mod.equals("m1") || mod.equals("m2") || mod.equals("me") ? SignedType.DEMETHYLATES : null);

				if (st == null) throw new RuntimeException("Shouldn't reach here. There is a bug!");

				SiteSpecificGraph graph = graphs.get(st);

				String sGenes = protNameToSym.get(sourceName);
				if (sGenes != null)
				{
					for (String source : sGenes.split(";"))
					{
						if (!source.equals(target))
						{

//							if (target.equals("RB1") && source.equals("ABL1"))
//							{
//								System.out.println(Arrays.toString(t));
//							}


							if (!graph.hasRelation(source, target))
							{
								graph.putRelation(source, target, pmids, site);
							}
							else
							{
								graph.addSite(source, target, site);
								graph.addMediators(source, target, pmids);
							}

						}
					}
				}
			}
		});

		SiteSpecificGraph graph = graphs.get(SignedType.PHOSPHORYLATES);
		getResourceAsStream(getLocalFilenames()[2], StandardCharsets.ISO_8859_1).skip(4).map(l -> l.split("\t"))
			.filter(t -> t[3].equals("human") && t[8].equals("human"))
			.filter(t -> !t[0].isEmpty() && !t[7].isEmpty()).forEach(t ->
			{
				String source = t[0];
				String target = t[7];
				String site = t[9];

				graph.putRelation(source, target, Collections.emptySet(), site);
			});

		return true;
	}

	public static void main(String[] args)
	{
//		writeAsSif();

		printUpstream("XRCC1");
	}

	private static void printUpstream(String target)
	{
		SiteSpecificGraph graph = get().getGraph(SignedType.DEPHOSPHORYLATES);
		Set<String> upstreamSet = graph.getUpstream(target);
		for (String s : upstreamSet)
		{
			System.out.println(s + "\t" + graph.getSites(s, target));
		}
	}

	private static void writeAsSif()
	{
		BufferedWriter writer = FileUtil.newBufferedWriter("/Users/ozgun/Documents/Data/PathwayCommonsV12/PhosphoSitePlus.sif");
		for (SignedType type : SignedType.values())
		{
			SiteSpecificGraph graph = get().getGraph(type);
			if (graph != null) graph.write(writer);
		}
		FileUtil.closeWriter(writer);
	}
}
